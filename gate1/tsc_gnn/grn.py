"""grn.py — 从真实卒中 scRNA 推断状态条件化 GRN（基因调控图）。

图结构：对称 k-NN 相关图（每基因出边齐全 → 可对**任意**基因做 KO 并传播）。
边权重：Pearson 相关（per-study pooled 子样本估计）。
状态门控：每条边带一个 state-affinity a_e —— 该边活性随炎症状态增强/抑制的程度
          （由 (X_u*X_v) 对 state 的回归斜率估计）。细胞 i 的有效边权重 =
          w_e * (1 + gamma * a_e * z_i)，z_i = (state_i - mean)/std。

所有估计都用子样本，内存安全。
"""
import numpy as np
from scipy import sparse as sp


def _gene_corr_knn(X_sub, k=15):
    """X_sub: (n_sub, G) → 对称 k-NN 相关邻接 (G,G) 稀疏 csr。"""
    n_sub, G = X_sub.shape
    Xz = (X_sub - X_sub.mean(0)) / (X_sub.std(0) + 1e-8)
    C = (Xz.T @ Xz) / n_sub  # (G, G) 相关矩阵
    rows, cols, vals = [], [], []
    for g in range(G):
        row = C[g]
        order = np.argsort(-np.abs(row))
        nbrs = order[1:k + 1]  # 排除自身
        for nb in nbrs:
            rows.append(g)
            cols.append(int(nb))
            vals.append(float(row[nb]))
    A = sp.csr_matrix((vals, (rows, cols)), shape=(G, G))
    # 对称化：让每基因同时拥有出/入边
    A = (A + A.T) / 2.0
    A.data = np.clip(A.data, -1.0, 1.0)
    return A


def _edge_state_affinity(X_sub, A, state_sub):
    """每条边 (u,v) 的 state-affinity：cov(X_u*X_v, state)/var(state)。"""
    G = X_sub.shape[1]
    # 取 A 的上三角边（避免重复），对称图只需算一次
    coo = A.tocoo()
    edges = sorted(set((min(u, v), max(u, v)) for u, v in zip(coo.row, coo.col)))
    if not edges:
        return sp.csr_matrix((G, G))
    u = np.array([e[0] for e in edges])
    v = np.array([e[1] for e in edges])
    Xprod = X_sub[:, u] * X_sub[:, v]  # (n_sub, E)
    var_s = state_sub.var() + 1e-12
    a = (Xprod.T @ state_sub) / X_sub.shape[0]
    a = a - Xprod.mean(0) * state_sub.mean()
    a = a / var_s  # (E,)
    rows = np.concatenate([u, v])
    cols = np.concatenate([v, u])
    vals = np.concatenate([a, a])
    Aaff = sp.csr_matrix((vals, (rows, cols)), shape=(G, G))
    return Aaff


def build_grn(X, state, k=15, n_sub=4000, seed=0):
    """推断状态条件化 GRN。

    返回 dict：
        A      : (G,G) 稀疏，基础边权重（对称）
        A_aff  : (G,G) 稀疏，每条边 state-affinity
        state_mean, state_std
        build_adj(z_i, gamma) : 返回细胞 i 的（行归一化、含自环）状态条件化邻接
    """
    rng = np.random.default_rng(seed)
    n, G = X.shape
    idx = rng.choice(n, size=min(n_sub, n), replace=False)
    X_sub = X[idx]
    state_sub = state[idx]
    print(f"  [grn] 子样本 {X_sub.shape[0]} 细胞估计 k-NN(k={k}) 相关图 ...")
    A = _gene_corr_knn(X_sub, k=k)
    print(f"  [grn] 计算 {A.nnz} 条边的状态亲和度 ...")
    A_aff = _edge_state_affinity(X_sub, A, state_sub)
    s_mean, s_std = float(state.mean()), float(state.std() + 1e-8)

    def build_adj(z_i, gamma=1.0, add_self=True):
        data = A.data * (1.0 + gamma * A_aff.data * z_i)
        data = np.clip(data, -2.0, 2.0)
        Ai = sp.csr_matrix((data, A.indices, A.indptr), shape=A.shape)
        if add_self:
            Ai = Ai + sp.identity(G, format="csr")
        rs = np.asarray(Ai.sum(1)).ravel()
        rs[rs == 0] = 1.0
        Ai = sp.diags(1.0 / rs) @ Ai
        return Ai

    return dict(A=A, A_aff=A_aff, state_mean=s_mean, state_std=s_std,
                build_adj=build_adj, G=G, k=k)


# ───────────────────────────────────────────────────────────────────────────
# DoRothEA 因果方向 GRN（TF -> target 有向，带符号置信权重）
# 返回与 build_grn 完全兼容的 dict，可直接喂 model.graph_propagate。
# ───────────────────────────────────────────────────────────────────────────
import os as _os
_DOROTHEA_DIR = _os.path.join(_os.path.dirname(__file__), "..", "data", "dorothea")


def _load_dorothea(species, confidence_levels, data_dir=None):
    if data_dir is None:
        data_dir = _DOROTHEA_DIR
    path = _os.path.join(data_dir, f"{species}_dorothea_regulon.tsv")
    if not _os.path.exists(path):
        raise FileNotFoundError(
            f"DoRothEA 文件缺失: {path}。请先运行导出步骤（decoupler.op.dorothea）。")
    import pandas as pd
    df = pd.read_csv(path, sep="\t")
    df = df[df["confidence"].isin(list(confidence_levels))]
    return df


def build_dorothea_grn(genes, species="mouse", confidence_levels=("A", "B", "C"),
                       id_to_symbol=None, X=None, state=None, n_sub=4000,
                       seed=0, data_dir=None, add_self=True):
    """构建 DoRothEA 因果方向 GRN（TF -> target 有向，边权含方向/置信）。

    返回 dict（与 build_grn 同接口）：
        A           : (G,G) 有向 csr，边权 = DoRothEA weight（已 clip [-1,1]）
        A_aff       : (G,G) 有向 csr，每条边 state-affinity（仅当 X,state 提供）
        state_mean, state_std
        build_adj(z_i, gamma) : 行归一化 + 自环 的状态条件化邻接
        is_directed : True   n_edges, species

    genes         : scRNA var 名列表（按数据顺序）
    id_to_symbol  : 可选 gene_id->symbol 映射；None 则 genes 本身当作 symbol
    X, state      : 可选，用于估计有向 state-affinity（与 k-NN 同协议）
    """
    df = _load_dorothea(species, confidence_levels, data_dir=data_dir)
    G = len(genes)
    # 每基因索引对应的 symbol（统一 upper，跨物种/大小写鲁棒）
    if id_to_symbol is not None:
        sym = [str(id_to_symbol.get(g, g)).upper() for g in genes]
    else:
        sym = [str(g).upper() for g in genes]
    sym2idx = {s: i for i, s in enumerate(sym) if s not in ("", "NAN", "NONE")}

    src = df["source"].astype(str).str.upper().values
    tgt = df["target"].astype(str).str.upper().values
    w = df["weight"].astype(float).values
    # 过滤两端都在 genes 内的边（向量化）
    keys = np.array(list(sym2idx.keys()))
    src_in = np.isin(src, keys)
    tgt_in = np.isin(tgt, keys)
    m = src_in & tgt_in
    su = np.array([sym2idx[s] for s in src[m]])
    tv = np.array([sym2idx[t] for t in tgt[m]])
    ww = w[m]
    A = sp.csr_matrix((ww, (su, tv)), shape=(G, G))
    A.data = np.clip(A.data, -1.0, 1.0)
    print(f"  [dorothea] {species}: 基因覆盖 {len(sym2idx)}/{G}；"
          f"保留边 {A.nnz}（confidence={confidence_levels}）")

    # 可选：有向 state-affinity（边活性 X_TF*X_target 对 state 回归）
    A_aff = sp.csr_matrix((G, G))
    s_mean, s_std = 0.0, 1.0
    if X is not None and state is not None:
        rng = np.random.default_rng(seed)
        n = X.shape[0]
        idx = rng.choice(n, size=min(n_sub, n), replace=False)
        X_sub = np.asarray(X[idx], dtype=np.float64)
        state_sub = np.asarray(state[idx], dtype=np.float64)
        coo = A.tocoo()
        s_idx = coo.row.astype(int)
        t_idx = coo.col.astype(int)
        Xprod = X_sub[:, s_idx] * X_sub[:, t_idx]          # (n_sub, E)
        var_s = state_sub.var() + 1e-12
        a = (Xprod.T @ state_sub) / X_sub.shape[0]
        a = a - Xprod.mean(0) * state_sub.mean()
        a = a / var_s
        A_aff = sp.csr_matrix((a, (s_idx, t_idx)), shape=(G, G))
        s_mean, s_std = float(state.mean()), float(state.std() + 1e-8)
        print(f"  [dorothea] 估计 {A_aff.nnz} 条边的状态亲和度")

    def edge_weights(z_i, gamma=1.0):
        """门控后的原始有向边权（无自环、未行归一化）——rewiring 解释用。"""
        if A_aff.nnz > 0:
            data = A.data * (1.0 + gamma * A_aff.data * z_i)
        else:
            data = A.data.copy()
        data = np.clip(data, -2.0, 2.0)
        return sp.csr_matrix((data, A.indices, A.indptr), shape=A.shape)

    def build_adj(z_i, gamma=1.0):
        data = edge_weights(z_i, gamma).data.copy()
        Ai = sp.csr_matrix((data, A.indices, A.indptr), shape=A.shape)
        if add_self:
            Ai = Ai + sp.identity(G, format="csr")
        rs = np.asarray(Ai.sum(1)).ravel()
        rs[rs == 0] = 1.0
        Ai = sp.diags(1.0 / rs) @ Ai
        return Ai

    return dict(A=A, A_aff=A_aff, state_mean=s_mean, state_std=s_std,
                build_adj=build_adj, edge_weights=edge_weights, G=G,
                species=species, n_edges=A.nnz, is_directed=True,
                confidence_levels=confidence_levels)
