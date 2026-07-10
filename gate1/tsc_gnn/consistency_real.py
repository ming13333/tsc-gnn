"""consistency_real.py — 真实数据一致性检查（held-out 基因预测）。

目的：把 Phase 2 半合成基准的"TSC-GNN 胜出"补成"真实数据也有边际增益"的更硬证据。
任务：给定非 held-out 基因表达 + 连续炎症状态，预测 held-out 基因表达。
     这是 Gate 1 held-out 任务（真实数据、真实基因、真实状态），只是把模型从
     线性/粗条件化换成 TSC-GNN。

防泄漏设计：
  - held-out 基因的节点表达特征被 mask=0；消息传递后该基因嵌入只含邻居
    （非 held-out 基因）的真实表达 → 预测 y_target 不泄漏自身值。
  - GRN 是全局基因-基因结构先验（非细胞特异），读头只用 80% 训练细胞拟合，
    测试细胞仅用于评估 → 与线性基线公平对照。

模型：
  - linear      : Ridge(y_target ~ X_input)   —— Ahlmann-Eltze 式扁平基线
  - coarse      : Ridge(y_target ~ [X_input | state])
  - tscgnn      : Ridge(y_target ~ [图嵌入_target | state])  —— 状态门控图消息传递
  - tscgnn_g0   : 同上但 gamma=0（纯图、无状态门控）→ 检验增益是否来自状态门控
  - tscgnn_perm : 同上但 GRN 边随机置换（结构无关）→ 检验增益是否来自真实图结构
"""
import numpy as np
from scipy import sparse as sp

from . import io_data, grn, model as M
from .model import graph_propagate, _standardize, _ridge_solve
from gate1 import evaluate as E


def _embed(Z, p, state, g, K, gamma, include_state):
    """返回 (n, G, d_in)：基因嵌入。偶数通道 = x 与图传播(Ax, A^2x...)；
    状态门控边已含在传播里。include_state 时追加 state 通道。"""
    H = graph_propagate(Z, p, state, g, K=K, gamma=gamma)  # (n, G, 2*(K+1))
    nb_chan = H.shape[2] // 2
    even = [H[:, :, 2 * h] for h in range(nb_chan)]         # 每 hop 的 x 通道
    Emat = np.stack(even, axis=-1)                          # (n, G, nb_chan)
    if include_state:
        st = state[:, None, None]
        Emat = np.concatenate(
            [Emat, np.broadcast_to(st, Emat.shape[:2] + (1,))], axis=-1)
    return Emat


def _permute_grn(g, rng):
    """随机置换节点顺序 → 边结构随机化但保持稀疏度/度数，作为结构无关对照。"""
    G = g["G"]
    perm = rng.permutation(G)
    A = g["A"]
    A_perm = A[perm][:, perm]
    Aaff_perm = g["A_aff"][perm][:, perm]
    sm, ssd = g["state_mean"], g["state_std"]

    def build_adj(z_i, gamma=1.0, add_self=True):
        data = A_perm.data * (1.0 + gamma * Aaff_perm.data * z_i)
        data = np.clip(data, -2.0, 2.0)
        Ai = sp.csr_matrix((data, A_perm.indices, A_perm.indptr),
                           shape=(G, G))
        if add_self:
            Ai = Ai + sp.identity(G, format="csr")
        rs = np.asarray(Ai.sum(1)).ravel()
        rs[rs == 0] = 1.0
        Ai = sp.diags(1.0 / rs) @ Ai
        return Ai

    return dict(A=A_perm, A_aff=Aaff_perm, state_mean=sm, state_std=ssd,
                build_adj=build_adj, G=G, k=g["k"])


def _ridge_predict(Ftr, Ytr_z, Fte, lam=1.0):
    """Ridge 读头（z 空间）。Ftr/Fte: (m, d)；Ytr_z: (m,1)。返回预测 (m,1) z空间。"""
    Ftr_z, pp = _standardize(Ftr)
    Fte_z, _ = _standardize(Fte, params=pp)
    W = _ridge_solve(Ftr_z, Ytr_z, lam=lam)
    return Fte_z @ W


def _ridge_predict_multi(Ftr, Ytr_mat, Fte, lam=1.0):
    """多输出 Ridge 读头（z 空间，逐细胞）。Ftr/Fte: (n, d)；Ytr_mat: (n, n_t)。
    返回预测 (n, n_t) z空间。避免把 (n*n_t, d) 摊平导致 OOM（d 大时）。"""
    Ftr_z, pp = _standardize(Ftr)
    Fte_z, _ = _standardize(Fte, params=pp)
    W = _ridge_solve(Ftr_z, Ytr_mat, lam=lam)  # (d, n_t)
    return Fte_z @ W


def run(data_root="data", cohorts=None, n_top_genes=1000, K=2, gamma=1.0,
        target_frac=0.3, n_eval=20000, n_boot=500, seed=2026):
    cohorts = cohorts or {}
    t0 = __import__("time").time()
    data = io_data.load_phase2_data(data_root, cohorts=cohorts,
                                     n_top_genes=n_top_genes)
    X, state = data["X"], data["state"]
    n, G = X.shape
    print(f"[consist] 加载 X={X.shape} time_dist="
          f"{dict(zip(*np.unique(data['time_label'], return_counts=True)))}",
          flush=True)

    rng = np.random.default_rng(seed)
    # 细胞下采样（控时）
    if n_eval and n_eval < n:
        cidx = rng.choice(n, n_eval, replace=False)
        X, state = X[cidx], state[cidx]
        n = n_eval
    # held-out 基因划分
    n_t = max(1, int(round(target_frac * G)))
    target_idx = rng.choice(G, size=n_t, replace=False)
    input_mask = np.ones(G, bool)
    input_mask[target_idx] = False
    X_input = X[:, input_mask]
    y_target = X[:, target_idx]
    print(f"[consist] held-out: n_t={n_t} 输入基因={int(input_mask.sum())} "
          f"y_target std={y_target.std():.3f}", flush=True)

    # 真实 GRN（全细胞）
    g = grn.build_grn(X, state, k=15, n_sub=4000, seed=seed)
    print(f"[consist] GRN edges={g['A'].nnz}", flush=True)

    # 节点特征：held-out 基因 mask=0
    Z = X.copy()
    Z[:, target_idx] = 0.0
    p = np.zeros((n, G), dtype=np.float32)

    # 细胞 80/20 划分
    cidx = np.arange(n)
    rng.shuffle(cidx)
    n_tr = int(0.8 * n)
    tr, te = cidx[:n_tr], cidx[n_tr:]
    n_te = len(te)

    # y 标准化（z 空间，公平对照）：先对全部扁平目标求 mu/sd，再按细胞切分重映射
    Yall = y_target.reshape(-1, 1)
    Yz, pp_y = _standardize(Yall)  # pp_y = (mu, sd)
    Yz_mat = Yz.reshape(n, n_t)            # (n, n_t) z空间
    Ytr_mat = Yz_mat[tr]                   # (n_tr, n_t)
    Yte_mat = Yz_mat[te]                   # (n_te, n_t)
    Ytr_flat = Yz_mat[tr].reshape(-1, 1)   # (n_tr*n_t, 1) 供 GNN 摊平读头

    # ---- 线性基线（逐细胞多输出 Ridge，避免大特征维 OOM）----
    yhat_lin = _ridge_predict_multi(X_input[tr], Ytr_mat, X_input[te])

    # ---- 粗条件化（X + state）----
    Ftr_c = np.concatenate([X_input[tr], state[tr][:, None]], axis=1)
    Fte_c = np.concatenate([X_input[te], state[te][:, None]], axis=1)
    yhat_coarse = _ridge_predict_multi(Ftr_c, Ytr_mat, Fte_c)

    # ---- TSC-GNN（真实 GRN，gamma=1；嵌入维度小，可摊平）----
    Eg = _embed(Z, p, state, g, K, gamma, include_state=True)
    E_tr = Eg[tr][:, target_idx, :].reshape(-1, Eg.shape[-1])
    E_te = Eg[te][:, target_idx, :].reshape(-1, Eg.shape[-1])
    yhat_gnn = _ridge_predict(E_tr, Ytr_flat, E_te).reshape(n_te, n_t)

    # ---- 消融1：gamma=0（纯图无状态门控）----
    Eg0 = _embed(Z, p, state, g, K, 0.0, include_state=True)
    E0_tr = Eg0[tr][:, target_idx, :].reshape(-1, Eg0.shape[-1])
    E0_te = Eg0[te][:, target_idx, :].reshape(-1, Eg0.shape[-1])
    yhat_gnn0 = _ridge_predict(E0_tr, Ytr_flat, E0_te).reshape(n_te, n_t)

    # ---- 消融2：permuted GRN（结构无关）----
    g_perm = _permute_grn(g, rng)
    Egp = _embed(Z, p, state, g_perm, K, gamma, include_state=True)
    Ep_tr = Egp[tr][:, target_idx, :].reshape(-1, Egp.shape[-1])
    Ep_te = Egp[te][:, target_idx, :].reshape(-1, Egp.shape[-1])
    yhat_gnnp = _ridge_predict(Ep_tr, Ytr_flat, Ep_te).reshape(n_te, n_t)

    # ---- 折叠为每细胞 MSE（z 空间，与预测同尺度）----
    def cell_mse(yhat):
        sqz = ((yhat - Yte_mat) ** 2)
        return sqz.mean(axis=1)

    cm_lin = cell_mse(yhat_lin)
    cm_coarse = cell_mse(yhat_coarse)
    cm_gnn = cell_mse(yhat_gnn)
    cm_gnn0 = cell_mse(yhat_gnn0)
    cm_gnnp = cell_mse(yhat_gnnp)

    # 全样本 MSE（z 空间均值）
    def mse_val(cm):
        return float(cm.mean())

    res = {
        "n_cells": n, "G": G, "n_target": n_t,
        "mse_linear": mse_val(cm_lin),
        "mse_coarse": mse_val(cm_coarse),
        "mse_tscgnn": mse_val(cm_gnn),
        "mse_tscgnn_g0": mse_val(cm_gnn0),
        "mse_tscgnn_perm": mse_val(cm_gnnp),
    }
    ri_v_lin = E.bootstrap_rel_improvement(cm_lin, cm_gnn, n_boot, seed)
    ri_v_coarse = E.bootstrap_rel_improvement(cm_coarse, cm_gnn, n_boot, seed)
    ri_v_g0 = E.bootstrap_rel_improvement(cm_lin, cm_gnn0, n_boot, seed)
    ri_v_perm = E.bootstrap_rel_improvement(cm_lin, cm_gnnp, n_boot, seed)
    ri_real_vs_perm = E.bootstrap_rel_improvement(cm_gnnp, cm_gnn, n_boot, seed)
    ri_g1_vs_g0 = E.bootstrap_rel_improvement(cm_gnn0, cm_gnn, n_boot, seed)
    res["ri_vs_linear"] = ri_v_lin
    res["ri_vs_coarse"] = ri_v_coarse
    res["ri_vs_g0"] = ri_v_g0
    res["ri_vs_perm"] = ri_v_perm
    res["ri_real_vs_perm"] = ri_real_vs_perm
    res["ri_g1_vs_g0"] = ri_g1_vs_g0
    res["elapsed"] = __import__("time").time() - t0
    return res


def format_report(res, title="真实数据一致性检查"):
    def f(x):
        return f"{x*100:.1f}%"
    L = []
    L.append(f"===== {title} =====")
    L.append(f"规模: n_cells={res['n_cells']} G={res['G']} held-out genes={res['n_target']}")
    L.append(f"MSE  linear      : {res['mse_linear']:.4f}")
    L.append(f"MSE  coarse      : {res['mse_coarse']:.4f}")
    L.append(f"MSE  tscgnn(g=1) : {res['mse_tscgnn']:.4f}")
    L.append(f"MSE  tscgnn(g=0) : {res['mse_tscgnn_g0']:.4f}")
    L.append(f"MSE  tscgnn(perm): {res['mse_tscgnn_perm']:.4f}")
    L.append("")
    m, lo, hi = res["ri_vs_linear"]
    L.append(f"TSC-GNN(g=1) vs LINEAR   : rel_imp={f(m)} (95% CI {f(lo)}..{f(hi)}) "
             f"{'✅显著' if lo>0 else '—不显著'}")
    m, lo, hi = res["ri_vs_coarse"]
    L.append(f"TSC-GNN(g=1) vs COARSE   : rel_imp={f(m)} (95% CI {f(lo)}..{f(hi)}) "
             f"{'✅显著' if lo>0 else '—不显著'}")
    m, lo, hi = res["ri_vs_g0"]
    L.append(f"TSC-GNN(g=0,纯图无门控) vs LINEAR: rel_imp={f(m)} "
             f"(95% CI {f(lo)}..{f(hi)}) "
             f"{'✅显著' if lo>0 else '—不显著'}")
    m, lo, hi = res["ri_g1_vs_g0"]
    L.append(f"TSC-GNN(g=1) vs g=0(纯图)        : rel_imp={f(m)} "
             f"(95% CI {f(lo)}..{f(hi)})  "
             f"({'状态门控增益' if m>0 else '纯图更优→稳态下门控非必需'})")
    m, lo, hi = res["ri_vs_perm"]
    L.append(f"随机图 vs LINEAR (通用平滑/降维) : rel_imp={f(m)} "
             f"(95% CI {f(lo)}..{f(hi)})")
    m, lo, hi = res["ri_real_vs_perm"]
    L.append(f"真实GRN vs 随机图 (结构特异性增益): rel_imp={f(m)} "
             f"(95% CI {f(lo)}..{f(hi)})  "
             f"{'✅真实结构额外贡献' if m>0 else '—无结构特异性'}")
    m, lo, hi = res["ri_g1_vs_g0"]
    L.append(f"TSC-GNN(g=1) vs g=0(纯图)        : rel_imp={f(m)} "
             f"(95% CI {f(lo)}..{f(hi)})  "
             f"({'状态门控增益' if m>0 else '纯图更优→稳态下门控非必需'})")
    L.append(f"耗时 {res['elapsed']:.1f}s")
    return "\n".join(L)
