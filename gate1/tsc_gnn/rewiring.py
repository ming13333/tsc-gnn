"""rewiring.py v2 — 时间×状态 GRN 重布线 (edge rewiring) 可解释性模块。

v2 修复（Bug Audit 2026-07-09 阻断级 5 条）：
1. p 分辨率：n_perm=200 → p_min=1/201≈0.005
2. 多重校正：BH-FDR + 置换 pooled FDR（不依赖 p 分辨率）
3. 细胞类型组成混淆：PC 回归（n_pc=10）残差化 Xu/Xv
4. 方向标签：改为"关联增强/关联减弱"（coupling change），DoRothEA 先验方向独立报告
5. 运行清单：rewiring_table 返回 manifest dict

性能优化（解决内存抖动 + CPU 空耗）：
- 预计算 Xu/Xv/Xuv/Xu2/Xv2 一次（subsample 后 ~2GB，无 swap）
- 置换内仅做 sgemv（m @ M），无 fancy indexing / 拷贝
- SVD 仅算前 n_pc 列（scipy.linalg.svd full_matrices=False + 截断）
- 存储 null ΔW (n_perm×E) 供 pooled FDR（<40MB）
"""
import sys
import numpy as np
from scipy import sparse as sp
from scipy import linalg as spla


# ──────────────────────── 内部工具 ────────────────────────

def _residualize_pcs(Xu, Xv, X_full, n_pc=10):
    """PC 回归残差化（去除细胞类型组成效应）。in-place 修改 Xu, Xv。

    流程：
    1. 中心化 X_full → SVD → 取前 n_pc 个 PC 分数 (n, n_pc)
    2. 将 Xu, Xv 对 [1, PCs] 做最小二乘回归，减去拟合值 → 残差
    """
    n = X_full.shape[0]
    Xc = np.ascontiguousarray(
        X_full - X_full.mean(0, keepdims=True), dtype=np.float32)
    # SVD: Xc = U S V^T；PC 分数 = U[:, :n_pc] * S[:n_pc]
    U_svd, S_svd, _ = spla.svd(
        Xc, full_matrices=False, overwrite_a=True, check_finite=False)
    PCs = (U_svd[:, :n_pc] * S_svd[:n_pc]).astype(np.float32)  # (n, n_pc)
    del U_svd, S_svd  # 立即释放 ~400MB
    # 回归残差：beta = (W^T W)^{-1} W^T M；M_res = M - W @ beta
    W = np.column_stack([np.ones(n, dtype=np.float32), PCs])  # (n, n_pc+1)
    WtW_inv = np.linalg.inv(W.T @ W)  # (n_pc+1, n_pc+1)
    for M in (Xu, Xv):
        beta = WtW_inv @ (W.T @ M)  # (n_pc+1, E)
        M -= W @ beta  # in-place 残差化
    return Xu, Xv


def _coupling_sums(Xu, Xv, Xuv, Xu2, Xv2, mask):
    """用 sgemv（m @ M）计算组内 Pearson r，无行拷贝。

    Pearson r = (n*sxy - sx*sy) / sqrt((n*sx2-sx^2)*(n*sy2-sy^2))
    """
    n_t = int(mask.sum())
    E = Xu.shape[1]
    if n_t < 5:
        return np.zeros(E, dtype=np.float32)
    m = mask.astype(np.float32)  # (n,)  ~32KB
    sx = m @ Xu    # (E,) sgemv
    sy = m @ Xv
    sxy = m @ Xuv
    sx2 = m @ Xu2
    sy2 = m @ Xv2
    cov = n_t * sxy - sx * sy
    vx = np.maximum(n_t * sx2 - sx * sx, 0.0)
    vy = np.maximum(n_t * sy2 - sy * sy, 0.0)
    denom = np.sqrt(vx * vy)
    denom[denom < 1e-8] = 1.0
    r = cov / denom
    return np.clip(r, -1.0, 1.0).astype(np.float32)


def _timepoint_coupling(Xu, Xv, Xuv, Xu2, Xv2, time_label, time_points):
    """对每个时间点 t 估计耦合 r_t（组内 Pearson）。返回 {t: (E,)}。"""
    lab = np.asarray(time_label)
    return {t: _coupling_sums(Xu, Xv, Xuv, Xu2, Xv2, (lab == t))
            for t in time_points}


def _rewiring_deltas(coupling, transitions):
    """{tr: coupling[tr[1]] - coupling[tr[0]]}。"""
    return {tr: coupling[tr[1]] - coupling[tr[0]] for tr in transitions}


# ──────────────────────── 置换检验 ────────────────────────

def permutation_test(Xu, Xv, Xuv, Xu2, Xv2, time_label, time_points,
                     transitions, n_perm=200, seed=2):
    """时间标签置换 → 逐边 p 值 + null ΔW。

    预计算矩阵 Xu/Xv/Xuv/Xu2/Xv2 常驻内存，置换内仅做 sgemv。
    p_e = (1 + #{perm: |dW_perm_e| >= |dW_obs_e|}) / (n_perm + 1)

    返回 (obs_deltas, pvals, null_deltas)：
      obs_deltas  : {tr: (E,)}  观测 ΔW
      pvals       : {tr: (E,)}  逐边 p 值
      null_deltas : {tr: (n_perm, E)}  全部 null ΔW（供 pooled FDR）
    """
    rng = np.random.default_rng(seed)
    lab = np.asarray(time_label)
    n, E = Xu.shape

    # 观测
    coupling_obs = _timepoint_coupling(
        Xu, Xv, Xuv, Xu2, Xv2, lab, time_points)
    obs = _rewiring_deltas(coupling_obs, transitions)

    null_ge = {tr: np.zeros(E, dtype=np.int64) for tr in transitions}
    null_dw = {tr: np.zeros((n_perm, E), dtype=np.float32) for tr in transitions}
    base = lab.copy()

    for p in range(n_perm):
        perm = rng.permutation(base)
        coupling_perm = _timepoint_coupling(
            Xu, Xv, Xuv, Xu2, Xv2, perm, time_points)
        for tr in transitions:
            d = coupling_perm[tr[1]] - coupling_perm[tr[0]]
            null_dw[tr][p] = d
            null_ge[tr] += (np.abs(d) >= np.abs(obs[tr])).astype(np.int64)
        if (p + 1) % 50 == 0:
            print(f"    [perm] {p+1}/{n_perm}", flush=True)

    pvals = {tr: (null_ge[tr] + 1) / (n_perm + 1) for tr in transitions}
    return obs, pvals, null_dw


# ──────────────────────── 多重校正 ────────────────────────

def benjamini_hochberg(p):
    """BH-FDR 校正（逐边 p → q）。"""
    p = np.asarray(p, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    bh = ranked * n / np.arange(1, n + 1)
    bh = np.minimum.accumulate(bh[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(bh, 0.0, 1.0)
    return out


def pooled_fdr(obs_dw, null_dw, transitions):
    """置换 pooled FDR（不依赖 p 分辨率）。

    1. 逐边标准化：z = |dW| / sd(null_dW)
    2. 池化 null z（n_perm × E → 一维）
    3. 对排序后的观测 z，FDR(k) = #{null z >= z_obs(k)} / (n_perm * k)
    4. 单调性校正 + clip

    返回 {tr: (E,) q-values}。
    """
    out = {}
    for tr in transitions:
        obs = np.abs(obs_dw[tr])           # (E,)
        null = null_dw[tr]                 # (n_perm, E)
        n_perm_val, E = null.shape

        # 逐边标准化
        null_sd = null.std(0)              # (E,)
        null_sd[null_sd < 1e-8] = 1.0
        z_obs = obs / null_sd              # (E,)
        z_null = np.abs(null) / null_sd[np.newaxis, :]  # (n_perm, E)
        z_null_flat = z_null.ravel()       # (n_perm * E,)

        # 排序
        null_asc = np.sort(z_null_flat)    # 升序
        obs_order = np.argsort(-z_obs)     # 降序
        z_sorted = z_obs[obs_order]        # 降序

        # 对每个 rank k (1-indexed)，count null >= z_sorted[k-1]
        cnt = len(z_null_flat) - np.searchsorted(
            null_asc, z_sorted, side='left')
        ranks = np.arange(1, E + 1)
        fdr = cnt / (n_perm_val * ranks)

        # 单调性校正（从末尾向前取最小）
        fdr = np.minimum.accumulate(fdr[::-1])[::-1]
        fdr = np.clip(fdr, 0, 1)

        q = np.ones(E)
        q[obs_order] = fdr
        out[tr] = q
    return out


# ──────────────────────── 端到端入口 ────────────────────────

def rewiring_table(grn, genes, X, time_label, state=None,
                   transitions=None, n_perm=200, seed=2,
                   a_aff_threshold=None, n_pc=10):
    """端到端 rewiring 分析 v2。

    参数
    ----
    grn   : build_dorothea_grn 返回的 dict（含 A, A_aff, edge_weights）。
    genes : (G,) 基因名（与 X 列对齐）。
    X     : (n_cells, G) 表达（已 subsample、已归一化）。
    state : (n_cells,) 可选；用于按 |A_aff| 筛选 state-conditioned 边。
    transitions : list[(t1,t2)]；默认连续时间轴 + 总体。
    n_perm : 置换次数（200 → p_min≈0.005）。
    seed   : 随机种子。
    a_aff_threshold : 筛选 |A_aff| 超过该分位的边（None=全部）。
    n_pc   : PC 回归维数（0=不做组成校正，10=默认）。

    返回 (df, transitions, time_points, manifest)。
    """
    import pandas as pd
    import scipy

    sym = [str(g) for g in genes]
    coo = grn["A"].tocoo()
    u_idx = coo.row.astype(int)
    v_idx = coo.col.astype(int)
    E = len(u_idx)
    base_w = np.asarray(grn["A"].data, dtype=float)
    if grn.get("A_aff") is not None and grn["A_aff"].nnz > 0:
        a_aff = np.asarray(grn["A_aff"].tocoo().data, dtype=float)
    else:
        a_aff = np.zeros(E)

    # state-conditioned 边筛选
    if a_aff_threshold is not None and grn.get("A_aff") is not None:
        if a_aff_threshold < 1.0:
            thr = np.quantile(np.abs(a_aff), a_aff_threshold)
        else:
            thr = float(a_aff_threshold)
        keep = np.abs(a_aff) >= thr
        print(f"  [rewire] state-conditioned 筛选: |A_aff|>={thr:.4f} -> "
              f"保留 {int(keep.sum())}/{E} 边")
        u_idx = u_idx[keep]
        v_idx = v_idx[keep]
        base_w = base_w[keep]
        a_aff = a_aff[keep]
        E = len(u_idx)

    time_points = sorted(set(np.asarray(time_label)))
    if transitions is None:
        tp = ["sham", "24h", "2d", "14d"]
        tp = [t for t in tp if t in time_points] + \
             [t for t in time_points if t not in tp]
        transitions = [(tp[i], tp[i + 1]) for i in range(len(tp) - 1)]
        if "sham" in tp and tp[-1] != "sham":
            transitions.append(("sham", tp[-1]))

    # ── 输入校验 ──
    assert not np.isnan(X).any(), "X contains NaN"
    assert len(u_idx) == len(v_idx), "edge index mismatch"
    n = X.shape[0]
    print(f"  [rewire] n_cells={n}, n_edges={E}, n_perm={n_perm}, "
          f"n_pc={n_pc}, p_min={1/(n_perm+1):.4f}")

    # ── 预计算 Xu, Xv（未残差化）──
    Xu = np.ascontiguousarray(np.asarray(X[:, u_idx], dtype=np.float32))
    Xv = np.ascontiguousarray(np.asarray(X[:, v_idx], dtype=np.float32))

    # 未残差化耦合（对比用）
    Xuv_raw = np.ascontiguousarray(Xu * Xv)
    Xu2_raw = np.ascontiguousarray(Xu * Xu)
    Xv2_raw = np.ascontiguousarray(Xv * Xv)
    coupling_raw = _timepoint_coupling(
        Xu, Xv, Xuv_raw, Xu2_raw, Xv2_raw, time_label, time_points)
    obs_raw = _rewiring_deltas(coupling_raw, transitions)
    del Xuv_raw, Xu2_raw, Xv2_raw  # 释放 ~1.1GB

    # ── PC 残差化（组成校正）──
    if n_pc > 0:
        print(f"  [rewire] PC 回归组成校正 (n_pc={n_pc}) ...")
        Xu, Xv = _residualize_pcs(
            Xu, Xv,
            np.ascontiguousarray(np.asarray(X, dtype=np.float32)),
            n_pc=n_pc)

    # ── 预计算残差化后的乘积矩阵 ──
    Xuv = np.ascontiguousarray(Xu * Xv)
    Xu2 = np.ascontiguousarray(Xu * Xu)
    Xv2 = np.ascontiguousarray(Xv * Xv)
    mem_gb = (Xu.nbytes + Xv.nbytes + Xuv.nbytes +
              Xu2.nbytes + Xv2.nbytes) / 1e9
    print(f"  [rewire] 常驻矩阵: {mem_gb:.2f}GB "
          f"(Xu+Xv+Xuv+Xu2+Xv2)")

    # ── 置换检验 ──
    print(f"  [rewire] 置换检验 n_perm={n_perm} ...")
    obs, pvals, null_dw = permutation_test(
        Xu, Xv, Xuv, Xu2, Xv2, time_label, time_points, transitions,
        n_perm=n_perm, seed=seed)

    # ── 多重校正 ──
    fdr_bh = {tr: benjamini_hochberg(pvals[tr]) for tr in transitions}
    pooled_q = pooled_fdr(obs, null_dw, transitions)

    # 残差化观测耦合
    coupling_obs = _timepoint_coupling(
        Xu, Xv, Xuv, Xu2, Xv2, time_label, time_points)

    # ── 组装 DataFrame ──
    data = {
        "tf": [sym[int(u_idx[e])] for e in range(E)],
        "target": [sym[int(v_idx[e])] for e in range(E)],
        "dorothea_weight": base_w.astype(float),
        "dorothea_sign": [
            "activation" if w > 0 else "repression" for w in base_w],
        "a_aff": a_aff.astype(float),
    }
    for t in time_points:
        data[f"r_{t}"] = coupling_obs[t].astype(float)      # PC-corrected
        data[f"r_raw_{t}"] = coupling_raw[t].astype(float)   # uncorrected
    for tr in transitions:
        t1, t2 = tr
        data[f"dW_{t1}_{t2}"] = obs[tr].astype(float)        # PC-corrected
        data[f"dW_raw_{t1}_{t2}"] = obs_raw[tr].astype(float)
        data[f"p_{t1}_{t2}"] = pvals[tr].astype(float)
        data[f"fdr_bh_{t1}_{t2}"] = fdr_bh[tr]
        data[f"q_pooled_{t1}_{t2}"] = pooled_q[tr]
    df = pd.DataFrame(data)

    # ── Manifest ──
    manifest = {
        "n_cells": int(n),
        "n_edges_tested": int(E),
        "n_perm": int(n_perm),
        "p_min": float(1 / (n_perm + 1)),
        "n_pc": int(n_pc),
        "seed_perm": int(seed),
        "a_aff_threshold": (float(a_aff_threshold)
                            if a_aff_threshold is not None else None),
        "transitions": [list(tr) for tr in transitions],
        "time_points": list(time_points),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "pandas": pd.__version__,
        "resident_matrices_GB": round(mem_gb, 2),
    }
    return df, transitions, time_points, manifest
