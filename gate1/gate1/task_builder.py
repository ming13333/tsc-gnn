"""
task_builder.py — 把表达矩阵 + 时间标签 + 状态 score 转换成 Gate 1 的回归任务。

任务定义（在 24h→2d→14d 轴上验证"条件化是否必要"）：
    目标 y_c = 细胞 c 从当前时间点 t 到下一时间点 t+1 的"表达位移"。
    - 合成数据：y 直接给定（已知真实位移）。
    - 真实数据：y = centroid_{t+1} − x_c（下一时间点伪 bulk 质心相对当前细胞的位移）。
      这一步刻意让目标依赖时间：时间盲模型无法知道该预测哪个未来质心 → 必然失败；
      条件化模型拿到 time/state 上下文即可还原。这正是论文的 load-bearing 假设。

特征：
    - 线性基线特征 X_feat = x_c（HVG 表达）
    - 粗条件化特征  C_feat = [x_c | onehot(time) | state_score]
"""
import numpy as np


def make_context(time_idx, state, n_times):
    """构造粗条件化上下文：时间 one-hot 拼接连续状态 score。"""
    n = len(time_idx)
    oh = np.zeros((n, n_times))
    oh[np.arange(n), time_idx] = 1.0
    if state is not None and state.shape[0] == n:
        return np.hstack([oh, state]).astype(float)
    return oh.astype(float)


def build_centroid_shift(X, time_idx, n_times):
    """真实数据：计算每个非最后时间点细胞的 centroid 位移目标。"""
    X = np.asarray(X, dtype=float)
    centroids = []
    for t in range(n_times):
        mask = time_idx == t
        if mask.sum() == 0:
            centroids.append(np.zeros(X.shape[1]))
        else:
            centroids.append(X[mask].mean(axis=0))
    centroids = np.array(centroids)  # (n_times, n_genes)
    y = np.full_like(X, np.nan)
    has_target = np.zeros(X.shape[0], dtype=bool)
    for t in range(n_times - 1):
        mask = time_idx == t
        y[mask] = centroids[t + 1][None, :] - X[mask]
        has_target[mask] = True
    return y, has_target


def chronological_splits(time_idx, has_target, n_times):
    """按时间顺序做"用早预测晚"的切分：对每个测试时间点 t（t>=1 且有 target），
    训练集 = 所有更早时间点中带 target 的细胞。返回 [(train_idx, test_idx), ...]。"""
    splits = []
    trainable = np.where(has_target)[0]
    testable_times = [t for t in range(1, n_times) if (has_target & (time_idx == t)).any()]
    for t_test in testable_times:
        test_idx = np.where(has_target & (time_idx == t_test))[0]
        train_idx = np.where(has_target & (time_idx < t_test))[0]
        if len(train_idx) > 0 and len(test_idx) > 0:
            splits.append((train_idx, test_idx))
    return splits


def random_splits(time_idx, has_target, test_frac=0.2, seed=1):
    """跨时间分层随机切分：在每个时间点内部做 80/20 划分，保证训练集
    覆盖所有时间点的 onehot（有方差），测试集的时间上下文模型已见过。
    这是 Gate 1 公平性切分——验证'条件化(独占 state + 显式 time)是否优于线性'。"""
    rng = np.random.default_rng(seed)
    idx = np.where(has_target)[0]
    train, test = [], []
    for t in np.unique(time_idx[idx]):
        ti = idx[time_idx[idx] == t]
        perm = rng.permutation(len(ti))
        k = max(1, int(len(ti) * test_frac))
        test.append(ti[perm[:k]])
        train.append(ti[perm[k:]])
    if not train or not test:
        return []
    return [(np.concatenate(train), np.concatenate(test))]


def build_real_task(adata, time_key, state_keys, gene_mask=None):
    """从真实 anndata 构造任务。
    adata.X : (n_cells, n_genes) 表达（建议已 log-normalize）
    adata.obs[time_key] : 时间标签字符串，如 '24h','2d','14d'
    adata.obs[state_keys] : 连续状态 score 列
    gene_mask : bool (n_genes,) 选 HVG；为 None 则用全部基因

    ⚠️ 已弃用（审计 2026-07-09）：本函数构造的 centroid_shift 任务 y=下一质心−x_c
    中，条件化特征含 onehot(time) 即"下一质心"答案键 → 100% 虚假 PASS（leakage）。
    不得用于任何方法学结论。请改用 build_heldout_task（run_gate1_timeseries.py）。
    """
    import warnings
    warnings.warn(
        "build_real_task (centroid_shift) 已被证明存在 leakage（onehot(time) 即答案键，"
        "产生 100% 假 PASS）。请勿用于方法学结论；请改用 build_heldout_task。",
        RuntimeWarning, stacklevel=2)
    import pandas as pd
    try:
        from scipy import sparse as _sp
        X = adata.X.toarray() if _sp.issparse(adata.X) else np.asarray(adata.X)
    except Exception:
        X = np.asarray(adata.X)
    X = np.asarray(X, dtype=float)
    if gene_mask is not None:
        X = X[:, gene_mask]
    labels = list(pd.unique(adata.obs[time_key].astype(str)))
    # 按时间顺序排序标签
    order = sorted(labels, key=lambda s: _time_sort_key(s))
    label2idx = {lab: i for i, lab in enumerate(order)}
    time_idx = adata.obs[time_key].astype(str).map(label2idx).to_numpy()
    state = adata.obs[state_keys].to_numpy(dtype=float) if state_keys else None
    n_times = len(order)
    y, has_target = build_centroid_shift(X, time_idx, n_times)
    ctx = make_context(time_idx, state, n_times)
    splits = random_splits(time_idx, has_target)
    return {
        "X": X, "ctx": ctx, "y": y, "time_idx": time_idx,
        "has_target": has_target, "splits": splits,
        "time_labels": order, "n_times": n_times,
    }


def _time_sort_key(s):
    """把 '24h','2d','14d' 等映射成可排序数值（小时）。"""
    s = str(s).lower().strip()
    try:
        if "h" in s and "d" not in s:
            return float(s.replace("h", ""))
        if "d" in s:
            return float(s.replace("d", "")) * 24
        if "w" in s:
            return float(s.replace("w", "")) * 24 * 7
        return float(s)
    except Exception:
        return 1e9


def build_heldout_task(adata, state_keys, holdout_frac=0.2, seed=0):
    """方案 A：同时间点内 held-out-gene 表达预测（in-silico perturbation/imputation 范式）。

    设计要点（避免质心位移任务的 trivial leakage）：
      - 目标 y = 抽出的 held-out 基因（**非状态 marker 基因**）的表达。
      - 线性特征 X = 其余非 marker 基因表达（共表达基线）。
      - 条件化上下文 ctx = **连续状态 score**（DAM/炎症），不再含 onehot(time) 答案键。
      - 时间维度完全退出特征，仅作为评估切分（见 random_splits / loo 由调用方决定）。

    检验命题：炎症状态是否提供超出共表达的边际预测力（状态条件化是否有真实价值）。
    若 rel_imp 显著 >0 → 状态条件化有用；若 ≈0 → 粗条件化浅（呼应 Ahlmann-Eltze 现象，
    反而说明深核 TSC-GNN 才有必要）。两种都是有价值的科学结论。
    """
    try:
        from scipy import sparse as _sp
        X = adata.X.toarray() if _sp.issparse(adata.X) else np.asarray(adata.X)
    except Exception:
        X = np.asarray(adata.X)
    X = np.asarray(X, dtype=float)
    # ---- 正确性守卫（审计 2026-07-09 增加）----
    # 1) 状态列必须存在，否则下面 adata.obs[state_keys] 会 KeyError 或返回错形状
    for k in state_keys:
        assert k in adata.obs.columns, f"[build_heldout_task] 状态列缺失: {k!r}"
    var = np.asarray(adata.var_names).astype(str)
    # 状态来源 marker 基因（从 preprocessing 常量取），需排除出 input/held-out，避免 state 泄漏 target
    try:
        from .preprocessing import DAM_GENES, INFLAM_GENES
        marker_genes = list(DAM_GENES) + list(INFLAM_GENES)
    except Exception:
        marker_genes = []
    marker_set = set(str(g).lower() for g in marker_genes)
    marker_mask = np.array([v.lower() in marker_set for v in var])
    non_marker = np.where(~marker_mask)[0]
    if len(non_marker) == 0:
        non_marker = np.arange(X.shape[1])
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(non_marker))
    n_ho = max(1, int(len(non_marker) * holdout_frac))
    ho_idx = non_marker[perm[:n_ho]]
    in_idx = non_marker[perm[n_ho:]]
    X_input = X[:, in_idx]
    y = X[:, ho_idx]
    # 状态 score 作条件化上下文（连续，2 维；无 time onehot）
    state = adata.obs[state_keys].to_numpy(dtype=float)
    ctx = state
    # 2) 状态 score 不能全零：若 marker 基因未进入 HVG，add_state_scores 会静默返回
    #    全零 → 条件化退化为线性 → 产生假阴性。此处显式拦截。
    if not np.any(np.abs(ctx) > 1e-9):
        raise RuntimeError(
            "[build_heldout_task] 状态 score 全为零（marker 基因可能未进入 HVG 子集）——"
            "条件化将退化为线性，产生假阴性。请检查 prep_for_gate 的 HVG 步骤。")
    # 跨细胞随机切分（训练/测试都覆盖所有时间点，公平比较）
    n = X.shape[0]
    r2 = np.random.default_rng(seed + 1)
    p = r2.permutation(n)
    k = max(1, int(n * 0.2))
    splits = [(p[k:], p[:k])]
    return {
        "X": X_input, "ctx": ctx, "y": y,
        "has_target": np.ones(n, dtype=bool),
        "splits": splits,
        "holdout_genes": ho_idx, "input_genes": in_idx,
        "n_input_genes": len(in_idx), "n_holdout_genes": len(ho_idx),
        "time_labels": None, "n_times": 1,
    }
