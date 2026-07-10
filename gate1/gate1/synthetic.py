"""
synthetic.py — 构造模拟的 24h→2d→14d 卒中时间轴数据，用于 Gate 1 冒烟测试。

目的：在不依赖真实 GEO 下载 / scanpy 的前提下，证明 Gate 1 的"线性基线 vs 粗条件化"
判定管线能正确运行、并能识别出"时间×状态条件化信号"是否存在。

═════════════════════════════════════════════════════════════════
设计要点（v3 — 确保条件化独占线性看不到的信息）：
  · X 表达 = base + time_effect[t] + noise
            → X 只编码"时间"(time_effect)，【完全不含】细胞状态信息。
  · shift(到下一时间点的位移) = signal · [ (s·W_priv)·time_gain(t)  + 0.2·Δtime ] + noise_s
            → shift 强依赖"私有状态成分"W_priv，且随时间放大(time_gain)。
              W_priv 不进入 X，仅通过 ctx(state score s) 提供给条件化模型。

  于是：
    · 线性基线(只看 X)：X 与状态 s 完全无关 → 无法还原 shift 的 state 驱动部分
      → MSE 高、corr 低 → 应失败。
    · 粗条件化(看 X + onehot(time) + state s)：s 含 W_priv 方向 → 能还原 → MSE 低
      → rel_imp 大 → 应通过。

  null 反例：shift 与 time/state 完全无关(纯随机) → 条件化也救不回 → FAIL。

  注：真实 scRNA 中状态(DAM/activated score)会部分编码进表达谱，此合成设为
      完全独占是为 Gate 1 提供清晰信号；真实数据验证时 state 信号强度会弱于
      此处，需在真实数据上重新评估阈值。
─────────────────────────────────────────────────────────────────
"""
import numpy as np


def generate_synthetic(
    n_per_time=(400, 400, 400),
    n_genes=500,
    n_state_dims=3,
    seed=0,
    signal_strength=1.0,
    noise_level=0.3,
):
    """
    返回 dict：
        X         : (N, n_genes) 细胞表达（仅编码时间，不含状态）
        time_idx  : (N,) 时间索引 0/1/2 对应 24h/2d/14d
        state     : (N, n_state_dims) 连续细胞状态 score（完整，供 ctx 使用，独占于条件化）
        shift_true: (N, n_genes) 真实的"到下一时间点"表达位移（最后时间点设为 NaN）
        has_target: (N,) bool，是否有 target（非最后时间点）
    """
    rng = np.random.default_rng(seed)
    n_times = len(n_per_time)
    N = sum(n_per_time)
    base = rng.normal(0, 1, size=(n_genes,))

    # 时间效应（进入 X，编码时间）
    time_effect = rng.normal(0, 0.8, size=(n_times, n_genes))
    # 状态的私有强效应（不进入 X，独占于条件化 ctx，驱动 shift）
    W_priv = rng.normal(0, 1.0, size=(n_state_dims, n_genes))

    X = np.zeros((N, n_genes))
    time_idx = np.zeros(N, dtype=int)
    state = np.zeros((N, n_state_dims))
    shift_true = np.full((N, n_genes), np.nan)
    has_target = np.zeros(N, dtype=bool)

    start = 0
    for t in range(n_times):
        n = n_per_time[t]
        idx = np.arange(start, start + n)
        # 连续状态 score（与 time 部分相关但不完全决定）
        s = rng.normal(0, 1, size=(n, n_state_dims))
        # X：仅编码时间 + 噪声（不含状态）
        expr = (base[None, :]
                + time_effect[t][None, :]
                + rng.normal(0, noise_level, size=(n, n_genes)))
        X[idx] = expr
        time_idx[idx] = t
        state[idx] = s
        if t < n_times - 1:
            # shift 强依赖私有状态成分，并随 time 放大（时间×状态交互）
            time_gain = 0.6 + 0.7 * t / (n_times - 1)
            shift = signal_strength * (
                (s @ W_priv) * time_gain
                + 0.2 * (time_effect[t][None, :] - time_effect[t + 1][None, :])
            ) + rng.normal(0, noise_level * 0.4, size=(n, n_genes))
            shift_true[idx] = shift
            has_target[idx] = True
        start += n

    return {
        "X": X,
        "time_idx": time_idx,
        "state": state,
        "shift_true": shift_true,
        "has_target": has_target,
        "time_labels": ["24h", "2d", "14d"][:n_times],
    }


def generate_null(
    n_per_time=(400, 400, 400),
    n_genes=500,
    n_state_dims=3,
    seed=1,
    noise_level=0.3,
):
    """反例：shift 与 time/state 无关，仅依赖纯随机噪声。
    用于验证管线在'无信号'时不会错误地宣告通过。"""
    rng = np.random.default_rng(seed)
    n_times = len(n_per_time)
    N = sum(n_per_time)
    base = rng.normal(0, 1, size=(n_genes,))
    time_effect = rng.normal(0, 0.8, size=(n_times, n_genes))

    X = np.zeros((N, n_genes))
    time_idx = np.zeros(N, dtype=int)
    state = np.zeros((N, n_state_dims))
    shift_true = np.full((N, n_genes), np.nan)
    has_target = np.zeros(N, dtype=bool)

    start = 0
    for t in range(n_times):
        n = n_per_time[t]
        idx = np.arange(start, start + n)
        s = rng.normal(0, 1, size=(n, n_state_dims))
        X[idx] = (base[None, :]
                  + time_effect[t][None, :]
                  + rng.normal(0, noise_level, size=(n, n_genes)))
        time_idx[idx] = t
        state[idx] = s
        if t < n_times - 1:
            # 纯随机 shift，与时间/状态无关
            shift_true[idx] = rng.normal(0, noise_level, size=(n, n_genes))
            has_target[idx] = True
        start += n
    return {
        "X": X, "time_idx": time_idx, "state": state,
        "shift_true": shift_true, "has_target": has_target,
        "time_labels": ["24h", "2d", "14d"][:n_times],
    }
