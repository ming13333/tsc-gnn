"""synthetic_pert.py — 半合成扰动基准（ground truth 由真实推断 GRN 传播生成）。

设计（标准 semi-synthetic perturbation benchmark 协议，类 GEARS/ChemCPA/Nat Methods 2025）：
  - p：每个细胞 KO 一个基因（onehot）。
  - ground truth 效应 Δ = β·tanh(A_i p) + β²·tanh(A_i² p) + 噪声，
    其中 A_i 是**与 TSC-GNN 完全相同的状态条件化邻接**（保证 GNN 可解，
    扁平 baseline 因无图传播而解不出未见过的扰动）。
  - tanh 引入非线性 + 观测噪声 → 非平凡学习任务（GNN 不能 100% 完美恢复，但远优于线性）。
  - null 模式：Δ = 纯高斯噪声（与 signal 同量级，独立于图）→ 检验 GNN 不会假阳性。

训练/测试按**扰动基因**划分（held-out perturbation）：训练集只用一部分 KO 基因，
测试集用完全未见过的 KO 基因 → 公正检验泛化，是方法 paper 的标准做法。
"""
import numpy as np
from tsc_gnn.model import graph_propagate


def build_benchmark(X, state, time_label, grn, pert_genes, K=2, gamma=1.0,
                    beta=0.6, noise=0.05, seed=1, mode="signal"):
    n, G = X.shape
    rng = np.random.default_rng(seed)
    pert_of_cell = rng.choice(np.asarray(pert_genes), size=n)
    p = np.zeros((n, G), dtype=np.float32)
    p[np.arange(n), pert_of_cell] = 1.0

    H = graph_propagate(X, p, state, grn, K=K, gamma=gamma)  # (n,G,2(K+1))
    # p 通道 hop h 的索引 = 2h+1
    Ap1 = H[:, :, 2 * 1 + 1]
    Ap2 = H[:, :, 2 * 2 + 1] if K >= 2 else np.zeros_like(Ap1)
    signal = (beta * np.tanh(Ap1) + (beta ** 2) * np.tanh(Ap2)).astype(np.float32)

    if mode == "signal":
        sig_std = float(signal.std()) + 1e-8
        delta = signal + (noise * sig_std * rng.standard_normal(signal.shape)).astype(np.float32)
    else:  # null：纯噪声，同量级
        delta = ((1.0 + noise) * signal.std() *
                 rng.standard_normal(signal.shape)).astype(np.float32)

    tmap = {"24h": 0, "2d": 1, "14d": 2, "sham": 3}
    tcode = np.array([tmap.get(str(t), 3) for t in time_label])
    time_onehot = np.eye(4)[tcode].astype(np.float32)
    return dict(p=p, pert_of_cell=pert_of_cell, delta=delta,
                time_onehot=time_onehot, signal=signal,
                Ap1=Ap1, Ap2=Ap2)


def split_by_perturbation(pert_of_cell, train_frac=0.8, seed=2):
    rng = np.random.default_rng(seed)
    uniq = np.unique(pert_of_cell)
    rng.shuffle(uniq)
    n_train = int(train_frac * len(uniq))
    train_pert = set(uniq[:n_train].tolist())
    train_mask = np.array([c in train_pert for c in pert_of_cell])
    return train_mask, ~train_mask
