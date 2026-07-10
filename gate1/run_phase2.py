"""run_phase2.py — TSC-GNN 深核端到端验证（Phase 2 主入口）。

流程：
  1. 加载整合卒中时间序列 → log 表达 + HVG + 连续状态。
  2. 推断状态条件化 GRN（真实数据）。
  3. 半合成扰动基准（signal 版）：ground truth 由 GRN 传播生成。
  4. 训练/评估：TSC-GNN vs 线性/粗条件化 + 两个消融。
  5. null 对照：Δ=纯噪声 → TSC-GNN 不应假阳性超越。
  6. 打印报告（同时写 gate1/phase2_report.log）。
"""
import os
import sys
import time

import numpy as np

_ROOT = os.path.dirname(os.path.abspath(__file__))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tsc_gnn import io_data, grn, synthetic_pert, train_eval

DATA_ROOT = "data"
COHORTS = {"sex": "male", "age": "W"}  # young (W8/W10) male MCAO 时间轴
N_HVG = 1000
K = 2
GAMMA = 1.0
N_PERT = 300
N_EVAL = 30000   # 评估下采样细胞数（控时；GRN 仍用全数据推断）
N_BOOT = 500
SEED = 2026


def main():
    t0 = time.time()
    print(f"[phase2] 加载数据 cohorts={COHORTS} n_hvg={N_HVG} ...", flush=True)
    data = io_data.load_phase2_data(DATA_ROOT, cohorts=COHORTS, n_top_genes=N_HVG)
    X, state, time_label = data["X"], data["state"], data["time_label"]
    n, G = X.shape
    print(f"[phase2] 数据: n_cells={n} G={G} time_dist="
          f"{dict(zip(*np.unique(time_label, return_counts=True)))}", flush=True)

    print(f"[phase2] 推断 GRN (k=15, 全数据) ...", flush=True)
    g = grn.build_grn(X, state, k=15, n_sub=4000, seed=SEED)
    print(f"[phase2] GRN 完成: edges={g['A'].nnz} G={g['G']}", flush=True)

    # 评估细胞下采样（确定性）
    rng = np.random.default_rng(SEED)
    ev_idx = rng.choice(n, size=min(N_EVAL, n), replace=False)
    Xs, ss, tls = X[ev_idx], state[ev_idx], time_label[ev_idx]

    # 扰动基因池（从 HVG 随机取）
    pert_genes = rng.choice(G, size=N_PERT, replace=False)
    bench = synthetic_pert.build_benchmark(
        Xs, ss, tls, g, pert_genes, K=K, gamma=GAMMA, mode="signal", seed=SEED)
    p, delta, time_onehot = bench["p"], bench["delta"], bench["time_onehot"]
    tr_mask, te_mask = synthetic_pert.split_by_perturbation(
        bench["pert_of_cell"], train_frac=0.8, seed=SEED)
    print(f"[phase2] 基准(signal): 训练细胞={int(tr_mask.sum())} 测试细胞={int(te_mask.sum())} "
          f"Δ std={delta.std():.3f}", flush=True)

    print(f"[phase2] 训练+评估 TSC-GNN vs baselines (bootstrap={N_BOOT}) ...", flush=True)
    res_sig, _ = train_eval.evaluate_all(
        Xs, p, ss, time_onehot, delta, tr_mask, te_mask, g,
        K=K, gamma=GAMMA, n_boot=N_BOOT, seed=SEED)
    print(train_eval.format_report(res_sig, "SIGNAL 基准（ground truth=GRN 传播）"), flush=True)

    # null 对照
    bench_null = synthetic_pert.build_benchmark(
        Xs, ss, tls, g, pert_genes, K=K, gamma=GAMMA, mode="null", seed=SEED + 1)
    p_n, delta_n = bench_null["p"], bench_null["delta"]
    tr_n, te_n = synthetic_pert.split_by_perturbation(
        bench_null["pert_of_cell"], train_frac=0.8, seed=SEED)
    print(f"[phase2] 基准(null): Δ=纯噪声 ...", flush=True)
    res_null, _ = train_eval.evaluate_all(
        Xs, p_n, ss, time_onehot, delta_n, tr_n, te_n, g,
        K=K, gamma=GAMMA, n_boot=N_BOOT, seed=SEED)
    print(train_eval.format_report(res_null, "NULL 对照（Δ 与图无关，应不超越）"), flush=True)

    print(f"[phase2] 总耗时 {time.time()-t0:.1f}s", flush=True)
    return res_sig, res_null


if __name__ == "__main__":
    main()
