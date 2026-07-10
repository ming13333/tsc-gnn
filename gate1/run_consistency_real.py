"""run_consistency_real.py — 真实数据一致性检查编排。

把 Phase 2 的 TSC-GNN 放在 Gate 1 held-out 任务（真实卒中数据）上，检验：
  TSC-GNN（状态条件化图消息传递）是否相比线性/Ridge 基线有真实边际增益。
"""
import time
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tsc_gnn import consistency_real as C

DATA_ROOT = "data"
COHORTS = {"sex": "male", "age": "W"}   # young male MCAO 时间轴
N_HVG = 1000
K = 2
GAMMA = 1.0
TARGET_FRAC = 0.3
N_EVAL = 20000
N_BOOT = 500
SEED = 2026


def main():
    t0 = time.time()
    print(f"[consist] 真实数据一致性检查 cohorts={COHORTS} "
          f"n_hvg={N_HVG} K={K} gamma={GAMMA} target_frac={TARGET_FRAC} "
          f"n_eval={N_EVAL}", flush=True)
    res = C.run(DATA_ROOT, COHORTS, N_HVG, K, GAMMA, TARGET_FRAC,
                N_EVAL, N_BOOT, SEED)
    print(C.format_report(res), flush=True)
    print(f"[consist] 总耗时 {time.time()-t0:.1f}s", flush=True)
    return res


if __name__ == "__main__":
    main()
