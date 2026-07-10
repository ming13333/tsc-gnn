"""
main.py — Gate 1 驱动脚本。

用法：
  # 合成数据冒烟测试（有信号，应 PASS）—— 验证管线逻辑正确
  python main.py --mode synthetic --signal

  # 合成反例（无信号，应 FAIL）—— 验证管线不会假阳性
  python main.py --mode synthetic --null

  # 真实数据（需先下载 + scanpy 预处理，见 preprocessing.py）
  python main.py --mode real --dataset gse174574

Gate 1 判定（见 evaluate.evaluate_gate）：
  条件化模型相对线性基线的 MSE 改善 ≥10% 且 bootstrap 95%CI 下界>0 且
  条件化每基因相关高于线性 → PASS（值得推进 TSC-GNN 深核）。
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gate1 import synthetic, task_builder, baselines, evaluate


def run_synthetic(mode, n_boot=1000):
    if mode == "signal":
        data = synthetic.generate_synthetic(seed=0, signal_strength=1.0)
        title = "Gate 1 [SYNTHETIC, 含 time×state 信号]"
    else:
        data = synthetic.generate_null(seed=1)
        title = "Gate 1 [SYNTHETIC, 无信号 / 反例]"

    X = data["X"]
    time_idx = data["time_idx"]
    state = data["state"]
    y = data["shift_true"]
    has_target = data["has_target"]
    n_times = len(data["time_labels"])

    ctx = task_builder.make_context(time_idx, state, n_times)

    # 跨时间分层随机切分（训练需见过所有时间点的 onehot 才有方差）
    splits = task_builder.random_splits(time_idx, has_target)
    if not splits:
        print("ERROR: 无可用切分（需要 ≥2 个时间点且各点有 target）")
        return

    pred_lin_all = np_full = y.copy() * 0.0  # placeholder, fill per-split
    pred_cond_all = y.copy() * 0.0
    true_all = y.copy()
    mask_all = has_target.copy()

    import numpy as np
    pred_lin_all = np.full_like(y, np.nan)
    pred_cond_all = np.full_like(y, np.nan)
    predicted = np.zeros(y.shape[0], dtype=bool)

    for tr, te in splits:
        # 线性基线：只用表达
        pl = baselines.fit_predict_linear(X[tr], y[tr], X[te])
        # 粗条件化：表达 + 时间/状态上下文
        feat_tr = np.hstack([X[tr], ctx[tr]])
        feat_te = np.hstack([X[te], ctx[te]])
        pc = baselines.fit_predict_conditional(feat_tr, y[tr], feat_te)
        pred_lin_all[te] = pl
        pred_cond_all[te] = pc
        predicted[te] = True

    res = evaluate.evaluate_gate(pred_lin_all, pred_cond_all, true_all, predicted, n_boot=n_boot)
    print(evaluate.format_report(res, title))
    return res


def run_real(dataset, n_boot=1000):
    import warnings
    warnings.warn(
        "main.py --mode real 原先走 build_real_task(质心位移)，已被证明 leakage 弃用；"
        "现自动改用正确的 held-out-gene 任务（与 run_gate1_timeseries.py 一致）。",
        RuntimeWarning, stacklevel=2)
    from gate1 import preprocessing
    adata = preprocessing.load_processed(dataset)
    # 只用高变基因(HVG)做特征，避免全基因稠密化爆内存，且聚焦信号基因
    hvg = adata.var.get("highly_variable", None)
    gene_mask = hvg.values if hvg is not None else None
    # 弃用 centroid_shift；改用方案 A（held-out-gene + 连续状态条件化，无 time leakage）
    task = task_builder.build_heldout_task(
        adata, state_keys=["dam_score", "infl_score"], holdout_frac=0.2
    )
    X, ctx, y, has_target = task["X"], task["ctx"], task["y"], task["has_target"]
    splits = task["splits"]
    if not splits:
        print("ERROR: 真实数据无可用时序切分。")
        return
    import numpy as np
    pred_lin = np.full_like(y, np.nan)
    pred_cond = np.full_like(y, np.nan)
    predicted = np.zeros(y.shape[0], dtype=bool)
    for tr, te in splits:
        pl = baselines.fit_predict_linear(X[tr], y[tr], X[te])
        feat_tr = np.hstack([X[tr], ctx[tr]])
        feat_te = np.hstack([X[te], ctx[te]])
        pc = baselines.fit_predict_conditional(feat_tr, y[tr], feat_te)
        pred_lin[te] = pl
        pred_cond[te] = pc
        predicted[te] = True
    res = evaluate.evaluate_gate(pred_lin, pred_cond, y, predicted, n_boot=n_boot)
    print(evaluate.format_report(res, f"Gate 1 [REAL {dataset}]"))
    return res


def main():
    ap = argparse.ArgumentParser(description="Gate 1 feasibility probe")
    ap.add_argument("--mode", choices=["synthetic", "real"], default="synthetic")
    ap.add_argument("--signal", action="store_true", help="synthetic with signal")
    ap.add_argument("--null", action="store_true", help="synthetic null (no signal)")
    ap.add_argument("--dataset", default="gse174574")
    ap.add_argument("--nboot", type=int, default=1000)
    args = ap.parse_args()

    if args.mode == "synthetic":
        mode = "signal" if args.signal else ("null" if args.null else "signal")
        run_synthetic(mode, n_boot=args.nboot)
    else:
        run_real(args.dataset, n_boot=args.nboot)


if __name__ == "__main__":
    main()
