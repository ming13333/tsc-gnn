"""测试 24h->2d->14d 整合 + Gate 1 时间条件化。"""
import numpy as np
import time
from gate1 import preprocessing, task_builder, baselines, evaluate

t0 = time.time()
cohorts = {"sex": "male", "age": "W8", "condition": "MCAO"}
print("=== 加载整合时间序列 (cohorts=%s) ===" % cohorts)
adata = preprocessing.load_integrated_timeseries("data", cohorts=cohorts)
print("整合后细胞数: %d, 基因数: %d (%.1fs)" % (adata.n_obs, adata.n_vars, time.time()-t0))
print("time_label 分布:", adata.obs["time_label"].value_counts().to_dict())
print("study 分布:", adata.obs["study"].value_counts().to_dict())

# 只保留时间轴（stroke 24h/2d/14d），排除 sham 基线
adata = adata[adata.obs["time_label"].isin(["24h", "2d", "14d"])].copy()
print("\n时间轴细胞数: %d" % adata.n_obs)
print("时间轴分布:", adata.obs["time_label"].value_counts().to_dict())

print("\n=== prep_for_gate (z-score 批次校正 + HVG) ===")
adata = preprocessing.prep_for_gate(adata, batch_key="study", n_top_genes=2000, method="zscore")
print("prep 后: cells=%d genes=%d (%.1fs)" % (adata.n_obs, adata.n_vars, time.time()-t0))

print("\n=== 构建 held-out-gene 任务 (方案 A: 状态条件化, 无 time leakage) ===")
task = task_builder.build_heldout_task(
    adata, state_keys=["dam_score", "infl_score"], holdout_frac=0.2
)
X, ctx, y, has_target = task["X"], task["ctx"], task["y"], task["has_target"]
print("task: X(input)=%s ctx(state)=%s y(heldout)=%s" % (X.shape, ctx.shape, y.shape))
print("input_genes=%d holdout_genes=%d" % (task["n_input_genes"], task["n_holdout_genes"]))
print("has_target: %d" % has_target.sum())
print("y stats: std=%.3f mean=%.3f" % (np.std(y), np.mean(y)))

splits = task["splits"]
if not splits:
    print("ERROR: 无可用时序切分")
else:
    print("\n=== 运行 Gate 1 (真实时间轴) ===")
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
    res = evaluate.evaluate_gate(pred_lin, pred_cond, y, predicted, n_boot=1000)
    print(evaluate.format_report(res, "Gate 1 [REAL 24h->2d->14d timeseries]"))
