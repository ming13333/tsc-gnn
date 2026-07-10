"""
audit_positive_control.py — held-out 任务「敏感性」正对照（审计用）。

问题：真实数据上 rel_imp≈0%（FAIL）。这是「真阴性」（确实无状态条件信号），
还是「假阴性」（任务构造本身不敏感，测不出信号）？

方法：复用与 run_gate1_timeseries.py 完全相同的真实数据 prep。
然后对 held-out 目标基因**人工注入**已知的、由状态 score 驱动的扰动：
    y_inj[:, j] += signal * ctx[:, 0]      (j ∈ held-out 基因集)
重跑 线性 vs 条件化：
  · 若条件化显著恢复 rel_imp>0 → 任务敏感，真实 FAIL 是真阴性（可信）。
  · 若仍 ≈0 → 任务构造本身「瞎」（假阴性），方法无效，必须重设。

设计合法性（非 leakage）：
  - 注入信号只进入 held-out 目标基因 y，不进入输入基因 X。
  - 线性模型只看到 X（不含 ctx）；条件化模型看到 [X | ctx]。
  - 只有条件化模型能利用 ctx 还原注入项 → 公平的能力测试。
"""
import numpy as np
import time
from gate1 import preprocessing, task_builder, baselines, evaluate

t0 = time.time()
cohorts = {"sex": "male", "age": "W8", "condition": "MCAO"}
print("=== 加载整合时间序列 ===")
adata = preprocessing.load_integrated_timeseries("data", cohorts=cohorts)
adata = adata[adata.obs["time_label"].isin(["24h", "2d", "14d"])].copy()
adata = preprocessing.prep_for_gate(adata, batch_key="study", n_top_genes=2000, method="zscore")
print("prep 后: cells=%d genes=%d (%.1fs)" % (adata.n_obs, adata.n_vars, time.time() - t0))

task = task_builder.build_heldout_task(
    adata, state_keys=["dam_score", "infl_score"], holdout_frac=0.2
)
X, ctx, y, splits = task["X"], task["ctx"], task["y"], task["splits"]
print("task: X(input)=%s ctx(state)=%s y(heldout)=%s" % (X.shape, ctx.shape, y.shape))


def run_evaluate(y_target):
    pred_lin = np.full_like(y_target, np.nan)
    pred_cond = np.full_like(y_target, np.nan)
    predicted = np.zeros(y_target.shape[0], dtype=bool)
    for tr, te in splits:
        pl = baselines.fit_predict_linear(X[tr], y_target[tr], X[te])
        feat_tr = np.hstack([X[tr], ctx[tr]])
        feat_te = np.hstack([X[te], ctx[te]])
        pc = baselines.fit_predict_conditional(feat_tr, y_target[tr], feat_te)
        pred_lin[te] = pl
        pred_cond[te] = pc
        predicted[te] = True
    return evaluate.evaluate_gate(pred_lin, pred_cond, y_target, predicted, n_boot=1000)


print("\n--- 基线（真实目标，无注入）---")
r0 = run_evaluate(y)
print(evaluate.format_report(r0, "POS-CTRL baseline (real)"))

print("\n--- 注入 state 驱动信号，逐强度 ---")
for signal in [0.5, 1.0, 2.0, 4.0]:
    y_inj = y + signal * ctx[:, 0:1]   # 注入 dam_score 驱动分量到 held-out 基因
    r = run_evaluate(y_inj)
    print(f"[signal={signal}] rel_imp={r['rel_improvement']*100:.1f}%  "
          f"bootCI=[{r['boot_ci_lo']*100:.1f}%,{r['boot_ci_hi']*100:.1f}%]  "
          f"verdict={'PASS' if r['verdict_pass'] else 'FAIL'}")
    if r["verdict_pass"]:
        print(f"  >>> 任务在 signal={signal} 即敏感：真实 FAIL 为真阴性，可信。")
        break
else:
    print("  >>> 所有强度均未恢复 rel_imp>0：任务构造疑似「假阴性」，需重设！")
