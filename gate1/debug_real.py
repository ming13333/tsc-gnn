import numpy as np
from gate1 import preprocessing, task_builder, baselines

adata = preprocessing.load_processed("gse174574")
hvg = adata.var.get("highly_variable", None)
gene_mask = hvg.values if hvg is not None else None
task = task_builder.build_real_task(
    adata, time_key="time_label", state_keys=["dam_score", "infl_score"], gene_mask=gene_mask
)
X, ctx, y, has_target = task["X"], task["ctx"], task["y"], task["has_target"]
print("X shape:", X.shape, "y shape:", y.shape)
print("n_times:", task["n_times"], "time_labels:", task["time_labels"])
print("has_target sum:", has_target.sum(), "/", has_target.shape[0])
print("y stats (target cells):",
      "min=%.4f mean=%.4f max=%.4f std=%.4f" % (
          np.nanmin(y[has_target]), np.nanmean(y[has_target]),
          np.nanmax(y[has_target]), np.nanstd(y[has_target])))
print("ctx shape:", ctx.shape, "ctx std per col (first 6):",
      np.round(ctx.std(axis=0)[:6], 3))
# 快速线性拟合残差
splits = task["splits"]
tr, te = splits[0]
pl = baselines.fit_predict_linear(X[tr], y[tr], X[te])
# 测试集真实 y
yt = y[te]
print("test mse (linear, raw):", float(np.mean((pl - yt) ** 2)))
print("test pred std:", float(np.std(pl)), "true std:", float(np.std(yt)))
# 条件化
feat_tr = np.hstack([X[tr], ctx[tr]])
feat_te = np.hstack([X[te], ctx[te]])
pc = baselines.fit_predict_conditional(feat_tr, y[tr], feat_te)
print("test mse (cond, raw):", float(np.mean((pc - yt) ** 2)))
