import pandas as pd, numpy as np
from scipy.stats import spearmanr

full = pd.read_csv("rewiring_full.csv")
c1 = pd.read_csv("cc_gse174_rewiring_full.csv")
c2 = pd.read_csv("cc_gse225_rewiring_full.csv")


def edge_dW_corr(a, b, col_a, col_b):
    a2 = a.copy(); a2["k"] = a2["tf"] + "->" + a2["target"]
    b2 = b.copy(); b2["k"] = b2["tf"] + "->" + b2["target"]
    m = a2.merge(b2, on="k", suffixes=("_a", "_b"))
    if len(m) < 3:
        return None
    x = m[col_a].values.astype(float)
    y = m[col_b].values.astype(float)
    ok = ~(np.isnan(x) | np.isnan(y))
    if ok.sum() < 3:
        return None
    r = np.corrcoef(x[ok], y[ok])[0, 1]
    # 也看 raw dW 相关
    xr = m["dW_raw_" + col_a.split("_", 1)[1]].values.astype(float) if "dW_raw_" + col_a.split("_", 1)[1] in m else np.nan
    return dict(n=len(m), ncommon=int(ok.sum()), r_pc=r)


def tf_rank_corr(a, b, col_a, col_b):
    ga = a.groupby("tf")[col_a].apply(lambda s: s.abs().max())
    gb = b.groupby("tf")[col_b].apply(lambda s: s.abs().max())
    common = ga.index.intersection(gb.index)
    if len(common) < 3:
        return None
    rho, p = spearmanr(ga[common], gb[common])
    return dict(ntf=len(common), rho=rho, p=p)


print("=== 边级 ΔW(PC-corrected) 相关性（同 transition 不同数据集）===")
r1 = edge_dW_corr(full, c1, "dW_sham_24h", "dW_sham_24h")
print(f"  sham->24h  完整集 vs GSE174574: {r1}")
r2 = edge_dW_corr(full, c2, "dW_2d_14d", "dW_2d_14d")
print(f"  2d->14d    完整集 vs GSE225948: {r2}")

print("\n=== TF 级 |ΔW|max rank 相关性（Spearman）===")
t1 = tf_rank_corr(full, c1, "dW_sham_24h", "dW_sham_24h")
print(f"  sham->24h  完整集 vs GSE174574: {t1}")
t2 = tf_rank_corr(full, c2, "dW_2d_14d", "dW_2d_14d")
print(f"  2d->14d    完整集 vs GSE225948: {t2}")

print("\n=== 完整集 Top 20 TF（按 |ΔW|max 跨4转移）是否在独立队列重现 ===")
trans_all = [("sham", "24h"), ("24h", "2d"), ("2d", "14d"), ("sham", "14d")]
full_max = {}
for (t1, t2) in trans_all:
    col = f"dW_{t1}_{t2}"
    g = full.groupby("tf")[col].apply(lambda s: s.abs().max())
    for tf, v in g.items():
        full_max[tf] = max(full_max.get(tf, 0), v)
top_full = sorted(full_max, key=lambda x: -full_max[x])[:20]
print("  完整集 Top20 TF:", top_full)
c1_max = c1.groupby("tf")["dW_sham_24h"].apply(lambda s: s.abs().max()).to_dict()
c2_max = {}
for (t1, t2) in [("sham", "2d"), ("2d", "14d")]:
    col = f"dW_{t1}_{t2}"
    g = c2.groupby("tf")[col].apply(lambda s: s.abs().max())
    for tf, v in g.items():
        c2_max[tf] = max(c2_max.get(tf, 0), v)
print("  Top20 在 GSE174574(sham24h) 出现:", sum(1 for tf in top_full if tf in c1_max))
print("  Top20 在 GSE225948(sham2d/2d14d) 出现:", sum(1 for tf in top_full if tf in c2_max))
print("  GSE174574 命中的 Top20 TF:", [tf for tf in top_full if tf in c1_max])
print("  GSE225948 命中的 Top20 TF:", [tf for tf in top_full if tf in c2_max])
