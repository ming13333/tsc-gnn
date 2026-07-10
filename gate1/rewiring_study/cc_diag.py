"""
Step C(1) 鼠内跨队列 rewiring 一致性 — 连续信号诊断（v2，修复 merge 后缀 bug + 去循环比较）
------------------------------------------------------------------------
关键设计修正：
1. 整合集(rewiring_full)的 sham 基线 = GSE174574 sham + GSE225948 sham 混合；
   单队列 c1 的 sham 基线只有 GSE174574 sham。故「整合 vs 单队列 同 transition」不干净
   （基线不同可直接翻转 dW 方向），只能作为估计量稳健性的弱附注，不能当独立重现证据。
2. 真正干净的跨队列重现 = 比较 c1 与 c2 各自「最大|ΔW| TF 排序」的重叠 / Spearman rank。
   因为 c1(sham→24h) 与 c2(sham→2d,2d→14d) 无共享 transition，边级比较无意义，
   只能在 TF 主调控因子层级跨 transition 比较。这是生物学上最有意义的重现性检验。
3. 边级 dW 相关（整合 vs 单队列 同 transition）保留为附注，但显式标注「基线混淆」限制。
"""
import pandas as pd, numpy as np
from scipy.stats import spearmanr
from collections import OrderedDict

full = pd.read_csv("rewiring_full.csv")
c1 = pd.read_csv("cc_gse174_rewiring_full.csv")   # GSE174574: sham->24h
c2 = pd.read_csv("cc_gse225_rewiring_full.csv")   # GSE225948: sham->2d, 2d->14d

# 各数据集覆盖的 transition 列
FULL_TRANS  = [("sham", "24h"), ("24h", "2d"), ("2d", "14d"), ("sham", "14d")]
C1_TRANS    = [("sham", "24h")]
C2_TRANS    = [("sham", "2d"), ("2d", "14d")]

def tf_maxdw(df, trans):
    """每个 TF 在其覆盖 transition 上的最大 |dW|（PC 校正列）。"""
    best = {}
    for (a, b) in trans:
        col = f"dW_{a}_{b}"
        if col not in df.columns:
            continue
        g = df.groupby("tf")[col].apply(lambda s: s.abs().max())
        for tf, v in g.items():
            best[tf] = max(best.get(tf, 0.0), v)
    return pd.Series(best)

def top_tfs(series, k=20):
    return set(series.sort_values(ascending=False).head(k).index)

def jac(a, b):
    if not a and not b:
        return float("nan")
    return len(a & b) / len(a | b)

s_full = tf_maxdw(full, FULL_TRANS)
s_c1   = tf_maxdw(c1, C1_TRANS)
s_c2   = tf_maxdw(c2, C2_TRANS)

top_full, top_c1, top_c2 = top_tfs(s_full), top_tfs(s_c1), top_tfs(s_c2)

print("=== Part 1. TF 主调控因子重现性（跨队列，无循环）===")
print(f"  各数据集 TF 数: full={len(s_full)} c1={len(s_c1)} c2={len(s_c2)}")
print(f"  Top20 重叠 Jaccard:")
print(f"    full vs c1 : {jac(top_full, top_c1):.3f}  (交集 {len(top_full&top_c1)})")
print(f"    full vs c2 : {jac(top_full, top_c2):.3f}  (交集 {len(top_full&top_c2)})")
print(f"    c1   vs c2 : {jac(top_c1, top_c2):.3f}  (交集 {len(top_c1&top_c2)})")
print(f"  full ∩ c1 ∩ c2 (三队列共通 Top20 TF): {sorted(top_full & top_c1 & top_c2)}")
print(f"  full ∩ c1 (共通): {sorted(top_full & top_c1)}")
print(f"  full ∩ c2 (共通): {sorted(top_full & top_c2)}")
print(f"  c1 ∩ c2   (共通): {sorted(top_c1 & top_c2)}")

print("\n=== Part 2. TF rank 相关性（Spearman，跨队列/整合）===")
def spearman_pair(sa, sb, name):
    common = sa.index.intersection(sb.index)
    if len(common) < 3:
        print(f"    {name}: 共同 TF 不足 ({len(common)})"); return
    rho, p = spearmanr(sa[common], sb[common])
    print(f"    {name}: nTF={len(common)} rho={rho:+.3f} p={p:.2e}")
spearman_pair(s_full, s_c1, "full vs c1")
spearman_pair(s_full, s_c2, "full vs c2")
spearman_pair(s_c1, s_c2, "c1   vs c2")

print("\n=== Part 3. 边级 ΔW 相关（整合 vs 单队列 同 transition，附注·基线混淆）===")
def edge_dW_corr(a, b, col):
    """同 transition 边级 dW(PC) 相关；a/b 为两个数据集；col 为该 transition 的 dW 列。"""
    a2 = a[["tf", "target", col]].copy(); a2["k"] = a2["tf"] + "->" + a2["target"]
    b2 = b[["tf", "target", col]].copy(); b2["k"] = b2["tf"] + "->" + b2["target"]
    m = a2.merge(b2, on="k", suffixes=("_a", "_b"))
    if len(m) < 3:
        print(f"    {col}: 共同边不足 ({len(m)})"); return
    xa = m[col + "_a"].values.astype(float)
    xb = m[col + "_b"].values.astype(float)
    ok = ~(np.isnan(xa) | np.isnan(xb))
    if ok.sum() < 3:
        print(f"    {col}: 有效共同边不足 ({ok.sum()})"); return
    r = np.corrcoef(xa[ok], xb[ok])[0, 1]
    agr = np.mean(np.sign(xa[ok]) == np.sign(xb[ok]))
    print(f"    {col}: n_common={int(ok.sum())} r(dW)={r:+.3f} 方向一致={agr:.2f}")
edge_dW_corr(full, c1, "dW_sham_24h")
edge_dW_corr(full, c2, "dW_2d_14d")

print("\n=== Part 4. 可解释性主轴小结（仅信号存在性，非预测精度）===")
print(f"  全整合集 Top15 TF (|ΔW|max): {sorted(s_full.sort_values(ascending=False).head(15).index)}")
print(f"  GSE174574 Top15 TF         : {sorted(s_c1.sort_values(ascending=False).head(15).index)}")
print(f"  GSE225948 Top15 TF          : {sorted(s_c2.sort_values(ascending=False).head(15).index)}")
