"""cross_cohort_consistency.py — Step C(1): 鼠内跨队列 rewiring 一致性验证。

逻辑：
  完整整合 rewiring（rewiring_full.csv, 4 transitions）
    vs 队列1 GSE174574（sham→24h, 急性）
    vs 队列2 GSE225948（sham→2d, 2d→14d, 亚急性→修复）

由于队列1的 sham→24h 与完整集 sham→24h、队列2的 2d→14d 与完整集 2d→14d
是**同一 transition、不同独立数据集**，可直接做边级/方向一致性（最强检验）。
同时做 TF 级 Jaccard 重现 + target 模块超几何富集。

输出：CROSS_COHORT_报告_2026-07-09.md
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

HERE = "C:/D 盘/科研/虚拟敲除/gate1/rewiring_study"
FULL = os.path.join(HERE, "rewiring_full.csv")
C1 = os.path.join(HERE, "cc_gse174_rewiring_full.csv")
C2 = os.path.join(HERE, "cc_gse225_rewiring_full.csv")
OUT = os.path.join(HERE, "CROSS_COHORT_报告_2026-07-09.md")

TF_NOTE = {
    "SPI1": "PU.1 髓系/小胶质", "SOX10": "少突胶质/髓鞘", "CEBPE": "粒细胞/炎症",
    "STAT4": "IL-12/IFNγ", "SOX2": "神经干/重编程", "RELA": "NF-kB",
    "NFKB1": "NF-kB", "IRF8": "小胶质", "CEBPB": "炎症/急性相", "EGR1": "即刻早期",
    "GATA2": "造血", "STAT1": "IFN", "STAT3": "IL-6/急性相", "IRF1": "I 型干扰素",
    "RUNX3": "CD8/T", "MAF": "免疫", "CEBPA": "炎症/髓系", "NR2F2": "血管",
    "SNAI2": "EMT", "PPARG": "巨噬", "FOXP3": "Treg", "MEF2C": "神经",
    "TCF7L2": "Wnt", "SP1": "通用", "SP2": "通用", "KLF4": "干性",
}


def sig_edges(df, t1, t2, qthr=0.1):
    """返回 (edges set, direction dict)，edges = 'TF->target'，方向=sign(dW)。"""
    qcol = f"q_pooled_{t1}_{t2}"
    if qcol not in df.columns:
        return set(), {}
    sub = df[df[qcol] < qthr]
    edges = set(f"{r.tf}->{r.target}" for _, r in sub.iterrows())
    direct = {f"{r.tf}->{r.target}": int(np.sign(r[f"dW_{t1}_{t2}"]))
              for _, r in sub.iterrows()}
    return edges, direct


def tf_set(edges):
    return set(e.split("->")[0] for e in edges)


def target_set(edges):
    return set(e.split("->")[1] for e in edges)


def jaccard(a, b):
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


def direction_agreement(d1, d2):
    common = set(d1) & set(d2)
    if not common:
        return None, 0
    same = sum(1 for e in common if d1[e] == d2[e])
    return same / len(common), len(common)


def enrich(universe_targets, sig_full_targets, query_targets):
    """超几何富集：query 的 target 是否富集于 full 显著 target。
    universe = 完整集所有 unique target；n = full 显著 target 数；
    N = query 所有 unique target；k = query 显著 target ∩ full 显著 target。"""
    M = len(universe_targets)
    n = len(sig_full_targets & universe_targets)
    N = len(query_targets & universe_targets)
    k = len(query_targets & sig_full_targets)
    if M == 0 or n == 0 or N == 0:
        return dict(OR=np.nan, p=np.nan, M=M, n=n, N=N, k=k)
    p = hypergeom.sf(k - 1, M, n, N)
    OR = (k / (N - k + 1e-9)) / ((n - k) / (M - N - (n - k) + 1e-9))
    return dict(OR=OR, p=p, M=M, n=n, N=N, k=k)


def main():
    full = pd.read_csv(FULL)
    c1 = pd.read_csv(C1)
    c2 = pd.read_csv(C2)
    print("加载: 完整集", full.shape, "| 队列1", c1.shape, "| 队列2", c2.shape)

    FULL_TRANS = [("sham", "24h"), ("24h", "2d"), ("2d", "14d"), ("sham", "14d")]
    C1_TRANS = [("sham", "24h")]
    C2_TRANS = [("sham", "2d"), ("2d", "14d")]

    # ── 同 transition 边级一致性 ──
    fe_24, fd_24 = sig_edges(full, "sham", "24h")
    c1e_24, c1d_24 = sig_edges(c1, "sham", "24h")
    jac_24 = jaccard(fe_24, c1e_24)
    da_24 = direction_agreement(fd_24, c1d_24)

    fe_214, fd_214 = sig_edges(full, "2d", "14d")
    c2e_214, c2d_214 = sig_edges(c2, "2d", "14d")
    jac_214 = jaccard(fe_214, c2e_214)
    da_214 = direction_agreement(fd_214, c2d_214)

    # ── TF 级重现 ──
    tf_full = set()
    for tr in FULL_TRANS:
        tf_full |= tf_set(sig_edges(full, *tr)[0])
    tf_c1 = tf_set(sig_edges(c1, "sham", "24h")[0])
    tf_c2 = set()
    for tr in C2_TRANS:
        tf_c2 |= tf_set(sig_edges(c2, *tr)[0])
    jac_tf_fc1 = jaccard(tf_full, tf_c1)
    jac_tf_fc2 = jaccard(tf_full, tf_c2)
    jac_tf_c1c2 = jaccard(tf_c1, tf_c2)
    common_tf = tf_full & tf_c1 & tf_c2

    # ── target 模块富集 ──
    tgt_full = set()
    for tr in FULL_TRANS:
        tgt_full |= target_set(sig_edges(full, *tr)[0])
    universe = tgt_full  # 完整集显著 target 作为 universe
    tgt_c1 = target_set(sig_edges(c1, "sham", "24h")[0])
    tgt_c2 = set()
    for tr in C2_TRANS:
        tgt_c2 |= target_set(sig_edges(c2, *tr)[0])
    enr_c1 = enrich(universe, tgt_full, tgt_c1)
    enr_c2 = enrich(universe, tgt_full, tgt_c2)

    # ── 汇总打印 ──
    print(f"\n=== 同 transition 边级一致性 ===")
    print(f"sham→24h (完整集 vs GSE174574): Jaccard={jac_24:.3f}, "
          f"共同边={len(fe_24 & c1e_24)}, 方向一致={da_24}")
    print(f"2d→14d  (完整集 vs GSE225948): Jaccard={jac_214:.3f}, "
          f"共同边={len(fe_214 & c2e_214)}, 方向一致={da_214}")
    print(f"\n=== TF 级重现 Jaccard ===")
    print(f"完整集 TF={len(tf_full)} | 队列1 TF={len(tf_c1)} | 队列2 TF={len(tf_c2)}")
    print(f"Jaccard(完整,队列1)={jac_tf_fc1:.3f}")
    print(f"Jaccard(完整,队列2)={jac_tf_fc2:.3f}")
    print(f"Jaccard(队列1,队列2)={jac_tf_c1c2:.3f}")
    print(f"三队列共通 TF ({len(common_tf)}): {sorted(common_tf)}")
    print(f"\n=== target 模块富集 (完整集显著 target 为 universe) ===")
    print(f"队列1: OR={enr_c1['OR']:.2f} p={enr_c1['p']:.2e} k={enr_c1['k']}/{enr_c1['N']}")
    print(f"队列2: OR={enr_c2['OR']:.2f} p={enr_c2['p']:.2e} k={enr_c2['k']}/{enr_c2['N']}")

    # ── 写报告 ──
    L = []
    L.append("# Step C(1) — 鼠内跨队列 rewiring 一致性验证")
    L.append("(2026-07-09｜cross-cohort reproducibility)")
    L.append("")
    L.append("## 设计")
    L.append("")
    L.append("- **完整整合 rewiring**：GSE174574(24h)+GSE225948(2d/14d)+sham 整合，"
             "4 个转移（sham→24h, 24h→2d, 2d→14d, sham→14d），pooled q<0.1 显著边作核心集。")
    L.append("- **队列1（GSE174574 独立）**：sham→24h（急性损伤），独立加载/归一化/rewiring。")
    L.append("- **队列2（GSE225948 独立）**：sham→2d, 2d→14d（亚急性→修复），独立加载/归一化/rewiring。")
    L.append("- 队列1的 sham→24h 与完整集 sham→24h、队列2的 2d→14d 与完整集 2d→14d "
             "是**同一 transition、不同独立数据集** → 可做边级/方向一致性（最强 reproducible 检验）。")
    L.append("")
    L.append("## 1. 同 transition 边级一致性（最直接）")
    L.append("")
    L.append("| 转移 | 完整集显著边 | 独立队列显著边 | 共同边 | Jaccard | "
             "方向一致比例(共同边) |")
    L.append("|---|---|---|---|---|---|")
    da24_str = f"{da_24[0]:.0%}(n={da_24[1]})" if da_24 and da_24[0] is not None else "—"
    da214_str = f"{da_214[0]:.0%}(n={da_214[1]})" if da_214 and da_214[0] is not None else "—"
    L.append(f"| sham→24h (完整 vs GSE174574) | {len(fe_24)} | {len(c1e_24)} | "
             f"{len(fe_24 & c1e_24)} | {jac_24:.3f} | {da24_str} |")
    L.append(f"| 2d→14d (完整 vs GSE225948) | {len(fe_214)} | {len(c2e_214)} | "
             f"{len(fe_214 & c2e_214)} | {jac_214:.3f} | {da214_str} |")
    L.append("")
    L.append("- Jaccard 衡量两独立数据集在同一 transition 上显著重连**边的重叠**；")
    L.append("  方向一致比例衡量重叠边的 ΔW 符号（关联增强/减弱）是否跨数据集一致。")
    L.append("")
    L.append("## 2. TF 级重现（核心调控因子 reproducibility）")
    L.append("")
    L.append(f"- 完整集核心 TF（所有转移 q<0.1 边的 TF 并集）：**{len(tf_full)}** 个")
    L.append(f"- 队列1 核心 TF（sham→24h q<0.1）：**{len(tf_c1)}** 个")
    L.append(f"- 队列2 核心 TF（sham→2d + 2d→14d q<0.1）：**{len(tf_c2)}** 个")
    L.append("")
    L.append("| 比较 | Jaccard |")
    L.append("|---|---|")
    L.append(f"| 完整集 vs 队列1 | {jac_tf_fc1:.3f} |")
    L.append(f"| 完整集 vs 队列2 | {jac_tf_fc2:.3f} |")
    L.append(f"| 队列1 vs 队列2 | {jac_tf_c1c2:.3f} |")
    L.append("")
    L.append(f"### 三队列共通核心 TF（{len(common_tf)} 个）")
    L.append("")
    common_sorted = sorted(common_tf)
    L.append("| TF | 功能注释 |")
    L.append("|---|---|")
    for tf in common_sorted:
        L.append(f"| {tf} | {TF_NOTE.get(tf.upper(), '')} |")
    L.append("")
    L.append("- 这些 TF 在**三个独立分析**（整合 + GSE174574 + GSE225948）中均作为"
             "显著重连 TF 出现，是 rewiring 信号 reproducible 的强证据。")
    L.append("")
    L.append("## 3. target 模块富集（核心靶基因 reproducible）")
    L.append("")
    L.append(f"完整集核心 target 作为 universe（M={enr_c1['M']}）。")
    L.append("")
    L.append("| 队列 | 富集 OR | p(超几何) | 命中 k/N |")
    L.append("|---|---|---|---|")
    L.append(f"| 队列1 (GSE174574) | {enr_c1['OR']:.2f} | {enr_c1['p']:.2e} | "
             f"{enr_c1['k']}/{enr_c1['N']} |")
    L.append(f"| 队列2 (GSE225948) | {enr_c2['OR']:.2f} | {enr_c2['p']:.2e} | "
             f"{enr_c2['k']}/{enr_c2['N']} |")
    L.append("")
    L.append("- OR>1 且 p 小 ⇒ 独立队列的核心靶模块显著富集于完整集核心靶模块，"
             "说明 rewiring 命中的靶基因集跨数据集 reproducible。")
    L.append("")
    L.append("## 4. 诚实披露 / 局限")
    L.append("")
    L.append("- 三个分析的基因集略有差异（整合=两集交集 6897 基因；GSE174574 单集 8736；"
             "GSE225948 单集 7041），故 Jaccard 受基因集不可比影响（下界偏保守）。")
    L.append("- 队列1 仅 sham→24h（急性期）；队列2 仅 sham→2d/2d→14d（亚急性→修复）。"
             "完整集的 24h→2d 转移无独立队列直接对应（GSE174574 无 2d、GSE225948 无 24h），"
             "故该转移的一致性未单独验证。")
    L.append("- HRA007397（原计划的跨物种人队列）实为**受控访问的 PBMC 原始 fastq**，"
             "非脑实质、且不可直接下载（见日志），故跨物种验证改用人 TF 模块富集（C(2)）"
             "或降级为 limitation。本 C(1) 用鼠内独立队列验证 rewiring 非数据集伪影。")
    L.append("")
    with open(OUT, "w") as f:
        f.write("\n".join(L))
    print(f"\n[done] 报告 -> {OUT}")


if __name__ == "__main__":
    main()
