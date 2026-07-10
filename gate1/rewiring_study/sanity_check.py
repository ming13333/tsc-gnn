"""sanity_check.py v2 — Step B: 已知卒中基因 sanity（方法可信度生物学校验）。

v2 修复：
- 用 pooled q<0.1 做显著性筛选（不依赖 p 分辨率）
- 方向标签修正为"关联增强/关联减弱"
- 加 PC 校正 vs 未校正对比
- 富集标注"描述性"

读 rewiring_full.csv（Step A v2 产出），检验：
1. 重连 target 对卒中基因集的超几何富集
2. 已知神经/炎症 TF 是否作为重连 TF 出现
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

CSV = "C:/D 盘/科研/虚拟敲除/gate1/rewiring_study/rewiring_full.csv"
OUT = "C:/D 盘/科研/虚拟敲除/gate1/rewiring_study/SANITY_报告_v2_2026-07-09.md"
TRANSITIONS = [("sham", "24h"), ("24h", "2d"), ("2d", "14d"), ("sham", "14d")]

# curated 卒中相关基因集（mouse symbol，大写）
STROKE_GENES = {
    "小胶质/巨噬细胞": ["CSF1R", "CD68", "ITGAM", "AIF1", "TNF", "IL1B", "CX3CR1",
                      "TREM2", "CTSS", "FCGR3", "MSR1", "SIGLECH", "CD53", "CD86",
                      "ENTPD1", "HK2", "BCL2A1B"],
    "少突/髓鞘": ["MAG", "PLP1", "MBP", "GJB1", "ANK3", "MOG", "CNP", "OPALIN"],
    "星形/炎症应答": ["GFAP", "VIM", "C3", "SERPINA3N", "LCN2"],
    "神经元/轴突": ["MAP1B", "CTNND2", "NREP", "SNAP25", "SYNGR1", "NEFM", "NEFH"],
    "内皮/血管": ["CLDN5", "PECAM1", "KDR", "NOS3", "VEGFA"],
}
# 已知神经/炎症/免疫 TF
KNOWN_TFS = ["SPI1", "SOX10", "SOX2", "CEBPE", "STAT4", "IRF8", "RELA", "NFKB1",
             "CEBPB", "EGR1", "RUNX3", "GATA2", "ERG", "MAF", "NFE2"]

TF_NOTE = {
    "SPI1": "PU.1 髓系/小胶质", "SOX10": "少突胶质/髓鞘", "CEBPE": "粒细胞/炎症",
    "STAT4": "IL-12/IFNγ", "SOX2": "神经干/重编程", "RELA": "NF-kB",
    "NFKB1": "NF-kB", "IRF8": "小胶质", "CEBPB": "炎症/急性相", "EGR1": "即刻早期",
    "RUNX3": "CD8/T", "GATA2": "造血", "ERG": "内皮", "MAF": "免疫", "NFE2": "红系",
}


def _enrich(targets_sig, all_targets, gene_set):
    """超几何富集（Haldane-Anscombe 校正）。"""
    gs = set(g.upper() for g in gene_set)
    pop = set(g.upper() for g in all_targets)
    gs_in_pop = gs & pop
    if len(gs_in_pop) == 0:
        return None
    M = len(pop)
    n = len(gs_in_pop)
    N = len(targets_sig)
    k = len(set(g.upper() for g in targets_sig) & gs_in_pop)
    if N == 0 or k == 0:
        return dict(hit=k, N=N, M=M, n=n, OR=0.0, p=1.0)
    p = hypergeom.sf(k - 1, M, n, N)
    k1, nk1 = k + 0.5, (N - k) + 0.5
    nk2, mnnk2 = (n - k) + 0.5, (M - N - (n - k)) + 0.5
    OR = (k1 / nk1) / (nk2 / mnnk2)
    return dict(hit=k, N=N, M=M, n=n, OR=OR, p=p)


def main():
    df = pd.read_csv(CSV)
    all_targets = df["target"].unique().tolist()
    L = ["# Step B v2 — 已知卒中基因 Sanity Check", ""]
    L.append(f"读入 Step A v2 的 rewiring_full.csv（{df.shape[0]} 条边）")
    L.append("")
    L.append("## 显著性筛选标准")
    L.append("")
    L.append("- **Pooled q < 0.1**：置换 pooled FDR（不依赖 p 分辨率），"
             "用于正式富集分析")
    L.append("- **候选边**（p<0.05 且 |ΔW|>0.15）：用于描述性补充，"
             "标注为探索性")
    L.append("")

    # ── 1. 富集分析 ──
    L.append("## 1. 重连 target 对卒中基因集的富集（pooled q<0.1）")
    L.append("")
    L.append("| 转移 | 基因集 | 命中/显著 | 背景集大小 | OR | p |")
    L.append("|---|---|---|---|---|---|")
    for (t1, t2) in TRANSITIONS:
        sig = df[df[f"q_pooled_{t1}_{t2}"] < 0.1]
        targets_sig = sig["target"].tolist()
        for gset_name, gs in STROKE_GENES.items():
            r = _enrich(targets_sig, all_targets, gs)
            if r is None:
                continue
            mark = " **" if (r["p"] < 0.05 and r["OR"] > 1) else ""
            L.append(
                f"| {t1}→{t2} | {gset_name} | {r['hit']}/{r['N']} | "
                f"{r['n']} | {r['OR']:.2f} | {r['p']:.1e} |{mark}")
    L.append("")
    L.append("_注：超几何 p 为描述性（候选 target 经 |ΔW|&p 选择后非随机抽样），"
             "须结合阳性/阴性对照解读。_")
    L.append("")

    # ── 1b. 候选边富集（描述性补充）──
    L.append("## 1b. 候选边 target 富集（p<0.05 & |ΔW|>0.15，描述性）")
    L.append("")
    L.append("| 转移 | 基因集 | 命中/候选 | OR | p |")
    L.append("|---|---|---|---|---|")
    for (t1, t2) in TRANSITIONS:
        cand = df[(df[f"p_{t1}_{t2}"] < 0.05) &
                  (df[f"dW_{t1}_{t2}"].abs() > 0.15)]
        targets_cand = cand["target"].tolist()
        for gset_name, gs in STROKE_GENES.items():
            r = _enrich(targets_cand, all_targets, gs)
            if r is None:
                continue
            mark = " **" if (r["p"] < 0.05 and r["OR"] > 1) else ""
            L.append(
                f"| {t1}→{t2} | {gset_name} | {r['hit']}/{r['N']} | "
                f"{r['OR']:.2f} | {r['p']:.1e} |{mark}")
    L.append("")

    # ── 2. 已知 TF 出现检查 ──
    L.append("## 2. 已知神经/炎症 TF 是否作为重连 TF 出现")
    L.append("")
    L.append("| TF | TF功能 | 在 TF 列 | 重连边数 | "
             "最显著转移(min pooled q) | |ΔW| max |")
    L.append("|---|---|---|---|---|---|")
    tfs_upper = set(df["tf"].str.upper())
    for g in KNOWN_TFS:
        in_tf = g in tfs_upper
        sub = df[df["tf"].str.upper() == g]
        n_edges = len(sub)
        if n_edges > 0:
            min_q = min(sub[f"q_pooled_{t1}_{t2}"].min()
                        for (t1, t2) in TRANSITIONS)
            best_tr = None
            for (t1, t2) in TRANSITIONS:
                if sub[f"q_pooled_{t1}_{t2}"].min() == min_q:
                    best_tr = f"{t1}→{t2}"
            max_dw = max(sub[f"dW_{t1}_{t2}"].abs().max()
                         for (t1, t2) in TRANSITIONS)
            note = TF_NOTE.get(g, "")
            L.append(f"| {g} | {note} | 是 | {n_edges} | "
                     f"{best_tr} (q={min_q:.1e}) | {max_dw:.2f} |")
        else:
            note = TF_NOTE.get(g, "")
            L.append(f"| {g} | {note} | 否 | 0 | - | - |")
    L.append("")

    # ── 3. PC 校正前后 Top 边方向一致性 ──
    L.append("## 3. PC 校正前后 ΔW 方向一致性")
    L.append("")
    L.append("| 转移 | 同向比例 | Pearson(raw, PC) |")
    L.append("|---|---|---|")
    for (t1, t2) in TRANSITIONS:
        raw = df[f"dW_raw_{t1}_{t2}"].values
        pc = df[f"dW_{t1}_{t2}"].values
        same_sign = np.mean(np.sign(raw) == np.sign(pc))
        corr = np.corrcoef(raw, pc)[0, 1]
        L.append(f"| {t1}→{t2} | {same_sign:.2%} | {corr:.2f} |")
    L.append("")
    L.append("- 高同向比例（>80%）说明 PC 校正未翻转主要信号方向，"
             "结果稳健。")
    L.append("")

    # ── 4. 解读 ──
    L.append("## 4. 解读")
    L.append("")
    L.append("- 若重连 target 富集于**小胶质/少突/神经元**等卒中核心细胞类型，"
             "说明 rewiring 方法能 recover 已知生物学（可信度证据）。")
    L.append("- 已知神经/炎症 TF（Spi1/Sox10/Sox2/RELA…）若作为重连 TF 出现，"
             "将方法与卒中神经生物学直接挂钩。")
    L.append("- PC 校正后方向一致性高 → 信号非纯组成驱动；"
             "若低 → 须降级为'关联重连'并列 limitation。")
    L.append("- 富集 p 标注为描述性，须补阳性/阴性对照（下一步）。")
    L.append("")
    with open(OUT, "w") as f:
        f.write("\n".join(L))
    print(f"[done] {OUT}")


if __name__ == "__main__":
    main()
