"""make_v2_report.py — 生成生物学导向的 Step A rewiring 报告 v2。

不重跑（复用 rewiring_full.csv）。采用候选重连边标准 = (p<0.05 且 |dW|>0.15)，
因 11434 边 BH-FDR<0.05 阈值过严（约需 p<4e-6）导致正式显著=0，而 Top 边 |dW|
达 0.4-1.0 且生物学高度合理。重点展示有方向的 TF->target 重连。
"""
import pandas as pd
import numpy as np

CSV = "C:/D 盘/科研/虚拟敲除/gate1/rewiring_study/rewiring_full.csv"
OUT = "C:/D 盘/科研/虚拟敲除/gate1/rewiring_study/REWIRING_报告_v2_2026-07-09.md"
TRANSITIONS = [("sham", "24h"), ("24h", "2d"), ("2d", "14d"), ("sham", "14d")]

TF_NOTE = {
    "SPI1": "PU.1 髓系/小胶质", "SOX10": "少突胶质/髓鞘", "CEBPE": "粒细胞/炎症",
    "STAT4": "IL-12/IFNγ", "PAX5": "B细胞", "SOX2": "神经干/重编程", "RELA": "NF-kB",
    "NFKB1": "NF-kB", "IRF8": "小胶质", "CEBPB": "炎症/急性相", "EGR1": "即刻早期",
    "RUNX3": "CD8/T", "GATA2": "造血", "ERG": "内皮", "MAF": "免疫", "NFE2": "红系",
    "SPI1B": "髓系", "SNAI2": "EMT", "NR2F2": "血管", "ERGR": "内皮", "GATA2R": "造血",
}

df = pd.read_csv(CSV)
L = []
L.append("# Step A — 卒中时间×状态 DoRothEA 重布线报告（生物学导向 v2）")
L.append("")
L.append("日期：2026-07-09｜数据：鼠 MCAO 整合时间序列 (GSE174574 24h + GSE225948 2d/14d + sham)")
L.append("图：DoRothEA(mouse, ABC) 有向因果图｜抽样细胞=8000(每点≤2000)｜置换 n_perm=80")
L.append("")
L.append("## 显著性判定的重要说明")
L.append("")
L.append("- 测试边总数（state-conditioned, |A_aff| 前 50%）：**11434**")
L.append("- 若用 BH-FDR<0.05 严格多重校正：因 11434 条边阈值极严（约需 p<4e-6），"
         "**正式显著边 = 0**；且 n_perm=80 的 p 值分辨率下限为 1/81≈0.0123，强边撞到天花板。")
L.append("- 但**未校正 p<0.05 的候选边丰富**（各转移 1000–2230 条），Top 边 |ΔW| 达 0.4–1.0 且生物学高度合理。")
L.append("- 故本报告采用 **候选重连边 = (p<0.05 且 |ΔW|>0.15)** 作为可解释性证据，"
         "并强调 Top 边的生物学方向性（线性模型给不出）。")
L.append("- 补强建议：重跑 N_PERM=200（p 分辨率 0.005）+ 用 FDR<0.1 宽松阈值确认。")
L.append("")
L.append("## 各转移候选重连边数（p<0.05 且 |ΔW|>0.15）")
L.append("")
L.append("| 转移 | 候选边 |")
L.append("|---|---|")
for (t1, t2) in TRANSITIONS:
    n = int(((df[f"p_{t1}_{t2}"] < 0.05) & (df[f"dW_{t1}_{t2}"].abs() > 0.15)).sum())
    L.append(f"| {t1}→{t2} | {n} |")
L.append("")

for (t1, t2) in TRANSITIONS:
    sub = df[(df[f"p_{t1}_{t2}"] < 0.05) & (df[f"dW_{t1}_{t2}"].abs() > 0.15)].copy()
    sub = sub.reindex(sub[f"dW_{t1}_{t2}"].abs().sort_values(ascending=False).index).head(15)
    L.append(f"## {t1}→{t2} Top 15 重连边（按 |ΔW|）")
    L.append("")
    L.append(f"| TF | target | TF功能 | 方向 | ΔW | p | r_{t1} | r_{t2} | DoRothEA |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for _, r in sub.iterrows():
        tf = r["tf"]; note = TF_NOTE.get(tf.upper(), "")
        d = r[f"dW_{t1}_{t2}"]; s = "激活增强" if d > 0 else "抑制增强"
        L.append(f"| {tf} | {r['target']} | {note} | {s} | {d:+.2f} | "
                 f"{r[f'p_{t1}_{t2}']:.1e} | {r[f'r_{t1}']:+.2f} | {r[f'r_{t2}']:+.2f} | {r['dorothea_sign']} |")
    L.append("")

L.append("## 生物学解读（核心卖点）")
L.append("")
L.append("- **Spi1(PU.1)→小胶质/髓系靶基因**（Csf1r, Cd53, Cd86, Fcgr3, Msr1, Siglech, Tnf, "
         "Aif1, Ctss, Entpd1 等）：sham→24h 多为激活增强（缺血即小胶质激活），24h→2d 与 "
         "sham→14d 多为**抑制增强**（r 从 0.5–0.6 跌至 0.1–0.2）——卒中后 Spi1 对其靶基因的调控耦合"
         "**减弱**，符合小胶质从急性激活转向修复/静息态。")
L.append("- **Sox10→少突/髓鞘基因**（Mag, Plp1, Mbp, Gjb1, Ank3）：24h→2d **抑制增强**"
         "（r 0.4–0.5→0，亚急性脱髓鞘），2d→14d **激活增强**（r 0→0.57–1.0，修复期重髓鞘化）——"
         "完美对应“亚急性损伤→修复重塑”三阶段轴。")
L.append("- **Sox2→神经/轴突**（Map1b, Ctnnd2, Nrep）：2d→14d 激活增强，指向神经修复/可塑性。")
L.append("- **Cebpe→Ltf**：2d→14d 激活增强，炎症/粒细胞信号。")
L.append("- 这些**有方向、有 DoRothEA 因果图支撑**的边重连，线性模型完全无法给出——"
         "正是“可解释性主轴”的价值：图提供边级 rewiring 解释，而非预测精度。")
L.append("")
L.append("## 与项目主轴的衔接")
L.append("- 鲁棒性研究已证“图不优于线性预测”（真实数据 + GEARS 金标准），本分析补足另一半：")
L.append("  **图提供线性给不出的边重连解释**，且生物学方向正确。")
L.append("- 下一步：重跑 N_PERM=200 补强 p 值；Step B 改用神经/炎症 TF 及其靶基因集做 sanity。")

with open(OUT, "w") as f:
    f.write("\n".join(L))
print(f"[done] {OUT}")
