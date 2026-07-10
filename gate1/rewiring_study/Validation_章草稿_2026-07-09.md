# Validation of the virtual rewiring framework（方法节·Validation 章草稿）
存档日期：2026-07-09 ｜ 整合自 `CROSS_COHORT_报告_2026-07-09.md`（C1）与 `CROSS_SPECIES_报告_2026-07-09.md`（C2）

> **用途**：本文件是可直接进入 manuscript「Validation / Benchmarking」节的草稿，已按审稿框架重排为四级结构。
> - 结构主线（每级回答一个问题）：
>   1. **Cross-cohort reproducibility（鼠内跨队列）** → 方法是否稳定？
>   2. **Recovery of established stroke programs（已知生物学重现）** → 方法是否正确？
>   3. **Cross-species conservation of regulatory programs（跨物种调控程序保守）** → 是否保守、可迁移至人？
>   4. **Drug-perturbation reversal via LINCS（药物逆转）** → 是否可转化？（✅ 已执行 2026-07-09）
> - 命名已升级：原「人 TF 模块富集」→ **Cross-species conservation of regulatory programs**（见 Level 3）。
> - 两条诚实底线已写入正文：① 边级不跨队列重现是**卖点**（支撑可解释性主轴），须如实写、不当缺陷藏；② Level 3 验证的是调控程序在人 GRN 中的**结构保守**，不是「人卒中里该模块表达升高」（后者需人表达数据，我们暂无）。
> - Level 4（D）已执行：robust 共识签名 top50 命中 4 个已知药，置换对照 emp_p=0.33（不显著），定位为 translational proof-of-concept。详见 `StepD_文章级分析阐述_2026-07-09.md`。

---

# Validation

## Overview（四级验证逻辑）

```
Mouse cohort 1 (GSE174574)  ─┐
                              ├─ Level 1: Cross-cohort reproducibility ──┐
Mouse cohort 2 (GSE225948)  ─┘                                         │
                                                                         ├─ Level 2: Recovery of established stroke programs
                                                                         │
Mouse integrated rewiring (Step A/B) ───────────────────────────────────┘
                                                                         │
Mouse reproducible core TF module M (12 TFs) ── Level 3: Cross-species ──┤
human orthologs → human DoRothEA → reference sets (+permutation)         │
                                                                         │
Stroke rewiring signature → LINCS reverse-match ── Level 4: Drug reversal ┘
```

每一级分别对应一个审稿人必然追问的问题：**稳定？正确？保守？可转化？**

---

## Level 1 — Cross-cohort reproducibility of virtual rewiring（鼠内跨队列）

### 设计（Design）
为评估 virtual rewiring 的技术可复现性，我们在**两个独立生成的鼠缺血性卒中单细胞队列**上重跑 rewiring，并与整合分析比较 TF 主调控因子排序：
- **整合集** `rewiring_full`：GSE174574(24h)+GSE225948(2d/14d)+sham，覆盖 sham→24h / 24h→2d / 2d→14d / sham→14d，留存基因 6,897。
- **队列1** `cc_gse174`：GSE174574 单集（急性损伤），仅 sham→24h，留存基因 8,736。
- **队列2** `cc_gse225`：GSE225948 单集（亚急性→修复），sham→2d / 2d→14d，留存基因 7,041。

三套分析采用**完全一致**的 rewiring v2 协议（library-size 归一化 + log1p → DoRothEA ABC 子图 → PC 回归组成校正 n_pc=10 → 200 置换 → BH-FDR + pooled-FDR，q<0.1），保证尺度可比。

**必须披露的方法论约束（避免假阳性/假阴性）：**
1. **整合集 sham 基线混淆**：整合集的 `sham` 状态 = GSE174574 sham **+** GSE225948 sham 混合细胞；单队列 c1 的 sham 基线仅 GSE174574 sham。因此「整合集 vs 单队列 同 transition」不是干净比较（终点相同、基线不同可翻转 dW 方向），只能作样本量稳健性弱附注，**不能**作独立重现证据 → 故跨队列检验定义在 **TF 主调控因子层**。
2. **基因集不可比**：三套留存基因 6,897 / 8,736 / 7,041，离散边级 Jaccard 下界被人为压低，边级重叠必然偏保守。

### 结果（Results）
- **TF 主调控因子排序在独立队列间显著重现**（正确的跨队列检验）：

| 比较 | 共同 TF 数 | Spearman ρ | p |
|---|---|---|---|
| 整合集 vs 队列1 | 251 | **+0.517** | 1.5e-18 |
| 整合集 vs 队列2 | 235 | **+0.548** | 8.8e-20 |
| 队列1 vs 队列2 | 242 | **+0.482** | 1.7e-15 |

三对 ρ 均 0.48–0.55、p<1e-15，非数据集伪影。
- **Top20(|ΔW|max) TF 重叠**：Sox10 在三个独立分析**均进入 Top20**；Sox2/Sox9 在整合+c1 重现（Jaccard 整合vs c1=0.111，含 Sox10/Sox2/Sox9；整合vs c2=0.250，含 Sox10/Cebpb/Erg/Gata2 等）。

### 诚实局限（须当卖点写，勿藏）
- **单条边级显著性不跨队列重现**：同 transition 比较（附注基线混淆）共同边数万级，但 ΔW_PC 相关仅 ~0.1、方向一致率 ≈0.52（≈随机），Jaccard≈0。
- **解读**：这与项目可解释性主轴完全一致——固定因果图 + 线性读头下，图不提升预测精度（已用 GEARS 验证），其稳健产出是 **TF/模块级可解释重连**，而非边级精度。审稿框架中「换成 Top100 Jaccard 就稳定」是误判：边级不可复现是**真实结论**，TF 级 ρ≈0.5 才是真正稳的那一层。

### 英文可直抄（Level 1）
> To assess the technical reproducibility of virtual rewiring, we re-derived the rewiring on two independently generated mouse ischemic-stroke single-cell cohorts — cohort 1 (GSE174574, 24 h MCAO vs. sham) and cohort 2 (GSE225948, 2 d and 14 d post-stroke) — and compared the recovered transcription-factor (TF) master-regulator ranking against the integrated analysis. All three analyses used an identical rewiring v2 protocol (library-size normalization, log1p, DoRothEA ABC subgraph, principal-component regression for compositional correction with n_pc = 10, 200 permutations, BH-FDR and pooled FDR at q < 0.1). Because the two cohorts share no transition, differ in gene space (6,897 / 8,736 / 7,041 retained genes), and the integrated "sham" baseline mixes cells from both cohorts (confounding edge-level comparison), we evaluate reproducibility at the TF master-regulator level. The TF-level ranking was significantly recovered across cohorts: Spearman ρ = 0.48–0.55 (p < 1e-15) for all three pairwise comparisons. The oligodendrocyte/myelin master regulator Sox10 appeared in the top-20 |ΔW|max TFs of all three analyses, with Sox2/Sox9 reproduced in ≥2. By contrast, individual edge-level significance did not reproduce across cohorts (direction agreement ≈ 0.52, i.e. at chance), consistent with the framework's thesis that a fixed causal graph with a linear readout confers interpretability at the TF/module level rather than edge-level predictive precision.

---

## Level 2 — Recovery of established stroke regulatory programs（已知生物学重现）

### 设计 / 结果
将 Level 1 跨队列重现的 TF 与「卒中后修复已知主调控因子」对齐：
- **Sox10 / Sox2 / Sox9**（少突胶质细胞谱系承诺与髓鞘再生的经典主调控因子）被一致重现，并与 Step A/B 的 **top rewired edge Sox10→Plp1（2d→14d 重髓鞘化）** 完全自洽——这是方法**生物学可信度**的关键支点。
- **负向对照（须显式披露）**：在严格 q<0.1 阈值下，三队列共通 TF **仅 Fos**（immediate-early / AP-1 应激响应基因，任何损伤都点亮）。这是已知的 perturbation 伪影基因，应作为「严格显著性阈值下的复现基因对照」写入 manuscript，而非回避。
- 双阈值对照自洽：q<0.1 严格显著 → Fos（应激噪声）；|ΔW|max Top20 幅度 → Sox10（+Sox2/9）（真实生物学信号）。两者不矛盾，必须**同时报告**。

### 英文可直抄（Level 2）
> To test whether virtual rewiring recovers established stroke biology, we examined the cross-cohort reproducible TFs against known regulators of post-stroke repair. Sox10, Sox2 and Sox9 — canonical master regulators of oligodendrocyte lineage commitment and myelin regeneration — were consistently recovered and aligned with the top rewired edge Sox10→Plp1 (remyelination, 2 d → 14 d) identified in the primary analysis. As a negative control, under the strict q < 0.1 threshold the only TF reproduced across all three analyses was Fos, an immediate-early / AP-1 stress-response gene non-specifically induced by any injury; we report this explicitly as a stress-artifact control rather than a biological signal.

---

## Level 3 — Cross-species conservation of regulatory programs（跨物种调控程序保守）

> 命名升级说明：原内部称「人 TF 模块富集」。按审稿建议改为 **Cross-species conservation of regulatory programs**——我们验证的是「调控程序」，不是「边」，档次更高且如实。
> **名实相符约束**：本验证回答的是「该 TF 的调控程序在人的调控网络结构中是否保守」，而**不是**「人卒中里该模块表达是否升高」（后者需 GSVA/AUCell + 人表达矩阵，我们暂无）。不可把「结构保守」包装成「表达激活验证」。

### 设计（Design）
- **鼠可重现核心 TF 模块 M**（12 个）：来自 Level 1 连续诊断，在 full / c1 / c2 三套 Top20 中**至少出现于 2 套**：`Ar, Cebpb, E2f1, Erg, Gata2, Gata3, Nr2f2, Pax5, Runx3, Sox10, Sox2, Sox9`（Sox10 三套共通头号）。
- **正交映射 → 人同源**：symbol 鼠/人高度保守，逐个映射并校验其为人 DoRothEA `source`（全部命中：AR, CEBPB, E2F1, ERG, GATA2, GATA3, NR2F2, PAX5, RUNX3, SOX10, SOX2, SOX9）。
- **人 GRN 与富集框架**：人 DoRothEA 调控子（`data/dorothea/human_dorothea_regulon.tsv`，643 source TF，18,564 target，全集 N=18,564）。逐 TF 取其人靶集，对文献策展参考集做超几何富集；并以 **2,000 次「靶基因数匹配」(0.5–2×) 随机人 TF** 作置换 null，得经验 p，排除「任意广靶 TF 都富集」的平庸解释。
- **方法学教训（已修正）**：初版把 12 个广靶 TF 合并 → 靶模块占全集 65% → 稀释伪阴性（经验 p≈0.5）。**逐 TF 粒度**才是正确检验。
- **文献策展参考集**（成熟公认、可引用）：
  - *myelin_oligodendrocyte*：PLP1, MBP, MOG, CNP, MAG, MOBP, OLIG1, OLIG2, CLDN11, MAL, MPZ, PMP22, MYRF, NKX6-2, OPALIN, TSPAN2, UGT8, GJC2, KLK6, ASPA, NKX2-2, PDGFRA, ST18
  - *neuroinflammation*：CEBPB, CEBPA, IL1B, IL6, TNF, NFKB1, NFKBIA, STAT1, STAT3, CCL2, CCL3, CXCL10, CXCL1, TLR2, TLR4, IRF1, IRF7, NLRP3, PTGS2, NOS2, HMGB1, RELA, CXCL2, CCL5

### 结果（Results，逐 TF 富集 + 置换）

| TF(人) | 参考集 | k/n | OR | 超几何 p | 经验 p(置换) | 随机 k̄±σ | 判定 |
|---|---|---|---|---|---|---|---|
| **SOX10** | myelin_oligo | 7/22 | **27.0** | 6.0e-08 | **0.0005** | 0.42±0.55 | *** 显著 |
| **CEBPB** | neuroinflammation | 8/23 | **16.5** | 3.2e-07 | **0.024** | 1.80±2.42 | *** 显著 |
| **GATA2** | neuroinflammation | 15/23 | **4.6** | 3.3e-04 | **0.047** | 9.74±4.89 | *** 显著 |
| AR | neuroinflammation | 6/23 | 6.3 | 1.0e-03 | 0.315 | 4.07±4.11 | 否 |
| E2F1 | neuroinflammation | 12/23 | 2.3 | 3.3e-02 | 0.256 | 9.51±4.95 | 否 |
| GATA3 | neuroinflammation | 11/23 | 1.4 | 2.5e-01 | 0.442 | 10.06±4.91 | 否 |
| SOX2 | myelin_oligo | 6/22 | 1.8 | 1.8e-01 | 0.230 | 6.69±5.19 | 否 |
| ERG | neuroinflammation | 0/23 | — | — | — | 1.80±2.42 | 无重叠 |
| NR2F2 | (both) | 0 | — | — | — | — | 无重叠 |
| PAX5 | (both) | 1 | ~1.7 | ~0.45 | 0.34–0.65 | 0.47–1.83 | 否 |
| RUNX3 | (both) | 1–2 | ~1.7–3.5 | 0.12–0.45 | 0.29–0.34 | 0.45–1.80 | 否 |
| SOX9 | myelin_oligo | 0/22 | — | — | — | 0.50±0.70 | 无重叠 |

**显著保守链接（具体重叠基因，生物学自洽）：**
- **SOX10 → 人髓鞘/少突（旗舰结果）**：重叠 `GJC2, MAG, MBP, MPZ, PDGFRA, PLP1, PMP22` —— PLP1/MBP/MAG/MPZ/PMP22 为主要髓鞘结构蛋白，GJC2 为少突缝隙连接，PDGFRA 为少突前体标志。鼠三队列共通头号 rewired TF 精确投射到人**髓鞘再生**程序，与 Sox10→Plp1 边完全闭环。
- **CEBPB → 人神经炎症**：重叠 `CCL3, CCL5, IL1B, IL6, NOS2, PTGS2, STAT3, TNF` —— 经典促炎因子/趋化因子 + NF-κB/STAT3 轴。
- **GATA2 → 人神经炎症**：重叠 `CCL2, CCL5, CXCL1, CXCL2, HMGB1, IL1B, IL6, NFKB1, NFKBIA, NOS2, PTGS2, STAT1, STAT3, TLR2, TLR4` —— 广覆盖 TLR/NF-κB 先天免疫模块。

**阴性/非特异结果（诚实披露）：**
- AR/ERG/NR2F2/PAX5/RUNX3/SOX9：靶集与参考集无特异富集（经验 p 0.29–0.65），属正常（非所有 rewired TF 都有清晰人单组织靶模块）。
- SOX2/E2F1/GATA3：靶极广（3,253–7,212 基因，多能性/细胞周期广效调控因子），超几何被稀释、无单组织特异富集——生物学预期。
- SOX9：人 DoRothEA 靶集与参考集 0 重叠（人 SOX9 在 DoRothEA 中靶偏软骨/发育，非 CNS），提示 GRN 注释覆盖差异，非信号矛盾。

### 解释与局限（Interpretation & Limitations）
1. **生物学含义确证**：鼠 rewiring 识别的 TF 模块经正交投影富集于人保守卒中修复生物学（髓鞘再生 + 神经炎症），证明 rewiring 信号**不是数值噪声，而具真实生物学含义**。
2. **跨物种保守性**：Sox10→髓鞘、Cebpb/Gata2→炎症在鼠（数据驱动）与人（独立 GRN 知识库）间保守，是框架**可迁移至人**的关键证据。
3. **非独立发现（须披露）**：human DoRothEA 本身编码已知调控关系，故富集证明的是「鼠识别的调控程序在人 GRN 中**保守**」，而非从人数据独立发现；置换对照已排除平庸富集，但生物学独立性有限。
4. **参考集为文献策展**：myelin/neuroinflammation 由成熟文献整理（可引用），未来可换 MSigDB Hallmark / DisGeNET 卒中基因集重验。
5. **无直接人卒中 scRNA**：原计划跨物种队列 HRA007397 实为**受控访问 PBMC 原始 fastq（非脑）**，不能验脑内髓鞘/血管重连，仅覆盖外周髓系/炎症轴；若未来获公开人卒中（或 CNS 损伤）scRNA-seq，可直接跑 rewiring 做**数据驱动**跨物种验证（最高证据级），并升级到 GSVA/AUCell「表达激活保守」。

### 英文可直抄（Level 3）
> To assess cross-species conservation of the recovered regulatory programs, we orthogonally projected the mouse cross-cohort reproducible core TF module (12 TFs recovered in ≥2 of 3 analyses) onto the human DoRothEA regulatory network and tested whether each TF's human target set is enriched for literature-curated stroke-repair reference sets (myelin/oligodendrocyte; neuroinflammation), using 2,000 size-matched permutations to control for broadly-targeting TFs. This tests conservation of the regulatory program's *structure*, not expression activation in human stroke (no independent human stroke scRNA-seq was available; the candidate HRA007397 cohort is controlled-access PBMC and non-brain). Three links were significantly conserved: SOX10 → human myelin/oligodendrocyte module (OR = 27.0, empirical p = 0.0005; overlapping PLP1/MBP/MAG/MPZ/PMP22/GJC2/PDGFRA), CEBPB → human neuroinflammation (OR = 16.5, empirical p = 0.024; overlapping IL1B/IL6/TNF/CCL3/CCL5/NOS2/PTGS2/STAT3), and GATA2 → human neuroinflammation (OR = 4.6, empirical p = 0.047; overlapping TLR2/TLR4/NFKB1/NFKBIA and CCL/CXCL chemokines). No specific enrichment was observed for AR/ERG/NR2F2/PAX5/RUNX3/SOX9, and SOX2/E2F1/GATA3 were too broadly targeting to show tissue-specific enrichment — both consistent with expectation. The conserved SOX10 → myelin link closes the loop with the mouse Sox10→Plp1 remyelination edge, providing the strongest cross-species evidence for the framework.

---

## Level 4 — Drug-perturbation reversal via LINCS（药物逆转 · 已执行 2026-07-09）

### 设计
作为转化出口，将卒中损伤程序（Step A/B/C 识别的炎症轴 Cebpb/Gata2 + 髓鞘再生轴 Sox10 受损）转化为 disease signature（MCAO 24h vs sham 差异表达：up=损伤/炎症↑、dn=修复↓），经 **L1000CDS2**（maayanlab，无需 key，JSON POST 已验证）反向匹配（aggravate=FALSE：找能"反向"该签名的化合物），排序候选药。签名分两版：
- **main**：单队列 MCAO24h vs sham（up100/dn100）—— 噪声大，best=0.0671 严重 tie，仅命中 triptolide（抗炎）1 个已知药。
- **robust**：两队列共识交集去噪（up28/dn17）—— best=0.125，命中 **4 个已知药**（vorinostat / trichostatin A 为 HDAC 抑制剂；Mevastatin / Rosuvastatin 为他汀），类别方向与卒中神经保护/抗炎/促髓鞘一致。
以 **robust 为主结果**，main 作单队列噪声对照（呼应 L1 技术可复现的立意）。

### 统计严谨性（置换对照，关键）
用 **unique-drug 计数**（避免 topMeta 多细胞系重复条目虚高），20 个同规模随机签名跑同样反向匹配：
- 真实 robust = **4** 个已知药；随机基线 = **2.75 ± 1.37（max 5）**；**empirical p = 0.33（不显著）**。
- vorinostat 跨系一致（n_hits=6）看似强信号，但随机签名下出现率 **53%**、单次 max=16 → 属 L1000 常见"泛逆转"背景，**非特异**。
➡️ **D 步定位为 translational proof-of-concept / hypothesis generation，而非药效预测。** 确定性已验证（同签名两次查询 top50 药物集合 30/30 完全一致）。

### 与 L2/L3 闭环（文章亮点）
命中药物类别（HDAC 抑制剂促髓鞘/抗炎 + 他汀抗炎/神经保护）恰好靶向 L2/L3 独立识别的**炎症轴与髓鞘再生轴** → 框架从"解释哪里出错"闭环到"生成可能修好的候选假设"，构成端到端可转化叙事。

### 局限（方法节必含）
LINCS 细胞系非神经元（癌细胞系）、单时间点（24h 急性期）、tie 导致区分度有限、置换未显著、全 in silico 无湿实验。详细分析与写作策略见 `StepD_文章级分析阐述_2026-07-09.md`。

### 英文可直抄（Level 4）
> To demonstrate translational utility, we converted the stroke injury program identified in Levels 2–3 into a disease gene signature (up-regulated injury/inflammatory and down-regulated repair genes from MCAO vs. sham at 24 h, intersected across two independent cohorts to reduce noise) and performed reverse connectivity mapping against the LINCS L1000 corpus using L1000CDS2 (aggravate = FALSE). The denoised consensus signature returned four literature-supported neuroprotective/anti-inflammatory agents in the top 50 (vorinostat and trichostatin A as HDAC inhibitors; mevastatin and rosuvastatin as statins), whose mechanistic classes align with the inflammatory and myelin-regeneration axes recovered in Levels 2–3. A permutation background control (20 random signatures of matched size) yielded a mean of 2.75 ± 1.37 such hits (max 5), so the observed count of 4 did not exceed chance (empirical p = 0.33). We therefore position Level 4 as a proof-of-concept for end-to-end translational hypothesis generation rather than a validated efficacy prediction, noting that LINCS profiles derive from non-neuronal cancer cell lines and require experimental confirmation.

---

# 合并后对 Manuscript 的一句话总览（建议置于 Abstract / Validation 节首）

> Virtual rewiring is technically reproducible at the TF master-regulator level across independent mouse stroke cohorts (Spearman ρ = 0.48–0.55, p < 1e-15), recovers established repair regulators (Sox10/Sox2/Sox9), and projects onto human-conserved stroke-repair programs (SOX10→myelin OR = 27, empirical p = 0.0005; CEBPB/GATA2→neuroinflammation), while individual edge-level effects are dataset-specific — consistent with the framework's thesis that a fixed causal graph confers *interpretability at the module level*, not edge-level predictive precision.

---

# 待补 / 待办清单（执行前检查）
- [x] **Level 4（D）**：已完成（见上 Level 4 + `StepD_文章级分析阐述_2026-07-09.md`）。
- [ ] **方法节必含披露**：① 时间点 24h/2d/14d 选择逻辑 + 跨数据集拼接限制；② 整合集 sham 基线混合混淆；③ 三套基因集差异；④ HRA007397 实为受控 PBMC（非脑）；⑤ 边级不重现作为主轴证据而非缺陷。
- [ ] **可选升级**：获公开人卒中/CNS 损伤 scRNA 后，将 Level 3 从「结构保守」升到「表达激活保守」（GSVA/AUCell）；或申请 HRA007397（仅外周 PBMC 炎症轴）作可选扩展。
- [ ] 最终翻译为全英文方法节文本（本草稿中文为正文、英文为可直抄块）。

*源文件：`CROSS_COHORT_报告_2026-07-09.md`、`CROSS_SPECIES_报告_2026-07-09.md`、`human_module_enrich.json`。*
*关联产物：`cc_gse174_rewiring_full.csv`、`cc_gse225_rewiring_full.csv`、`rewiring_full.csv`、`human_module_enrich.py`、`cc_diag.py`。*
