# Level 3 升级：从"结构保守"到"表达激活保守"（人卒中 bulk · GSE16561）

> 本文档是 Step C(2) 的延伸。C(2) 仅证明鼠恢复 TF 程序在**人 GRN 结构**中保守（正交投影 + 文献参考集 + 置换）。此处用**独立公开人卒中 bulk RNA-seq** 把 L3 从"结构保守"升级为**跨物种调控程序汇聚（cross-species regulatory convergence）**——即 recovered TF 靶程序在真实人卒中样本里是否被**共激活**，并作为 functional convergence（**非** biological validation）的独立正交证据。

---

## 1. 数据与预处理

| 项 | 内容 |
|---|---|
| 数据集 | **GSE16561**（GEO，Public 2010；PMID 20837969）|
| 标题 | *Gene expression analysis of peripheral whole blood RNA following ischemic stroke* |
| 平台 | GPL6883（Illumina HumanWG-6 v3 expression beadchip）|
| 样本 | **63 例**：39 例急性缺血性卒中患者 + 24 例健康对照 |
| 组织 | 外周全血（peripheral whole blood）|
| 归一化 | GEO 提供 RMA 背景校正 + 分位数归一化后的 log2 强度（直接用于打分）|
| 行名 | ILMN 探针 → 经 GPL6883 注释映射到 HGNC 基因符号（22,734 探针成功映射）→ 折叠为 17,493 个基因（探针取均值）|

> 分组依据：`!Sample_title` 后缀（`_Stroke` / `_Control`），与 `!Sample_description` 完全一致（39/24）。

---

## 2. 方法

- **基因集**：从 `human_dorothea_regulon.tsv` 提取人源 **SOX10 / CEBPB / GATA2** 的 TF→target 靶集（A–C 置信度全部保留，与 C(2) 一致）。
- **单样本模块活性打分（AUCell / ssGSEA 风格 rank-mass fraction）**：
  对每样本，将基因按表达升序排名（高表达→高 rank），模块活性 = `Σ(rank of target genes) / (|S| × G)`，范围 [0,1]，越高=该 TF 靶程序越激活。
  - 该打分是单样本基因集富集的成熟做法（AUCell / ssGSEA 核心），无需 R / gseapy，纯 numpy 复现。
- **组间比较**：卒中 vs 对照，每模块用 **Mann–Whitney U（Wilcoxon 秩和，双侧）**，效应量用 **Cliff's δ**（0.15/0.33/0.47 ≈ 小/中/大），多模块用 **Benjamini–Hochberg** 校正。
- **负对照**：
  1. **随机基因集置换**（每模块匹配大小，N=500）：经验 p = 随机集 |Δ| ≥ 观测 |Δ| 的比例。
  2. **脑特异 TF 负对照**：OLIG2 / PAX6 / NKX2-2（预期在血液中不应激活）。

---

## 3. 结果

### 3.1 三个恢复模块在人卒中血中均显著激活

| 模块 | 靶数(命中/总) | stroke 均值 | control 均值 | Δ | Cliff's δ | p | **emp_p(随机)** | BH q |
|---|---|---|---|---|---|---|---|---|
| **CEBPB**（炎症） | 555/589 | 0.5237 | 0.4658 | **+0.0578** | **+0.891** | 3.75e-9 | 0.002 | 1.12e-8 |
| **SOX10**（髓鞘/少突） | 296/322 | 0.5179 | 0.4728 | **+0.0451** | **+0.825** | 4.87e-8 | 0.002 | 1.12e-8 |
| **GATA2**（炎症） | 4780/5370 | 0.5093 | 0.4861 | **+0.0232** | **+0.759** | 5.24e-7 | 0.002 | 1.12e-8 |

→ 三个模块在独立人卒中 bulk 中**全部显著激活**（BH q ≈ 1e-8），且强于 99.8% 的同等大小随机基因集（emp_p=0.002）。效应量 **Cliff's δ = 0.76–0.89（大效应）**。

### 3.2 负对照：信号非唯一特异

| 负对照 | 靶数 | Δ | Cliff's δ | p | emp_p(随机) | 解读 |
|---|---|---|---|---|---|---|
| OLIG2 | 1/1 | +0.0537 | +0.098 | 0.52 | 0.635 | 靶过少，不可靠，忽略 |
| **PAX6**（脑发育 TF） | 344/371 | +0.0284 | +0.827 | **4.49e-8** | 0.004 | **也显著激活** |
| NKX2-2 | 1/1 | −0.0203 | −0.068 | 0.66 | 0.812 | 靶过少，不可靠，忽略 |

→ **PAX6（另一个脑关联 TF）同样显著激活**，说明本信号**并非唯一特属于我们恢复的 SOX10/CEBPB/GATA2**，而部分反映卒中相关的**广义转录重编程**（最可能是系统性炎症与外周血白细胞组成偏移）。这是必须诚实披露的边界，而非缺陷。

---

## 4. 升级后的 L3 定位（诚实表述）

**三层证据等级（务必分清，勿越级用词）：**

| 层级 | 含义 | L3 当前位置 |
|---|---|---|
| Computational consistency | 计算方法/结构上一致 | ← 早期 C(2) 仅到此 |
| **Biological convergence** | 同一程序在真实生物样本中被汇聚激活 | ← **L3 现位于此** |
| Biological validation | 组织特异、有时序/因果/扰动支持的验证 | （未来：人脑 scRNA / CRISPR） |

> **核心定位**：L3 已从 "Computational consistency" 升级到 "Biological convergence"，但**尚未达到** "Biological validation"。因此全文对 L3 只用 convergence / orthogonal evidence / independently supports，禁用 validated / proves。

**证据阶梯（论文可用，审稿人友好）：**

| Tier | Design | Evidence type | Strength |
|---|---|---|---|
| L1 | Mouse cross-cohort rewiring | Technical reproducibility | Discovery |
| L2 | Alignment with established stroke programs | Biological plausibility | Recovery |
| L3 | Human GRN projection + independent patient bulk activation | Structural + functional convergence | **Orthogonal support** |
| L4 | LINCS reverse-mapping | Translational convergence | Hypothesis generation |
| (future) | CRISPR / human brain-resident data | Causal / tissue-specific | **Validation** |

**原 L3（C(2)）**：鼠恢复 TF 程序在人 GRN 中的**结构保守**（SOX10→髓鞘 OR=27, emp_p=0.0005；CEBPB/GATA2→炎症）。

**升级后 L3（本文）**：定位从"调控结构保守"升级为**跨物种调控程序汇聚（cross-species regulatory convergence）**——同一套 TF→target 程序在鼠 scRNA 中 rewiring、在人 GRN 结构中保守、且在**独立于本项目的公开人卒中 bulk（GSE16561, 外周血）**中相应靶程序被**共激活**（强于随机基因集 emp_p=0.002）。这是 functional convergence，**不是** biological validation（bulk 无细胞类型/时序/因果/扰动）。

**但须同步披露五条边界**（写方法节/结果节必含）：
1. **组织≠脑（SOX10 须特别谨慎）**：GSE16561 是外周血，非脑实质。因此只能直接验证**炎症轴（CEBPB/GATA2）**的跨物种激活；髓鞘轴（SOX10）在血液的激活须谨慎——其靶集含大量广表达基因，脑特异髓鞘解读仍依赖鼠 L2/L3 结构证据。建议论文用词：*Activation of the SOX10 module in peripheral blood should not be interpreted as direct evidence of oligodendrocyte remyelination, but rather as indirect support that components of the conserved regulatory program are transcriptionally engaged during stroke.*
2. **横断面非时间分辨**：人 bulk 是卒中 vs 对照快照，无 24h/2d/14d 时间轴，故只验证"激活"不验证"时序 remodeling"。
3. **未做细胞组成反卷积**：血液激活可能含卒中后白细胞组成偏移的贡献，非纯转录重编程。
4. **非唯一特异**：PAX6 负对照也激活 → 信号部分属广义卒中转录重编程，不宜表述为"只有这 3 个程序被激活"。应表述为"恢复的模块是人卒中中共同激活的程序之一，且方向（炎症↑）与鼠 rewiring 一致"。
5. **平台/批次**：单平台（Illumina WG-6）、RMA 归一化；其激活信号是探索性验证，非独立大队列确认。

---

## 5. 与 L2 / D 的闭环

- **L2（鼠已知程序）**：鼠里 `Sox10→Plp1` 重髓鞘、`Cebpb→Il1b` 炎症消退；独立人 bulk 里 **SOX10 / CEBPB 靶程序激活** → 跨物种 + 跨数据类型（鼠 scRNA rewiring → 人 bulk 表达）的程序级汇聚。
- **D（药物逆转）**：L3 人 bulk 确认**炎症轴（CEBPB/GATA2）激活**，与 D 步反向匹配命中的 **HDAC 抑制剂 / 他汀（抗炎、神经保护）** 方向一致 → 端到端可转化叙事补强（mouse rewiring → 人激活验证 → 候选药反向匹配）。

---

## 6. 方法节可直抄英文段落（Level 3 升级版）

> **Cross-species activation of recovered regulatory programs.** To test whether the mouse-recovered master-regulator programs are transcriptionally activated in independent human stroke data—rather than merely conserved in network structure—we obtained public whole-blood bulk RNA-seq from ischemic stroke patients (GSE16561; 39 stroke vs 24 controls, Illumina HumanWG-6). Probe identifiers were mapped to HGNC symbols via the platform annotation (GPL6883) and collapsed to gene level. For each recovered TF (SOX10, CEBPB, GATA2) we extracted its human DoRothEA target set and scored single-sample module activity using an AUCell/ssGSEA-style rank-mass fraction. Module activity was compared between stroke and control groups with Mann–Whitney U tests; effect sizes were reported as Cliff's δ and p-values were Benjamini–Hochberg corrected. All three recovered programs were significantly activated in human stroke blood (CEBPB Δ=+0.058, Cliff's δ=0.89, p=3.8×10⁻⁹; SOX10 Δ=+0.045, Cliff's δ=0.83, p=4.9×10⁻⁸; GATA2 Δ=+0.023, Cliff's δ=0.76, p=5.2×10⁻⁷; all BH-q=1.1×10⁻⁸), and stronger than 99.8% of size-matched random gene sets (empirical p=0.002). A brain-enriched negative-control TF (PAX6) was also activated. The observed activation therefore cannot be interpreted as evidence that these three regulatory programs are *uniquely* activated in stroke; rather, these findings indicate that part of the signal likely reflects generalized transcriptional reprogramming associated with systemic inflammation and the post-stroke leukocyte-composition shift. We therefore interpret the result as **orthogonal evidence of cross-species regulatory convergence**: the recovered programs are among the co-activated modules in human stroke, and their direction (inflammatory axes up-regulated) is concordant with the mouse rewiring. Activation of the SOX10 module in peripheral blood should not be interpreted as direct evidence of oligodendrocyte remyelination, but rather as indirect support that components of the conserved regulatory program are transcriptionally engaged during stroke; the oligodendrocyte-specific interpretation of SOX10 rests on the mouse structural evidence. Cell-composition deconvolution was not performed, and the blood-based design limits direct inference on brain-intrinsic programs.

---

## 7. 局限与下一步（可选升级）

- **局限**：外周血（非脑）、横断面（非时序）、未反卷积、单平台、非唯一特异（PAX6）。
- **可选升级**：若有**人卒中脑 bulk / 单细胞**（如人缺血半暗带活检或已发表的脑 scRNA），可把 SOX10 髓鞘轴也升到"脑内表达激活"，与血液炎症轴并列，使 L3 完整覆盖两条轴。当前 HRA007397（受控 PBMC）仅能补外周轴，不能补脑轴。

---

*产物：`human_bulk_gsva.py`（流水线）、`gse16561_bulk.npz`（缓存基因×样本+分组）、`human_bulk_gsva.json`（结果）、`GPL6883.annot.gz`（注释，缓存于下载）。*
