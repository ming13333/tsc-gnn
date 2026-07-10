# Manuscript Draft — Methods & Results（统一初稿 · 2026-07-09）

> **用途**：A→B→C→D 全部分析完成后统一撰写的方法节 + 结果节初稿。可直接进入 manuscript 的「Methods / Results」骨架。
> **作者红线（务必保留，勿为迎合叙事软化）**：
> 1. 图（固定因果 GRN + 线性读头）**不提升预测精度**——已用 GEARS 真实 Perturb-seq 与 90 点鲁棒性研究证实；主轴必须是**可解释性**（边级 ΔW + 模块级 master-regulator），不是预测精度。
> 2. 时间轴 24h/2d/14d 由 **跨数据集拼接**（GSE174574 + GSE225948），非均匀采样，存在 batch×time 混淆；须在方法/图注显式披露。
> 3. 边级显著性**不**跨队列重现（Jaccard≈0，方向一致≈0.52）——这是主轴卖点，须如实写、当证据而非缺陷藏。
> 4. Level 3 = 调控程序在**人 GRN 结构保守**，非「人卒中表达激活」（无独立人 scRNA 表达矩阵；HRA007397 是受控访问 PBMC 非脑）。
> 5. Level 4 药物逆转 **置换未显著（emp_p=0.33）**，定位为 translational proof-of-concept，非药效预测；vorinostat 跨系一致是 L1000 背景而非特异信号。
> 6. 整合集 sham 基线 = 两队列混合；三套基因集（6,897 / 8,736 / 7,041）不可直接比边级。

---

# Methods

## 2.1 Conceptual framework: temporal-state-conditioned virtual perturbation and the interpretability thesis

We frame *in-silico perturbation* as an operation on a latent graph space rather than on a single-gene knockout. Our framework, **TSC-GNN** (Temporal and cell-State-Conditioned Graph Neural Network), takes a multi-timepoint, multi-condition single-cell transcriptome and (i) builds a time- and state-conditioned gene-regulatory graph, (ii) propagates a perturbation along this graph to obtain a rewired latent embedding, and (iii) reads out *edge-level rewiring* — the change in transcription-factor (TF)→target coupling across a transition — together with its permutation significance and the resulting master-regulator ranking. Cell-cell communication edges and pseudotime-resolved state vectors condition the propagation; a LINCS/L1000 reverse-match converts the injury embedding into candidate-repurposing hypotheses (Fig. 1, conceptual).

> **作者红线 #1（写进 Introduction/Discussion）**：We do **not** claim that the graph-based readout predicts perturbation responses more accurately than a linear baseline. On the GEARS Perturb-seq benchmark, our graph readout yielded −8.3% (single) and −4.2% (double) improvement with confidence intervals spanning zero; a 90-configuration robustness study (10 seeds × 5 graphs × 2 tasks) found **zero** cases in which the graph beat the linear baseline and 86 with no significant difference. We therefore position TSC-GNN deliberately as an *interpretable* perturbation framework: its value is the **edge-level rewiring ΔW with permutation significance and module-level master-regulator recovery**, which a black-box linear or foundation-model readout does not expose. This positioning is the paper's central thesis, not a weakness to hide.

## 2.2 Data

**Mouse ischemic-stroke single-cell transcriptomes.** We used two independently generated mouse middle-cerebral-artery-occlusion (MCAO) scRNA-seq cohorts:
- **Cohort 1 — GSE174574** (Li et al., 2021): 3 MCAO (24 h) + 3 sham, 57,528 cells, with a BEAM pseudotime trajectory. Single acute time point.
- **Cohort 2 — GSE225948** (Anrather et al., *Nature Immunology* 2024): genuine two-timepoint post-stroke time course (2 d and 14 d, brain and blood).

**Constructed time axis and its disclosure.** No single public scRNA-seq cohort covers all three stages. We therefore assembled a 24 h → 2 d → 14 d axis spanning the acute → sub-acute peak → repair/remodelling phases by stitching cohort 1 (24 h + sham) with cohort 2 (2 d / 14 d). *This is a limitation that must be stated in the Methods and figure legends:* the axis is **non-uniformly sampled** (the 3–7 d window is absent) and carries a **batch×time confound** (the two time segments come from different studies/platforms). We mitigate — but do not eliminate — the confound by treating each time point's within-state coupling independently and by PC-regression composition correction (§2.4); we do **not** interpret discrete ΔW between stitched segments as a continuous *rate*.

**Pre-processing.** Count matrices were loaded offline (no scanpy dependency) and processed with library-size normalization (×(1e4)) followed by log1p, in a fixed conda environment (`bbb_gnn`; numpy 2.2.6 / scipy 1.15.3 / pandas 2.3.3). Cell-type composition was retained per cell for the composition-correction step. For cross-cohort analyses, each cohort was processed under an identical protocol so that the retained gene spaces (6,897 integrated / 8,736 cohort-1 / 7,041 cohort-2) are comparable, with the gene-space discrepancy itself reported as a limitation on edge-level comparison.

## 2.3 Gene-regulatory network construction

We used the **DoRothEA** consensus regulons (confidence levels A–C) as a *directed causal* TF→target graph — direction and sign (activation/repression) are taken from the literature-curated prior, so that any rewiring we report is interpretable as "whose regulation is enhanced/weakened", not merely as undirected co-expression change. Mouse and human regulons were read from local TSV exports (`gate1/data/dorothea/`), making the pipeline fully offline and reproducible (SHA-256 manifest per run). For each TF we additionally computed a **state-affinity** vector `A_aff` (subsample of n=4,000 cells) capturing how strongly the TF's targets are expressed in each cell state; edges in the top 50% by |A_aff| are retained as *state-conditioned* edges for rewiring testing.

> **为什么用 DoRothEA 而非 k-NN 相关图**：两者在预测上表现几乎相同（都不超线性），但 DoRothEA 是有向因果图，rewiring 可解释为 TF→target 的可解释调控变化——这正是「图价值在解释不在预测」的支点。鲁棒性研究已确认结论对图类型（相关 vs 因果）稳健。

## 2.4 State-conditioned edge rewiring

For each tested directed edge *e* = (TF *u* → target *v*) and each time point *t*, we compute the **within-state Pearson coupling**
$$r_{e,t} = \mathrm{corr}\big(x_u^{(t)},\, x_v^{(t)}\big),$$
where $x_u^{(t)}, x_v^{(t)}$ are the log-normalized expressions of *u* and *v* restricted to cells in state/time *t*. The **rewiring effect** for a transition $t_1\!\to\!t_2$ is the coupling change
$$\Delta W_{e,\,t_1\to t_2} = r_{e,t_2} - r_{e,t_1}.$$
Positive ΔW denotes *coupling enhancement* (the TF–target co-expression link strengthens across the transition); negative denotes *coupling weakening*. The DoRothEA prior direction (activation/repression) is reported independently and is not conflated with the empirical coupling change.

**Composition correction.** Because stroke alters cell-type composition across time, raw coupling can reflect composition shift rather than regulatory rewiring. We residualize $x_u, x_v$ against the top **n_pc = 10** principal components of the full expression matrix (least-squares regression of each gene on [1, PCs]), so that ΔW reflects coupling change after removing the major axes of compositional variation. We report both raw and PC-corrected ΔW; PC correction can flip edge direction when composition masked the true regulatory change (e.g. *Sox10→Ank3*: raw ΔW = −0.38 → PC-corrected +0.73).

**Permutation test and multiple testing.** For each transition we permute the time labels (n_perm = 200, seed = 2) and recompute ΔW to build a null per edge; the per-edge two-sided p-value is $p_e = (1 + \#\{\text{perm}: |\Delta W^{\text{null}}_e| \ge |\Delta W^{\text{obs}}_e|\})/(n_{\text{perm}}+1)$ (p_min ≈ 0.005). We apply two corrections: (i) **Benjamini–Hochberg FDR** per transition; and (ii) a **permutation pooled-FDR** that standardizes each edge's |ΔW| by its null standard deviation, pools null z-scores across edges and permutations, and computes q-values without depending on p-resolution — the preferred correction given the limited n_perm. Edges at pooled q < 0.05 (typically 8–36 per transition) are reported as significant rewiring; q < 0.1 is used for broader TF-level recovery.

**Reproducibility.** All runs emit a manifest (command, library versions, seeds, SHA-256 of the cache, resident-matrix size). The full analysis is deterministic under fixed seeds.

## 2.5 Four-level validation framework

Because a fixed-graph method confers interpretability rather than accuracy, we validate it along four escalating levels, each answering one question a reviewer must ask:

| Level | Question | Design |
|---|---|---|
| **L1** Technical reproducibility | Is the method stable? | Re-derive rewiring on two independent mouse cohorts; compare TF master-regulator ranking. |
| **L2** Recovery of established programs | Is it biologically correct? | Align recovered TFs with known post-stroke repair regulators; report stress-artifact controls. |
| **L3** Cross-species conservation | Does it transfer to human? | Orthogonally project the mouse core TF module onto the human DoRothEA network; test enrichment of literature-curated repair reference sets with size-matched permutation. |
| **L4** Drug-perturbation reversal | Is it translatable? | Convert the stroke injury signature into a LINCS L1000 reverse-connectivity map; rank candidate-repurposing agents. |

> **作者红线 #3/#4/#5 已内嵌于各 Level 设计与 Results 表述。**

---

# Results

## 3.1 Virtual rewiring captures time-resolved stroke regulatory rewiring at the module level

Applying state-conditioned rewiring across the four transitions (sham→24 h, 24 h→2 d, 2 d→14 d, sham→14 d) recovered a coherent, literature-aligned program of stroke recovery. Under pooled q < 0.05, 8–36 edges per transition were significant.

- **Acute injury onset (sham→24 h):** enrichment of damage/inflammatory coupling, consistent with the known acute ischaemic response.
- **Repair initiation (24 h→2 d):** a sharp surge in *Sox2/Sox10/Sox9* coupling — the oligodendrocyte-lineage commitment programme that initiates remyelination.
- **Active remyelination (2 d→14 d):** the top rewired edge was **Sox10→Plp1**, the canonical oligodendrocyte myelin structural gene — i.e. the method independently finds the remyelination programme at the repair peak.
- **Inflammation resolution (sham→14 d):** **Cebpb→Il1b** coupling weakens, tracking resolution of the acute inflammatory programme.

**Sanity of recovered targets (Step B).** Targets of significant rewired edges were enriched for oligodendrocyte (OR 36–61), neuron (OR 110) and microglia (OR 15) markers; 14/15 known stroke-relevant TFs appeared, and 77–83% of edges preserved their direction after PC composition correction — confirming the rewiring signal is not a compositional artefact.

> **PC 校正的意义（写进结果）**：raw→PC 方向翻转（如 Sox10→Ank3）证明组成混淆会掩盖真实调控增强，校正后是更可靠的 rewiring 估计。这是方法可靠性的内部证据。

## 3.2 Predictive accuracy is not superior to linear baselines — and this motivates the interpretability axis

Before presenting the biological results, we report an honest scoping result that defines the paper's contribution. We benchmarked the graph readout against a linear baseline under matched inputs (identical perturbation vector *p*, identical training split, identical readout head) across **10 random seeds × 5 graph types** (k-NN, DoRothEA, random, permuted-DoRothEA, 0-hop) **× 2 tasks = 90 configurations**. The graph conferred **no** predictive advantage: it beat the linear baseline in **0** cases and showed no significant difference in **86**. (A non-linear KernelRidge probe overfit, rel ≈ −126%, indicating that any accuracy gain would require a trained deep non-linear readout — out of scope here.) On the external GEARS Perturb-seq benchmark the graph readout gave −8.3% (single) / −4.2% (double) improvement, both with CIs including zero.

**Interpretation.** A fixed causal graph with a linear readout is mathematically contained in the linear hypothesis space, so it cannot *out-predict* linear. We therefore do not sell TSC-GNN as a more accurate predictor. Its deliverable is precisely what the linear/foundation-model baselines omit: **edge-level rewiring ΔW with permutation significance and module-level master-regulator rankings** — the substrate for the interpretability claims in §3.3–3.6. We state this trade-off explicitly rather than obscuring it.

## 3.3 Level 1 — Cross-cohort reproducibility

To test technical stability we re-derived the rewiring on the two independent cohorts and compared TF master-regulator rankings against the integrated analysis. *Crucially*, because the integrated "sham" baseline mixes cells from both cohorts (and the two cohorts share no transition and differ in gene space), we evaluate reproducibility at the **TF master-regulator level**, not at the edge level.

The TF-level ranking reproduced significantly across cohorts:

| Comparison | Common TFs | Spearman ρ | p |
|---|---|---|---|
| Integrated vs cohort 1 | 251 | **+0.517** | 1.5e-18 |
| Integrated vs cohort 2 | 235 | **+0.548** | 8.8e-20 |
| Cohort 1 vs cohort 2 | 242 | **+0.482** | 1.7e-15 |

**Sox10** appeared in the top-20 |ΔW|max TFs of **all three** analyses; **Sox2/Sox9** reproduced in ≥2. By contrast, individual edge-level significance did **not** reproduce (direction agreement ≈ 0.52, i.e. at chance; Jaccard ≈ 0). This is the expected consequence of a fixed-graph/linear-readout method and is reported as *evidence for* the interpretability thesis (module-level signal is robust; edge-level noise is dataset-specific), not as a defect.

> **作者红线 #3**：边级不跨队列重现 = 卖点。审稿框架中「换 Top100 Jaccard 就稳定」是误判——边级不可复现是真实结论，TF 级 ρ≈0.5 才是真正稳的那一层。

## 3.4 Level 2 — Recovery of established stroke programs

The cross-cohort-reproducible TFs aligned with established post-stroke repair biology. **Sox10, Sox2 and Sox9** — canonical master regulators of oligodendrocyte lineage commitment and myelin regeneration — were consistently recovered and closed the loop with the top rewired edge **Sox10→Plp1** (remyelination, 2 d→14 d) from §3.1.

As a negative control we report explicitly that, under the strict q < 0.1 threshold, the *only* TF reproduced across all three analyses was **Fos** — an immediate-early / AP-1 stress-response gene non-specifically induced by any injury. We surface this as a stress-artifact control rather than a biological signal, and we report it alongside the magnitude-based Top-20 recovery (Sox10 etc.) so the two thresholds are not confused.

## 3.5 Level 3 — Cross-species conservation of regulatory programs

> **Naming & scope (作者红线 #4).** We title this "Cross-species conservation of regulatory programs" — we validate *programs*, not edges. The test answers whether the mouse-recovered TF's regulatory program is **conserved in the structure of the human GRN** (via the independent human DoRothEA network), *not* whether the module is up-regulated in human stroke expression (no independent human stroke scRNA-seq was available; the candidate HRA007397 cohort is controlled-access PBMC and non-brain). We do not conflate structural conservation with expression-activation validation.

We took the 12-TF cross-cohort core module (recovered in ≥2 of 3 analyses; **Sox10** common to all three) and orthogonally projected each TF onto its human DoRothEA target set. For each TF we tested enrichment of two literature-curated reference sets — *myelin/oligodendrocyte* and *neuroinflammation* — by hypergeometric test, with **2,000 size-matched (0.5–2× target-count) random human TF permutations** as the null. Testing per-TF (not as a merged broad module) was essential: merging the 12 broadly-targeting TFs gave a target set covering 65% of the network and diluted every signal to non-significance.

Three links were significantly conserved:

| TF (human) | Reference set | k/n | OR | Hypergeom. p | Emp. p (perm) |
|---|---|---|---|---|---|
| **SOX10** | myelin/oligo | 7/22 | **27.0** | 6.0e-08 | **0.0005** |
| **CEBPB** | neuroinflammation | 8/23 | **16.5** | 3.2e-07 | **0.024** |
| **GATA2** | neuroinflammation | 15/23 | **4.6** | 3.3e-04 | **0.047** |

**SOX10 → human myelin/oligodendrocyte (flagship):** overlapping *PLP1, MBP, MAG, MPZ, PMP22* (major myelin structural proteins) plus *GJC2* (oligodendrocyte gap junction) and *PDGFRA* (oligodendrocyte-precursor marker). The mouse top rewired TF projects precisely onto the human myelin-regeneration program, closing the loop with Sox10→Plp1.

**CEBPB → human neuroinflammation:** overlapping *IL1B, IL6, TNF, CCL3, CCL5, NOS2, PTGS2, STAT3*. **GATA2 → neuroinflammation:** overlapping *TLR2, TLR4, NFKB1, NFKBIA* and *CCL/CXCL* chemokines — the TLR/NF-κB innate-immune module.

Non-significant / expected-negative results are reported honestly: AR/ERG/NR2F2/PAX5/RUNX3/SOX9 showed no tissue-specific enrichment (emp. p 0.29–0.65); SOX2/E2F1/GATA3 are too broadly targeting to show single-tissue enrichment (biologically expected). The enrichment uses the *human* DoRothEA knowledge base, so it confirms conservation of the program's structure rather than an independent human discovery; the permutation control excludes the trivial "any broad TF enriches" explanation.

## 3.6 Level 4 — Drug-perturbation reversal via LINCS

As a translational demonstration we converted the stroke injury program (inflammatory axis up, myelin-repair axis down, from §3.4–3.5) into a disease signature — up-regulated injury/inflammatory and down-regulated repair genes from MCAO vs. sham at 24 h — and performed **reverse connectivity mapping** against the LINCS L1000 corpus using L1000CDS2 (aggravate = FALSE), retrieving compounds whose transcriptional response opposes the stroke signature.

We built two signatures: a single-cohort **main** signature (up100/dn100, best score 0.0671 with severe ties, only 1 known agent hit) and a **robust** consensus signature intersected across the two cohorts (up28/dn17, best 0.125). We report **robust** as the primary result; main serves as a single-cohort noise control, echoing the L1 reproducibility theme.

The robust signature returned **four literature-supported neuroprotective/anti-inflammatory agents in the top 50**: **vorinostat** and **trichostatin A** (HDAC inhibitors) and **mevastatin** and **rosuvastatin** (statins) — mechanistic classes that align exactly with the inflammatory and myelin-regeneration axes recovered in Levels 2–3.

**Statistical rigour (permutation background).** Using unique-drug counting (to avoid inflating the count by L1000's multi-cell-line duplicate entries), 20 random signatures of matched size returned a mean of **2.75 ± 1.37** such hits (max 5). The observed count of **4** did **not** exceed chance (**empirical p = 0.33**). A preliminary duplicate-counting scheme had spuriously suggested p = 0.048; we corrected this and do **not** claim significant enrichment. Further, vorinostat's cross-cell-line consistency (n_hits = 6) proved to be an L1000 background effect: under random signatures vorinostat appeared in 53% of runs (max n_hits = 16).

> **作者红线 #5（必写进 Limitation）**：Level 4 is a **translational proof-of-concept / hypothesis-generation** demonstration, *not* a drug-efficacy prediction. LINCS profiles derive from non-neuronal cancer cell lines; the signature is single-timepoint (24 h acute); L1000CDS2 shows heavy ties; the permutation did not reach significance; and the entire pipeline is in silico with no wet-lab confirmation. We state these limits explicitly.

**Narrative closure (the paper's highlight).** The candidate classes surfaced by reverse mapping — HDAC inhibitors (pro-myelin / anti-inflammatory) and statins (anti-inflammatory / neuroprotective) — target *precisely the programs the framework independently recovered as dysregulated in stroke* (Levels 2–3). This closes the loop from interpretable virtual perturbation to a testable therapeutic hypothesis, which is the end-to-end translational story current GRN / virtual-perturbation methods papers largely lack.

---

# Abstract one-liner (suggested)

> Virtual rewiring is technically reproducible at the TF master-regulator level across independent mouse stroke cohorts (Spearman ρ = 0.48–0.55, p < 1e-15), recovers established repair regulators (Sox10/Sox2/Sox9), and projects onto human-conserved stroke-repair programs (SOX10→myelin OR = 27, empirical p = 0.0005; CEBPB/GATA2→neuroinflammation), while individual edge-level effects are dataset-specific — consistent with the framework's thesis that a fixed causal graph confers *interpretability at the module level*, not edge-level predictive precision; the recovered program further reverse-maps onto known neuroprotective drug classes, demonstrating end-to-end translational hypothesis generation.

---

# 待办 / 诚信核对清单（提交前逐条确认）

- [ ] Methods 必含：①时间轴 24h/2d/14d 跨数据集拼接 + 非均匀采样 + batch×time 混淆披露；②DoRothEA 有向因果图 + state-affinity 边筛选；③PC 组成校正(n_pc=10) + 置换(n_perm=200) + pooled-FDR；④图不优于线性的鲁棒性结论（§3.2）。
- [ ] Results 必含：①主 rewiring 四转移拓扑（Sox10→Plp1 等）；②L1 TF rank ρ 表 + 边级不重现当卖点；③L2 Sox10/2/9 + Fos 负向对照；④L3 三显著链接 OR/经验p/重叠基因；⑤L4 robust 4 药 + 置换 emp_p=0.33（不显著）+ vorinostat 背景 + 端到端闭环。
- [ ] 红线勿删：#1 图不优于线性；#3 边级不重现=卖点；#4 Level3=结构保守非表达激活；#5 Level4=假设生成非药效预测。
- [ ] 可选升级：获公开人卒中/CNS 损伤 scRNA 后，L3 可升到 GSVA/AUCell「表达激活保守」；或申请 HRA007397（仅外周 PBMC 炎症轴）作可选扩展。

*源产物：`rewiring_full.csv`、`cc_gse174_rewiring_full.csv`、`cc_gse225_rewiring_full.csv`、`human_module_enrich.json`、`drug_reversal_result.json`、`drug_perm_result.json`、`CROSS_COHORT_报告_2026-07-09.md`、`CROSS_SPECIES_报告_2026-07-09.md`、`Validation_章草稿_2026-07-09.md`、`StepD_文章级分析阐述_2026-07-09.md`。*
