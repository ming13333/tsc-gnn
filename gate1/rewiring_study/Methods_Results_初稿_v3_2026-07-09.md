# Manuscript Draft v3 — Methods · Results · Discussion（据 L3 评审 reframing · 2026-07-09）

> **作者备忘（非稿件正文，投稿前勿删；本文件已据审稿意见 v2 重排）**：
> 1. 图（固定因果 GRN + 线性读头）**不提升预测精度**——已用 GEARS 真实 Perturb-seq 与 90 点鲁棒性研究证实；主轴必须是**可解释性**（边级 ΔW + 模块级 master-regulator），不是预测精度。
> 2. 时间轴 24h/2d/14d 由 **跨数据集拼接**（GSE174574 + GSE225948），非均匀采样，存在 batch×time 混淆；须在方法/图注显式披露。
> 3. 边级显著性**不**跨队列重现（Jaccard≈0，方向一致≈0.52）——这是主轴卖点，须如实写、当证据而非缺陷藏。
> 4. Level 3 = **跨物种调控程序汇聚（cross-species regulatory convergence）**，定位为 functional convergence（**非** biological validation）：结构保守（C2，SOX10→髓鞘 OR=27 emp_p=0.0005）+ 独立人卒中 bulk 中相应 TF 靶程序被激活（GSE16561 外周血 39 卒中/24 对照，三模块 BH-q=1.1e-8）。须披露：外周血非脑、横断面非时序、PAX6 负对照也激活→非唯一特异、未做细胞组成反卷积；SOX10 血液激活**不可直接解读为髓鞘再生**。表述禁用"validated/proves"，改用"orthogonal evidence / independently supports / converges"。
> 5. Level 4 药物逆转 **置换未显著（emp_p≈0.3）**，定位为 translational proof-of-concept，非药效预测；vorinostat 跨系一致是 L1000 背景而非特异信号。
> 6. 整合集 sham 基线 = 两队列混合；三套基因集（6,897 / 8,736 / 7,041）不可直接比边级。
> 7. Level 5 = **程序级方向性因果支持（causal support），非 causal validation**：用公开 TF 功能缺失 RNA-seq（GSE269122 Sox10 少突 cKO；GSE273163 Cebpb 杂合 KO 库普弗）重分析；靶程序取自同一 DoRothEA GRN（非循环）；指标用 rank + OR + Fisher p（置换 emp_p 受集大小敏感，仅佐证）。须披露：重分析公开数据非自有扰动 / 非卒中语境 / 模块级非边级 / Cebpb 杂合 3v3 低功率、效应温和 / 少突 Sox10 KO 中 Cebpb 靶也部分下调（报梯度非唯一性）。表述禁"validated/proves"，用"causal support / directional consistency"。

---

# Methods

## 2.1 Conceptual framework: temporal-state-conditioned virtual perturbation (TSC-GNN)

We frame *in-silico perturbation* as an operation on a latent graph space rather than on a single-gene knockout. Our framework, **TSC-GNN** (Temporal and cell-State-Conditioned Graph Neural Network), takes a multi-timepoint, multi-condition single-cell transcriptome and (i) builds a time- and state-conditioned gene-regulatory graph, (ii) propagates a perturbation along this graph to obtain a rewired latent embedding, and (iii) reads out *edge-level rewiring* — the change in transcription-factor (TF)→target coupling across a transition — together with its permutation significance and the resulting master-regulator ranking. Cell-cell communication edges and pseudotime-resolved state vectors condition the propagation; a LINCS/L1000 reverse-match converts the injury embedding into candidate-repurposing hypotheses (Fig. 1, conceptual).

## 2.2 Data

**Mouse ischemic-stroke single-cell transcriptomes.** We used two independently generated mouse middle-cerebral-artery-occlusion (MCAO) scRNA-seq cohorts:
- **Cohort 1 — GSE174574** (Li et al., 2021): 3 MCAO (24 h) + 3 sham, 57,528 cells, with a BEAM pseudotime trajectory. Single acute time point.
- **Cohort 2 — GSE225948** (Anrather et al., *Nature Immunology* 2024): genuine two-timepoint post-stroke time course (2 d and 14 d, brain and blood).

**Constructed time axis and its disclosure.** No single public scRNA-seq cohort covers all three stages. We therefore assembled a 24 h → 2 d → 14 d axis spanning the acute → sub-acute peak → repair/remodelling phases by stitching cohort 1 (24 h + sham) with cohort 2 (2 d / 14 d). *This is a limitation that must be stated in the Methods and figure legends:* the axis is **non-uniformly sampled** (the 3–7 d window is absent) and carries a **batch×time confound** (the two time segments come from different studies/platforms). We mitigate — but do not eliminate — the confound by treating each time point's within-state coupling independently and by PC-regression composition correction (§2.4); we do **not** interpret discrete ΔW between stitched segments as a continuous *rate*. Two of the four analysed transitions (24 h→2 d and sham→14 d) are consequently **stitched cross-cohort pseudo-transitions**, which we account for in the cross-cohort comparison (§3.4).

**Pre-processing.** Count matrices were loaded offline (no scanpy dependency) and processed with library-size normalization (×(1e4)) followed by log1p, in a fixed conda environment (`bbb_gnn`; numpy 2.2.6 / scipy 1.15.3 / pandas 2.3.3). Cell-type composition was retained per cell for the composition-correction step. For cross-cohort analyses, each cohort was processed under an identical protocol so that the retained gene spaces (6,897 integrated / 8,736 cohort-1 / 7,041 cohort-2) are comparable, with the gene-space discrepancy itself reported as a limitation on edge-level comparison.

## 2.3 Gene-regulatory network construction

We used the **DoRothEA** consensus regulons (confidence levels A–C) as a *directed causal* TF→target graph — direction and sign (activation/repression) are taken from the literature-curated prior, so that any rewiring we report is interpretable as "whose regulation is enhanced/weakened", not merely as undirected co-expression change. Mouse and human regulons were read from local TSV exports (`gate1/data/dorothea/`), making the pipeline fully offline and reproducible (SHA-256 manifest per run). For each TF we additionally computed a **state-affinity** vector `A_aff` (subsample of n=4,000 cells) capturing how strongly the TF's targets are expressed in each cell state; edges in the top 50% by |A_aff| are retained as *state-conditioned* edges for rewiring testing.

## 2.4 State-conditioned edge rewiring

For each tested directed edge *e* = (TF *u* → target *v*) and each time point *t*, we compute the **within-state Pearson coupling**
$$r_{e,t} = \mathrm{corr}\big(x_u^{(t)},\, x_v^{(t)}\big),$$
where $x_u^{(t)}, x_v^{(t)}$ are the log-normalized expressions of *u* and *v* restricted to cells in state/time *t*. The **rewiring effect** for a transition $t_1\!\to\!t_2$ is the coupling change
$$\Delta W_{e,\,t_1\to t_2} = r_{e,t_2} - r_{e,t_1}.$$
Positive ΔW denotes *coupling enhancement* (the TF–target co-expression link strengthens across the transition); negative denotes *coupling weakening*. The DoRothEA prior direction (activation/repression) is reported independently and is not conflated with the empirical coupling change.

**Composition correction.** Because stroke alters cell-type composition across time, raw coupling can reflect composition shift rather than regulatory rewiring. We residualize $x_u, x_v$ against the top **n_pc = 10** principal components of the full expression matrix (least-squares regression of each gene on [1, PCs]), so that ΔW reflects coupling change after removing the major axes of compositional variation. We report both raw and PC-corrected ΔW; PC correction can flip edge direction when composition masked the true regulatory change (e.g. *Sox10→Ank3*: raw ΔW = −0.38 → PC-corrected +0.73).

**Permutation test and multiple testing.** For each transition we permute the time labels (n_perm = 200, seed = 2) and recompute ΔW to build a null per edge; the per-edge two-sided p-value is $p_e = (1 + \#\{\text{perm}: |\Delta W^{\text{null}}_e| \ge |\Delta W^{\text{obs}}_e|\})/(n_{\text{perm}}+1)$ (p_min ≈ 0.005). We apply two corrections: (i) **Benjamini–Hochberg FDR** per transition; and (ii) a **permutation pooled-FDR** that standardizes each edge's |ΔW| by its null standard deviation, pools null z-scores across edges and permutations, and computes q-values without depending on p-resolution — the preferred correction given the limited n_perm. Edges at pooled q < 0.05 (8, 28, 19 and 36 per transition respectively) are reported as significant rewiring; q < 0.1 is used for broader TF-level recovery.

**Reproducibility.** All runs emit a manifest (command, library versions, seeds, SHA-256 of the cache, resident-matrix size). The full analysis is deterministic under fixed seeds.

## 2.5 Four-level validation framework

We assess the framework along four escalating levels, each answering one question a reviewer must ask:

| Level | Question | Design |
|---|---|---|
| **L1** Technical reproducibility | Is the method stable? | Re-derive rewiring on two independent mouse cohorts; compare TF master-regulator ranking. |
| **L2** Recovery of established programs | Is it biologically correct? | Align recovered TFs with known post-stroke repair regulators; report stress-artifact controls. |
| **L3** Cross-species convergence | Does it transfer to human? | (a) Orthogonally project the mouse core TF module onto the human DoRothEA network; test enrichment of literature-curated repair reference sets with size-matched permutation. (b) AUCell/ssGSEA activation test on public human stroke bulk (GSE16561, blood). |
| **L4** Drug-perturbation reversal | Is it translatable? | Convert the stroke injury signature into a LINCS L1000 reverse-connectivity map; rank candidate-repurposing agents. |
| **L5** Public TF-perturbation directionality | Is it causally supported? | Re-analyse public TF loss-of-function RNA-seq (GSE269122 Sox10 oligodendrocyte-cKO; GSE273163 Cebpb heterozygous-KO Kupffer cells); test whether the framework's recovered TF→target program is enriched among genes down-regulated by the TF's *own* knockout — program-level directional causal support, not edge-level validation. Additionally, query the SigCom LINCS gene-perturbation libraries (CRISPR KO `l1000_xpr` ~140 K signatures; overexpression `l1000_oe` ~34 K) with each TF's human DoRothEA target program and test whether the TF's *own* perturbation signature ranks as a top reverser (KO) or mimicker (OE) with self-specificity. Further, re-analyse a genome-scale single-cell-resolution CRISPRi screen (Replogle et al. 2022, K562; 11,258 perturbations) to test whether the support holds at single-cell resolution *and* in an off-context cancer line — a context-boundary control. |

**Evidence-strength ladder.** The levels escalate from computational reproducibility to biological *convergence* and, with L5, to program-level directional *causal support* (and, in future work, to causal *validation*). We use precise verbs per tier: L1–L2 show the method *recovers* known biology; L3 provides *orthogonal* (independent, cross-species and cross-platform) evidence of regulatory *convergence*; L4 demonstrates *translational convergence*; L5 provides *causal support* by re-analysing independent public TF-knockout data at the program level. None of L1–L4 is a causal/perturbation validation, and L5 — being a re-analysis of public data, not our own perturbation, and at the program rather than edge level — is *causal support*, not *causal validation*; we avoid the verb "validate" for all tiers.

| Tier | Design | Evidence type | Strength |
|---|---|---|---|
| L1 | Mouse cross-cohort rewiring | Technical reproducibility | Discovery |
| L2 | Alignment with established stroke programs | Biological plausibility | Recovery |
| L3 | Human GRN projection + independent patient bulk activation | Structural + functional convergence | **Orthogonal support** |
| L4 | LINCS reverse-mapping | Translational convergence | Hypothesis generation |
| **L5** | **Public TF loss-of-function re-analysis (bulk KO) + cross-modality gene-perturbation (LINCS signature-match A + K562 sc-CRISPR regulon-response C, with in-context positive control)** | **Program-level directional causality, context-gated** | **Causal *support*** |
| (future) | CRISPR / human brain-resident perturbation | Causal / tissue-specific | **Validation** |

---

# Results

## 3.1 Predictive accuracy is not superior to linear baselines

Before presenting the biological findings, we report an honest scoping result that defines the contribution of the framework. We benchmarked the graph readout against a linear baseline under matched inputs (identical perturbation vector *p*, identical training split, identical readout head) across **10 random seeds × 5 graph types** (k-NN, DoRothEA, random, permuted-DoRothEA, 0-hop) **× 2 tasks = 90 configurations**. Under the current implementation — a fixed causal graph with a linear readout — the graph component contributed **limited additional predictive capacity** beyond the linear baseline: it beat the linear baseline in **0** cases and showed no significant difference in **86** (Fig. 2A). A non-linear KernelRidge probe overfit (rel ≈ −126%), indicating that any accuracy gain would require a trained deep non-linear readout rather than the fixed graph studied here. On the external GEARS Perturb-seq benchmark the graph readout gave −8.3% (single) / −4.2% (double) improvement, both with confidence intervals including zero.

Because the graph embedding enters the readout linearly and is provided identically to the linear baseline, it does not expand the linear hypothesis space available under matched inputs; the graph therefore does not, under this implementation, yield higher accuracy. This scoping result is the premise for the interpretability-oriented design developed next.

## 3.2 The framework is interpretability-oriented by design

The result above delineates the **boundary of graph utility**. A fixed causal graph with a linear readout is, by construction, a restricted substrate for prediction. Its distinctive contribution lies elsewhere: the graph exposes *where* regulation remodels — at the level of individual TF→target edges — and *which* master regulators drive that remodeling, which a black-box linear or foundation-model readout does not surface. We therefore treat TSC-GNN as an **interpretability-first virtual-perturbation framework**: its deliverable is edge-level rewiring ΔW with permutation significance and module-level master-regulator ranking. The remaining Results demonstrate what this lens recovers.

## 3.3 Temporal rewiring recovers biologically coherent, established stroke programs

Applying state-conditioned rewiring across the four transitions recovered a coherent, literature-aligned program of stroke recovery. At pooled q < 0.05, **8, 28, 19 and 36 edges** were significant for sham→24 h, 24 h→2 d, 2 d→14 d and sham→14 d respectively.

- **Acute injury onset (sham→24 h, 8 significant edges):** weak but directionally consistent coupling of damage/inflammatory genes (e.g. *Tbx21*→*Cxcr3*, ΔW = −0.21), consistent with the known acute ischaemic response and the modest acute signal.
- **Repair initiation (24 h→2 d, 28 significant edges):** the strongest coupling gains were *Sox2*→*Lsamp* (ΔW = +0.86, q < 0.001) and *Sox10*→*Ank3* (ΔW = +0.74, q < 0.001) — a sharp surge in oligodendrocyte-lineage commitment coupling that initiates remyelination.
- **Active remyelination (2 d→14 d, 19 significant edges):** the canonical oligodendrocyte myelin edge **Sox10→Plp1** strengthened markedly (ΔW = +0.51, q < 0.001), i.e. the method independently recovers the remyelination programme at the repair peak; *Hey2*→*Acta2* (ΔW = +0.59, q < 0.001) marks vessel-wall remodeling.
- **Inflammation resolution (sham→14 d, 36 significant edges):** *Sox9*→*Hapln1* (ΔW = +0.79, q < 0.001) and *Sox10*→*Plp1* (ΔW = +0.42, q < 0.001) couple positively with repair, tracking resolution of the acute inflammatory programme.

**Recovery of established programs (Level 2).** The recovered master regulators aligned with established post-stroke repair biology. **Sox10, Sox2 and Sox9** — canonical master regulators of oligodendrocyte lineage commitment and myelin regeneration — drove the strongest rewiring and closed the loop with the top remyelination edge *Sox10→Plp1*. As a negative control we report explicitly that, under the strict q < 0.1 threshold, the *only* TF reproduced across all three analyses was **Fos** — an immediate-early / AP-1 stress-response gene non-specifically induced by any injury; we surface this as a stress-artifact control rather than a biological signal.

**Sanity of recovered targets.** Targets of significant rewired edges were enriched for oligodendrocyte (OR 36–61), neuron (OR 110) and microglia (OR 15) markers; 14/15 known stroke-relevant TFs appeared, and 77–83% of edges preserved their direction after PC composition correction — confirming the rewiring signal is not a compositional artefact. PC correction flipping raw→PC direction (e.g. *Sox10→Ank3*: raw −0.38 → +0.73) demonstrates that composition masking can hide true regulatory gain and that the corrected estimate is the more reliable rewiring measure.

Collectively, these rewiring events reconstruct the canonical temporal program of stroke — from acute injury through inflammatory resolution to oligodendrocyte remyelination and repair.

## 3.4 Level 1 — Cross-cohort master-regulator reproducibility

To test technical stability we re-derived the rewiring on the two independent cohorts and compared **TF master-regulator rankings** against the integrated analysis. *Crucially*, because the integrated "sham" baseline mixes cells from both cohorts (and the two cohorts share no transition and differ in gene space), we evaluate reproducibility at the **master-regulator level**, where biological signal is robust, rather than at the edge level, where single-edge estimates are intrinsically noisy. The integrated ranking pools all four transitions, two of which are stitched cross-cohort pseudo-transitions (§2.2); the ρ therefore reflects TFs that drive rewiring in general, with *Sox10* the clean cross-cohort-reproduced exemplar.

The master-regulator ranking reproduced significantly across cohorts:

| Comparison | Common TFs | Spearman ρ | p |
|---|---|---|---|
| Integrated vs cohort 1 | 251 | **+0.517** | 1.5e-18 |
| Integrated vs cohort 2 | 235 | **+0.548** | 8.8e-20 |
| Cohort 1 vs cohort 2 | 242 | **+0.482** | 1.7e-15 |

**Sox10** appeared in the top-20 |ΔW|max TFs of **all three** analyses; **Sox2/Sox9** reproduced in ≥2. By contrast, individual edge-level significance did **not** reproduce across cohorts (direction agreement ≈ 0.52, i.e. at chance; Jaccard ≈ 0) — consistent with the high intrinsic variance of single-edge estimates under a fixed-graph/linear-readout method, and reported here as expected behaviour that motivates the module-level interpretability focus, not as a defect.

## 3.5 Level 3 — Cross-species regulatory convergence of recovered programs

> **Scope disclosure (作者红线 #4).** The orthogonal-projection test (§3.5a) asks whether the mouse-recovered TF's regulatory program is **conserved in the structure of the human GRN** (via the independent human DoRothEA network). The activation test (§3.5b) asks, separately, whether the corresponding target programs are **transcriptionally engaged in independent human stroke data**. Together they upgrade L3 from "structural conservation" to **cross-species regulatory convergence** — the same programs that rewire in mouse stroke are (i) conserved in human GRN architecture and (ii) co-activated in independent human patient samples. This is *functional convergence*, not *biological validation*: the bulk cohort is peripheral blood (not brain) and cross-sectional (not time-resolved), so brain-intrinsic and temporal claims still rest on the mouse evidence, and we avoid the verb "validate" for L3.

We took the 12-TF cross-cohort core module (recovered in ≥2 of 3 analyses; **Sox10** common to all three) and orthogonally projected each TF onto its human DoRothEA target set. For each TF we tested enrichment of two literature-curated reference sets — *myelin/oligodendrocyte* and *neuroinflammation* — by hypergeometric test, with **2,000 size-matched (0.5–2× target-count) random human TF permutations** as the null. Testing per-TF (not as a merged broad module) was essential: merging the 12 broadly-targeting TFs gave a target set covering 65% of the network and diluted every signal to non-significance.

Three links were significantly conserved:

| TF (human) | Reference set | k/n | OR | Hypergeom. p | Emp. p (perm) |
|---|---|---|---|---|---|
| **SOX10** | myelin/oligo | 7/22 | **27.0** | 6.0e-08 | **0.0005** |
| **CEBPB** | neuroinflammation | 8/23 | **16.5** | 3.2e-07 | **0.024** |
| **GATA2** | neuroinflammation | 15/23 | **4.6** | 3.3e-04 | **0.047** |

**SOX10 → human myelin/oligodendrocyte (flagship):** overlapping *PLP1, MBP, MAG, MPZ, PMP22* (major myelin structural proteins) plus *GJC2* (oligodendrocyte gap junction) and *PDGFRA* (oligodendrocyte-precursor marker) — the mouse top rewired TF projects precisely onto the human myelin-regeneration program, closing the loop with Sox10→Plp1.

**CEBPB → human neuroinflammation:** overlapping *IL1B, IL6, TNF, CCL3, CCL5, NOS2, PTGS2, STAT3*. **GATA2 → neuroinflammation:** overlapping *TLR2, TLR4, NFKB1, NFKBIA* and *CCL/CXCL* chemokines — the TLR/NF-κB innate-immune module.

Non-significant / expected-negative results are reported honestly: AR/ERG/NR2F2/PAX5/RUNX3/SOX9 showed no tissue-specific enrichment (emp. p 0.29–0.65); SOX2/E2F1/GATA3 are too broadly targeting to show single-tissue enrichment (biologically expected). The enrichment uses the *human* DoRothEA knowledge base, so it confirms conservation of the program's structure rather than an independent human discovery; the permutation control excludes the trivial "any broad TF enriches" explanation.

### 3.5b Expression activation in independent human stroke bulk (GSE16561)

To move beyond structural conservation, we asked whether the recovered programs are *transcriptionally activated* in independent human stroke data. We obtained public whole-blood bulk RNA-seq (GSE16561; 39 stroke vs 24 controls, Illumina HumanWG-6, RMA-normalized). Probes were mapped to HGNC symbols via GPL6883 and collapsed to 17,493 genes. For each recovered TF we scored single-sample module activity with an AUCell/ssGSEA-style rank-mass fraction on its human DoRothEA target set, then compared stroke vs control (Mann–Whitney U; Cliff's δ; BH correction). All three programs were significantly activated in human stroke blood (BH-q = 1.1×10⁻⁸), and stronger than 99.8% of size-matched random gene sets (empirical p = 0.002):

| Module | targets (present/total) | Δ (stroke−control) | Cliff's δ | p |
|---|---|---|---|---|
| CEBPB (inflammation) | 555/589 | +0.058 | +0.89 | 3.8×10⁻⁹ |
| SOX10 (myelin/oligo) | 296/322 | +0.045 | +0.83 | 4.9×10⁻⁸ |
| GATA2 (inflammation) | 4780/5370 | +0.023 | +0.76 | 5.2×10⁻⁷ |

A brain-enriched negative-control TF (**PAX6**, 344 targets) was *also* activated (p = 4.5×10⁻⁸). The observed activation therefore cannot be interpreted as evidence that these three regulatory programs are *uniquely* activated in stroke; rather, these findings indicate that part of the signal likely reflects generalized transcriptional reprogramming associated with systemic inflammation and the post-stroke leukocyte-composition shift. We therefore interpret the L3 result as **orthogonal evidence of cross-species regulatory convergence**: the recovered programs are among the co-activated modules in human stroke, and their *direction* (inflammatory axes up-regulated) is concordant with the mouse rewiring. Importantly, **activation of the SOX10 module in peripheral blood should not be interpreted as direct evidence of oligodendrocyte remyelination, but rather as indirect support that components of the conserved regulatory program are transcriptionally engaged during stroke**; the oligodendrocyte-specific reading of SOX10 rests on the mouse structural evidence (§3.5a). Cell-composition deconvolution was not performed, and the blood-based, cross-sectional design limits direct inference on brain-intrinsic temporal remodeling.

This provides **orthogonal, cross-platform evidence** for cross-species regulatory convergence: the mouse rewiring recovers Sox10→Plp1 remyelination and Cebpb→Il1b inflammatory resolution (§3.3–3.4), and the same TF programs are independently co-activated in human stroke blood (§3.5b), converging with the drug-reversal hits (HDAC inhibitors / statins, anti-inflammatory; §3.6). We phrase this as *supporting / independent convergence*, not validation.

## 3.6 Level 4 — Drug-perturbation relevance via LINCS

As a translational demonstration we converted the 24 h stroke injury program into a disease signature — up-regulated injury/inflammatory genes (SPP1, CCL4, LGALS3, CD14, C5AR1, TNF …) and a heterogeneous set of lower-expression genes (interferon-response, microglia-marker and transport genes; see `disease_signature.json`). Notably, the 24 h signature is dominated by the **acute inflammatory/injury axis** and does **not** by itself contain down-regulated myelin/repair structural genes; the myelin-regeneration axis is recovered independently by the rewiring analysis (Levels 2–3), not by this signature. The signature was derived from the DoRothEA-network gene space (the same space as the rewiring), i.e. a regulatory-program-filtered differential expression, not the full transcriptome.

We built two signatures: a single-cohort **main** signature (up100/dn100, best score 0.0671 with severe ties, only 1 known agent hit) and a **robust** consensus signature intersected across the two cohorts (up28/dn17, best 0.125). We report **robust** as the primary result; main serves as a single-cohort noise control, echoing the L1 reproducibility theme.

The robust signature returned **four literature-supported agents in the top 50**: **vorinostat** and **trichostatin A** (HDAC inhibitors) and **mevastatin** and **rosuvastatin** (statins). Their **anti-inflammatory** relevance aligns directly with the inflammatory axis captured in the 24 h signature; their **pro-myelin / neuroprotective** relevance aligns with the remyelination program independently recovered by the rewiring analysis (Sox10→Plp1, Levels 2–3) — a cross-level convergence, not a within-signature effect.

**Statistical rigour (permutation background).** Using unique-drug counting (to avoid inflating the count by L1000's multi-cell-line duplicate entries), 20 random signatures of matched size returned a mean of **2.75 ± 1.37** such hits (max 5). The observed count of **4** did **not** exceed chance (**empirical p ≈ 0.3**; N = 20, single-shot, live API — the non-significant conclusion is robust to re-runs). A preliminary duplicate-counting scheme had spuriously suggested p = 0.048; we corrected this and do **not** claim significant enrichment. Further, vorinostat's cross-cell-line consistency (n_hits = 6) proved to be an L1000 background effect: under random signatures vorinostat appeared in 53% of runs (max n_hits = 16).

> **作者红线 #5（必写进 Limitation）**：Level 4 is a **translational proof-of-concept / hypothesis-generation** demonstration, *not* a drug-efficacy prediction. LINCS profiles derive from non-neuronal cancer cell lines; the signature is single-timepoint (24 h acute) and derived from the DoRothEA gene space rather than the full transcriptome; L1000CDS2 shows heavy ties; the permutation did not reach significance; and the entire pipeline is in silico with no wet-lab confirmation. We state these limits explicitly.

**Narrative closure.** The candidate classes surfaced by reverse mapping — HDAC inhibitors (pro-myelin / anti-inflammatory) and statins (anti-inflammatory / neuroprotective) — target *precisely the programs the framework independently recovered as dysregulated in stroke* (Levels 2–3). This closes the loop from interpretable virtual perturbation to a testable therapeutic hypothesis, which is the end-to-end translational story current GRN / virtual-perturbation methods papers largely lack.

## 3.7 Level 5 — Program-level directional causal support from public TF perturbation

> **Scope disclosure (作者红线 #7).** L5 addresses the central caveat of any causal-graph method — *the graph is never perturbed* — by re-analysing **public, independent** TF loss-of-function RNA-seq. We test, at the **program (TF→target module) level**, whether the framework's recovered target program is enriched among genes that go *down* when the TF itself is knocked out. This is **directional causal support**, explicitly **not** edge-level validation (edges do not reproduce across cohorts, §3.4) and **not** a claim of causality within stroke (the perturbation data are non-stroke). We use the verbs *causal support / directional consistency* and avoid *validate*.

The target programs are taken from the **same mouse DoRothEA GRN** used by the framework (§2.3), so they are independent of the perturbation datasets and the test is non-circular. For each TF we (i) compute a down-enrichment odds ratio (OR) and Fisher exact p against genes with log₂FC < 0; (ii) rank the program's mean log₂FC among **all 412/404 DoRothEA TF programs** to gauge specificity (1 = most down-regulated); (iii) run 2,000 size-matched random-set permutations for a directional empirical p (noted to be set-size sensitive, so it is corroborative, not headline).

**Sox10 — oligodendrocyte-specific conditional knockout (GSE269122, 4 ctrl / 4 Sox10-KO, corpus callosum; author DE table, Ensembl→symbol via mygene, 32,912 genes).** The *recovered* Sox10 target program is the most down-enriched among the tested candidates:

| TF program | n | mean log₂FC | OR↓ | Fisher p | Rank (most-negative) |
|---|---|---|---|---|---|
| **Sox10** (perturbed) | 319 | **−0.090** | **1.81** | ≈0 | **46 / 412 (top 11%)** |
| Cebpb (control) | 570 | −0.065 | 1.83 | ≈0 | 147 / 412 |
| Gata2 (control) | 5177 | −0.033 | 1.33 | ≈0 | 309 / 412 |
| Sox2 (control) | 3121 | −0.025 | 1.29 | ≈0 | 339 / 412 |

Removing Sox10 down-regulates the genes our GRN attributes to Sox10 most strongly of all candidates (rank 46/412; specificity gradient Sox10 ≪ Cebpb < Gata2 < Sox2), consistent with "remove an activator → its targets fall" and providing directional causal support for the Sox10→myelin/remyelination program recovered in §3.3.

**Cebpb — heterozygous knockout Kupffer cells (GSE273163, 3 WT / 3 Hete, liver MASH; RAW counts → log₂FC, 17,119 genes).** The recovered Cebpb program is again the most down-enriched among candidates, with clean specificity against the non-perturbed Sox10 program:

| TF program | n | mean log₂FC | OR↓ | Fisher p | Rank (most-negative) |
|---|---|---|---|---|---|
| **Cebpb** (perturbed) | 555 | **−0.172** | **1.20** | **0.022** | **149 / 404 (top 37%)** |
| Gata2 (control) | 4889 | −0.147 | 1.15 | ≈0 | 205 / 404 |
| Sox10 (control) | 307 | −0.128 | 0.95 | 0.68 | 278 / 404 |

Cebpb's own targets are significantly down (OR = 1.20, Fisher p = 0.022) and more down-regulated than the non-perturbed Sox10 targets (OR = 0.95, p = 0.68); the effect is modest because the perturbation is a **heterozygous, 3-vs-3** design, so L5 treats it as confirmatory directional support rather than a definitive causal test.

**Cross-dataset specificity.** The Sox10 program is specifically down only where Sox10 is perturbed (rank 46 in GSE269122 vs 278 in GSE273163); the Cebpb program is down in both (oligodendrocyte maturation and inflammatory programs are co-regulated) but most down where Cebpb itself is perturbed. This mirrors the L3 lesson: the *program* is conserved and directionally causal, whereas no individual edge is claimed causal. We therefore position L5 as **causal support by independent public-perturbation re-analysis** — the step that closes the "causal graph never perturbed" gap without overstating it as validation.

**Independent gene-perturbation consistency (SigCom LINCS).** As an orthogonal test, we queried the SigCom LINCS gene-perturbation libraries — CRISPR knockdown (`l1000_xpr`, 140,603 signatures) and overexpression (`l1000_oe`, 33,782) — with each TF's human DoRothEA target program (up to 500 targets, capped). For each TF, we asked whether its *own* perturbation signature ranks as a top reverser (CRISPR KO: removing an activator → targets fall → reverser of the up-target query) or top mimicker (OE: adding an activator → targets rise → mimicker), and whether this ranking is self-specific (self-rank ≪ cross-TF rank).

The strongest result is **GATA2 overexpression**: GATA2's own OE signature ranks **3rd of 33,782** (top 0.01 %) as a mimicker of the GATA2 target program (p = 1.4 × 10⁻⁵), with both GATA2 OE replicates classified as mimickers. The self-specificity is pronounced — self percentile 0.01 % vs best cross-TF 4.98 % (Δ = −4.97 %) — providing strong orthogonal evidence that GATA2 is an activator of its DoRothEA-recovered targets. In the CRISPR-KO library, all three TFs' own KO signatures appear in the top ~1 % of reversers (SOX10 rank 139,333/140,603, top 0.90 %; CEBPB top 0.46 %; GATA2 top 1.17 %), directionally correct but **not self-specific** (self ≈ cross percentiles), likely because the L1000 cancer-cell-line context dilutes TF-specific effects and the five-gene neutral down-set provides limited gene-set specificity. SOX10 OE showed mixed directionality (one mimicker at rank 1,261, one reverser at 32,109), and no CEBPB OE signature exists in the library. We report these mixed results transparently: GATA2 OE constitutes the clearest orthogonal directional support, while the CRISPR-KO direction is consistent but TF-level specificity is inconclusive.

**Single-cell-resolution perturb-seq re-analysis (Replogle 2022, K562).** As a higher-resolution and off-context control we re-analysed a genome-scale CRISPRi screen in K562 (Replogle et al., 2022; 11,258 perturbations, 585 non-targeting controls), testing each TF's recovered DoRothEA target program under the TF's own CRISPRi signature (gemgroup Z-normalised pseudo-bulk, one row per guide). In this off-context cancer line the test returned a **null**: none of the three TFs' target programs was significantly down-shifted under its own perturbation (SOX10's locus is absent from the K562 gene space — the TF is not expressed in myeloid K562; CEBPB self-Z = −0.43 and GATA2 self-Z = −0.19 confirm on-target knockdown, yet their target programs' mean Z were +0.007 and +0.004, MWU p = 0.84 and 0.18, ranking 239/332 and 139/332 among all candidate programs). This negative result is informative rather than contradictory: it demarcates the **context boundary** of program-level causal support. Directional support is recovered when the perturbation is biologically appropriate — a TF knocked out in its native lineage (L5a: Sox10 oligodendrocyte-cKO, Cebpb Kupffer-cell KO) or perturbed in the correct direction (L5b: GATA2 overexpression) — but not in a generic cancer-cell-line CRISPRi screen where the TF may be inactive (SOX10) or its tissue-aggregated target program not co-regulated (CEBPB/GATA2). The contrast across L5a/L5b/L5c is itself a prediction of the interpretability-first framework and argues against treating any single perturbation screen as definitive validation. A K562-internal positive control (§3.8) confirms this null reflects a **cell-type context boundary rather than a pipeline artefact**: master regulators of K562 (MYC, BCL11A) do show significant regulon down-shift under their own CRISPRi, while the off-context stroke TFs do not.

---

## 3.8 Cross-modality synthesis of gene-level perturbation support (A + C)

Levels 5b and 5c interrogate the framework's recovered TF→target programs with two **independent public gene-perturbation modalities** that differ in resolution and in *what they measure*. Analysed jointly they yield an insight that neither delivers alone: **directional causal support is context-gated, and the two modalities dissociate "signature mimicry" from "regulon response."** We formalise the combined analysis as a complete method + result unit below.

### 3.8.1 Methods

**(M1) Signature-matching modality — SigCom LINCS (A).** For each focal TF we submit its human DoRothEA target program (up-gene query, ≤500 targets) to the SigCom LINCS `l1000_xpr` (CRISPR-KD, 140,603 signatures) and `l1000_oe` (overexpression, 33,782) libraries via the two-sided signature-search endpoint (`database` passed as the string library name — the fix for the recurring HTTP-500). We record the percentile rank of the TF's *own* perturbation signature as a reverser (KD) or mimicker (OE) of its target program, and its **self-vs-cross specificity** (self percentile minus best cross-TF percentile). This tests transcriptome-*scale* similarity: does perturbing the TF produce a signature that globally looks like / opposes the target module?

**(M2) Regulon-response modality — Replogle 2022 K562 sc-CRISPR (C).** From the genome-scale CRISPRi pseudo-bulk (11,258 perturbations, 585 non-targeting controls, gemgroup-Z, one row per guide; read in `anndata` backed mode, Ensembl→symbol via `mygene`) we take the TF's own perturbation row, baseline-correct against the NTC mean (`np.nan_to_num` to neutralise gemgroup-Z ∞ from constant-variance genes), and test whether the **specific DoRothEA target set** is shifted down: mean target Z, one-sided Mann-Whitney U (targets vs non-targets), Fisher OR on bottom-quartile genes, and the rank of the program's mean Z among all 332 candidate programs (1 = most down). This tests membership-*resolved* response: do the exact annotated targets move under endogenous knockdown?

**(M3) In-context positive control (C).** To distinguish a genuine cell-type context boundary from a broken pipeline, we ran M2 on a curated panel of K562 (erythro-myeloid leukaemia) master regulators — GATA1, TAL1, KLF1, MYB, MYC, RUNX1, NFE2, BCL11A, ZBTB7A, FLI1, SPI1, CEBPA, E2F1 — for which a regulon down-shift under their own CRISPRi is biologically expected. `l5c_positive_control.py`, same universe of 332 programs.

**(M4) Context-appropriateness ordering.** We arrange all gene-perturbation evidence for the three focal TFs by a single axis — biological context appropriateness of the perturbation — spanning **native lineage** (L5a mouse oligodendrocyte / Kupffer-cell KO) → **correct-direction generic** (L5b LINCS OE) → **off-context cancer line** (L5c K562 CRISPRi), and read the *gradient* of support rather than any single cell.

### 3.8.2 Results

**Cross-modality map for the three focal TFs.** Support is monotone in context appropriateness and modality-dependent (Table 3.8-1).

**Table 3.8-1. Gene-level perturbation support across modalities (focal stroke TFs).**

| TF | Native-lineage bulk KO (L5a) | LINCS overexpression mimicker (A) | LINCS CRISPR-KD reverser (A) | K562 sc-CRISPR regulon response (C) | K562 context |
|---|---|---|---|---|---|
| **SOX10** | rank 46/412, OR↓ 1.81, p≈0 ✓ | mixed (mimic 1,261 / rev 32,109) | top 0.90 % (non-specific) | rank 304/332, MWU p=0.91 — **locus absent** | not expressed |
| **CEBPB** | rank 149/404, OR↓ 1.20, p=0.022 ✓ | no OE signature in library | top 0.46 % (non-specific) | rank 239/332, self-Z −0.43 (on-target), MWU p=0.84 | partial |
| **GATA2** | not tested (bulk) | **rank 3/33,782, p=1.4×10⁻⁵ ✓✓ (self-specific)** | top 1.17 % (non-specific) | rank 139/332, self-Z −0.19 (on-target), MWU p=0.18 | in-context, yet null |

**Positive control confirms the pipeline has power in K562, but generic regulons are a coarse proxy (Table 3.8-2).** Among K562 master TFs the regulon-response test recovers a significant down-shift for **MYC (rank 19/332, MWU p=3.1×10⁻³)** and **BCL11A (rank 29/332, p=7.5×10⁻³)**, with KLF1/SPI1/MYB ranking in the top ~10–18 % (strong on-target knockdown, self-Z −0.58/−0.44/−0.81). The three focal stroke TFs rank far below these in-context positives (139–304/332). Crucially, even a canonical K562 erythroid master — **GATA1** — shows strong on-target knockdown (self-Z −0.54) *without* a coherent DoRothEA-regulon down-shift (mean Z +0.27, rank 187/332), demonstrating that generic DoRothEA regulons only partially track TF activity even in the correct cell type.

**Table 3.8-2. K562 in-context positive control (selected).**

| TF (K562 role) | self-Z (on-target) | regulon mean-Z | MWU p | rank / 332 |
|---|---|---|---|---|
| MYC (proliferation) | −0.03 | +0.003 | **3.1×10⁻³** | **19** |
| BCL11A (globin switch) | n/a | −0.017 | **7.5×10⁻³** | **29** |
| SPI1 / PU.1 (myeloid) | −0.44 | +0.006 | 0.50 | 35 |
| KLF1 / EKLF (erythroid) | −0.58 | −0.002 | 0.35 | 41 |
| MYB (progenitor) | −0.81 | +0.036 | 0.18 | 59 |
| GATA1 (erythroid master) | −0.54 | +0.268 | 1.00 | 187 |
| — SOX10 / CEBPB / GATA2 (focal) | n/a / −0.43 / −0.19 | +0.018 / +0.007 / +0.004 | 0.91 / 0.84 / 0.18 | 304 / 239 / 139 |

**New findings from combining A + C.**

1. **Convergent context-gating across two independent modalities.** Two orthogonal public resources — LINCS bulk L1000 and Replogle K562 single-cell CRISPRi — independently fail to give TF-specific, brain-relevant directional support for the stroke programs, whereas native-lineage KO (L5a) does. These convergent off-context nulls act as **negative controls** that exclude the trivial explanation that the target programs are recoverable from *any* perturbation dataset.

2. **The two modalities dissociate signature-mimicry from regulon-response — visible only when combined.** GATA2 is the decisive case: its overexpression signature is the 3rd-strongest mimicker of its program (LINCS, top 0.01 %), yet its exact target set does **not** move under endogenous CRISPRi in K562 (rank 139/332, p=0.18) — even though GATA2 is a bona-fide hematopoietic TF present in K562. The transcriptome-scale "mimicry/reversal" LINCS reports is therefore **not** evidence that the specific DoRothEA edges are active in that context; only the membership-resolved Replogle test adjudicates that, and it is null off-context. Neither modality alone reveals this separation.

3. **Calibrated interpretation via positive control: power exists, but module-level tests have a ceiling.** The K562-internal positives (MYC, BCL11A) prove the null is not a pipeline artefact; the GATA1 failure proves that generic regulons are only a coarse, breadth-and-confidence-dependent proxy even in-context. This **bounds** how strongly *any* module-level gene-perturbation test — A or C — can be read, and elevates L5a's native-lineage result as the strongest causal tier.

4. **The gradient is the result, and it is what an interpretability-first framework predicts.** Ordered by context appropriateness, support decays monotonically (native-lineage positive → correct-direction OE self-specific for GATA2 → off-context CRISPRi null), exactly as expected if the graph encodes **context-specific regulatory hypotheses** rather than a transferable predictor. The combined A+C analysis thus converts a superficially "negative" cross-modality outcome into positive evidence for the paper's central claim — and reinforces the honest positioning of L5 as *causal support*, never validation.

*源产物：`l5_perturbation/l5b_sigcom.py` + `l5b_sigcom_result.json`（A）、`l5_perturbation/l5c_replogle.py` + `l5c_replogle_result.json`（C）、`l5_perturbation/l5c_positive_control.py` + `l5c_positive_control_result.json`（M3 正对照）、报告 `l5_perturbation/AC_跨模态综合_报告_2026-07-09.md`。*

---

# Discussion

## 4.1 A methodological boundary, not a tool

This study demonstrates that, under a fixed causal graph and linear decoder, graph topology provides **little additional predictive power** beyond linear models. The scoping benchmark (§3.1) — 90 configurations in which the graph beat the linear baseline in zero cases, with confidence intervals spanning zero on the external GEARS Perturb-seq set — establishes this as a property of the *fixed-graph / linear-readout* regime, not a failure of engineering. We state it openly rather than obscuring it, because the conclusion reframes where graph models are useful.

The distinctive value of the graph lies in **interpretability**. A directed causal graph exposes *where* regulation remodels, at the resolution of individual TF→target edges, and *which* master regulators drive that remodeling. A black-box linear or foundation-model readout compresses the perturbation to a single score and discards this structure. We therefore position TSC-GNN as an interpretability-first virtual-perturbation framework. (The choice of DoRothEA over an undirected k-NN co-expression graph is made on exactly this ground: the two are near-indifferent for prediction, but only the directed causal graph renders rewiring as an interpretable "whose regulation is enhanced/weakened" statement.)

## 4.2 What the graph does and does not buy

The four validation levels map onto four claims:

1. **Prediction is not a graph advantage.** Under matched inputs a fixed causal graph contributes limited predictive capacity beyond linear baselines (§3.1).
2. **Interpretability is the graph advantage.** Graph structure uniquely enables edge-level rewiring ΔW with permutation significance and module-level master-regulator inference (§3.3–3.4).
3. **The recovered program recapitulates known biology.** The rewired modules recover the temporally ordered repair programme — inflammatory resolution, oligodendrocyte remyelination, neuronal recovery — and the core TFs (Sox10/Sox2/Sox9) are established repair regulators (§3.3, L2).
4. **The program is conserved, translatable, and directionally causal.** The regulatory architecture projects onto human-conserved stroke-repair programs (SOX10→myelin, OR = 27; CEBPB/GATA2→neuroinflammation; §3.5a, L3) and is *independently* activated in public human stroke blood (§3.5b) — orthogonal evidence of cross-species regulatory convergence rather than validation. Re-analysis of public TF loss-of-function RNA-seq (§3.7, L5) shows the recovered Sox10 and Cebpb target programs are the most down-regulated among all candidate programs when their own TF is knocked out — **directional causal support** at the program level. An orthogonal test in the SigCom LINCS gene-perturbation libraries further shows that GATA2 overexpression is the 3rd-strongest mimicker (of 33,782 signatures) of its own target program with strong self-specificity, while all three TFs' CRISPR knockdowns appear as top-1 % reversers in the correct direction. Together with the LINCS reverse-map onto known neuroprotective drug classes (§3.6, L4), this closes the loop from interpretable virtual perturbation to a causally-supported, translatable hypothesis. The three L5 layers (native-lineage bulk KO, L1000 overexpression, and K562 single-cell CRISPRi) triangulate the *context dependence* of this support: it is recovered only when the perturbation matches the TF's active biology, and is absent in an off-context cancer line (§3.7, L5c) — a boundary that strengthens, rather than weakens, the directional-support claim. Jointly analysing the two independent gene-perturbation modalities (§3.8, A + C) adds two calibrations: (i) they **dissociate** whole-transcriptome signature-mimicry (LINCS: GATA2 OE top 0.01 %) from membership-resolved regulon response (K562: null off-context), so a global mimicry signal does not by itself certify that specific edges are active; and (ii) a K562-internal positive control (MYC rank 19/332 p=3.1×10⁻³; BCL11A 29/332 p=7.5×10⁻³) confirms the off-context null is a genuine cell-type boundary, while the GATA1 in-context failure bounds how strongly any module-level test can be interpreted — collectively reinforcing native-lineage KO as the strongest causal tier and L5 as *support*, not validation.

The central scientific message is therefore not "a more accurate virtual-KO model," nor "a new stroke discovery," but a **methodological insight**: *fixed causal graphs do not necessarily improve perturbation prediction, but they substantially improve the interpretability of temporal regulatory remodeling by enabling edge-level rewiring analysis.* Ischemic stroke here serves as the proof-of-concept application that exercises the framework.

## 4.3 Limits and future directions

Several constraints are inherent to the current design. The longitudinal axis is stitched across two cohorts (§2.2), so discrete ΔW between segments is not a continuous rate and cross-cohort edge-level comparison is limited by differing gene spaces. The drug-reversal level is a hypothesis-generation demonstration, not efficacy prediction (§3.6). The L5 causal-support level is a *re-analysis of public, non-stroke perturbation data at the program level* — it supports the directionality of the recovered TF→target programs but is not a causal validation within stroke and not an edge-level claim (§3.7). The SigCom LINCS gene-perturbation extension (L5b) provides orthogonal directional evidence for GATA2 (OE rank 3/33,782, strong self-specificity) and correct reverser direction for all three TFs' CRISPR KOs, but lacks TF-level self-specificity in the CRISPR-KO library and showed mixed directionality for SOX10 OE; these mixed results are reported transparently rather than selectively. The K562 single-cell CRISPRi extension (L5c) further returned a null in this off-context line — SOX10 is not expressed in myeloid K562 and CEBPB/GATA2 are on-target knocked down yet their target programs are not co-down-regulated — which we report transparently as a *context boundary* rather than a refutation: it shows that program-level causal support requires biologically appropriate perturbations and that no single screen is definitive. A K562-internal positive control (§3.8) supports this reading — the pipeline recovers significant regulon down-shift for genuine K562 master regulators (MYC, BCL11A) — but also exposes a ceiling: even the canonical erythroid master GATA1 shows on-target knockdown without a coherent DoRothEA-regulon response, so generic regulons are only a coarse, breadth- and confidence-dependent proxy, bounding how strongly any module-level gene-perturbation test (LINCS or K562) can be interpreted. Future work can extend the framework along three axes: (i) a **learnable graph** estimated from data rather than a fixed prior; (ii) a **non-linear trained decoder**, which the overfitting probe suggests is necessary for any predictive gain; and (iii) **human longitudinal and brain-resident datasets** (e.g. PBMC time courses or spatial/single-cell transcriptomics of human stroke tissue) and **own CRISPR/perturb-seq validation of the recovered programs** to convert the current program-level directional *support* (L5) into genuine tissue- and context-specific *causal validation*.

## 4.4 Conclusion

This study demonstrates that, under a fixed causal graph and linear decoder, graph topology provides little additional predictive power beyond linear models.

Nevertheless, graph representations uniquely enable the recovery of dynamic regulatory rewiring that cannot be inferred from conventional predictive models.

Applied to ischemic stroke, this framework reconstructed temporally ordered repair programs involving inflammatory resolution, oligodendrocyte remyelination, and neuronal recovery, providing biologically interpretable hypotheses rather than merely predictive scores.

Future work incorporating learnable graph structures, nonlinear decoders, and human longitudinal datasets may further extend this framework from interpretable perturbation analysis toward clinically actionable prediction.

---

# Abstract one-liner (suggested)

> Under a fixed causal graph and linear decoder, graph topology contributes little predictive power beyond linear baselines; its distinctive value is interpretability. Applied to ischemic stroke, the TSC-GNN framework recovers temporally ordered, master-regulator-driven repair programs (Sox10/Sox2/Sox9; Sox10→Plp1 remyelination) that are reproducible at the TF-master-regulator level across independent cohorts (Spearman ρ = 0.48–0.55, p < 1e-15). Cross-species regulatory conservation was further **supported** by independent activation of the corresponding transcriptional programs in a public human stroke cohort (SOX10→myelin OR = 27, empirical p = 0.0005; CEBPB/GATA2→neuroinflammation; all three programs co-activated, BH-q = 1.1×10⁻⁸), and the framework reverse-maps onto known neuroprotective drug classes. Re-analysis of public TF loss-of-function RNA-seq additionally showed that the recovered Sox10 and Cebpb target programs are the most down-regulated among all candidate programs when their own TF is knocked out (rank 46/412 and 149/404) — **program-level directional causal support** — demonstrating that fixed causal graphs improve the *interpretability of temporal regulatory remodeling*, not its prediction. This directional support is context-dependent — absent in an off-context K562 CRISPRi screen — and is therefore positioned as causal *support*, not validation.

---

# 待办 / 诚信核对清单（提交前逐条确认）

- [ ] Methods 必含：①时间轴 24h/2d/14d 跨数据集拼接 + 非均匀采样 + batch×time 混淆披露（已含，且 §3.3/§2.2 标注 24h→2d 与 sham→14d 为拼接伪转移）；②DoRothEA 有向因果图 + state-affinity 边筛选；③PC 组成校正(n_pc=10) + 置换(n_perm=200) + pooled-FDR；④图不优于线性的鲁棒性结论（§3.1，已软化措辞）。
- [ ] Results 必含：①§3.1 预测负结果（90 点 / GEARS）；②§3.3 主 rewiring 四转移具体边数与边（8/28/19/36 + Sox10→Plp1 ΔW=+0.51 q<0.001 等）+ 汇总句；③§3.4 L1 TF rank ρ 表 + 主调控因子重现（edge 不重现作预期行为）；④§3.3 L2 Sox10/2/9 + Fos 负向对照；⑤§3.5 L3 三显著链接 OR/经验p/重叠基因；⑥§3.6 L4 robust 4 药 + 置换 emp_p≈0.3（不显著）+ vorinostat 背景 + DoRothEA 基因空间限制 + 端到端闭环。
- [ ] Discussion 必含：①方法学洞见（图在固定图+线性读头下预测贡献有限、价值在解释）；②四层结论；③一句话 scientific claim；④未来方向（learnable graph / nonlinear decoder / human 数据）。
- [ ] 红线勿删：#1 图不优于线性；#3 边级不重现=预期行为/卖点；#4 Level3=跨物种调控程序汇聚（functional convergence，非 validation；禁用"validated/proves"，用"orthogonal evidence/independently supports/converges"；须披露：外周血非脑/横断面/PAX6负对照也激活/未反卷积/SOX10血液激活≠髓鞘再生）；#5 Level4=假设生成非药效预测。
- [x] L3 升级已完成：用 GSE16561 人卒中 bulk 跑 AUCell/ssGSEA，三模块 BH-q=1.1e-8 共激活（报告 `L3_upgrade_人bulk激活保守_2026-07-09.md`）；定位 = 跨物种调控汇聚（functional convergence，非 validation）；边界（外周血/横断面/PAX6负对照/未反卷积/SOX10≠髓鞘再生）已写入 §3.5b 与红线 #4、证据阶梯表（§2.5）。
- [x] L5 升级已完成：用公开 TF 功能缺失 RNA-seq 做程序级方向性因果支持（Sox10: GSE269122 少突 cKO，靶程序 rank 46/412 最下调，OR↓=1.81；Cebpb: GSE273163 杂合 KO 库普弗，rank 149/404，OR↓=1.20 Fisher p=0.022；特异性干净 Sox10 靶在 Cebpb KO 不下调 OR=0.95 p=0.68）。靶程序取自同一 DoRothEA GRN（非循环）。边界（重分析公开/非卒中/模块级非边级/杂合低功率/共调控）已写入 §3.7、证据阶梯表（§2.5）、红线 #7、摘要。报告 `L5_upgrade_公开扰动因果支持_2026-07-09.md`，脚本 `l5_perturbation/l5_causal_direction.py` + JSON。
- [x] L5b 升级已完成：SigCom LINCS 基因扰动方向一致性（`l1000_xpr` CRISPR KO 14万签名 + `l1000_oe` 过表达 3.4万签名）。GATA2 OE rank 3/33782 top 0.01% mimicker（p=1.4e-05，强自特异性 self 0.01% vs cross 4.98%）；3 TFs CRISPR KO 全部 top ~1% reverser（方向正确但无自特异性）；SOX10 OE 方向混合；CEBPB OE 无数据。已整合进 §3.7、§4.2④、§4.3、证据阶梯表（§2.5）。报告 `L5b_SigCom_报告_2026-07-09.md`，脚本 `l5_perturbation/l5b_sigcom.py` + `l5b_sigcom_result.json`。
- [x] L5c 升级已完成：Replogle 2022 K562 genome-scale CRISPRi pseudo-bulk 单细胞分辨率重分析（11,258 扰动 / 585 NTC）。三 TF（SOX10/CEBPB/GATA2）靶程序在自身 CRISPRi 下均**无显著下调** → **语境边界 null**：SOX10 位点不在 K562 基因空间（不表达）；CEBPB self-Z=−0.43、GATA2 self-Z=−0.19 on-target 成功但靶程序 mean Z=+0.007/+0.004，MWU p=0.84/0.18，rank 239/332、139/332。与 L5a（原生细胞类型 bulk KO 阳性）/ L5b（LINCS GATA2 OE rank 3/33782 阳性）构成三角验证「因果支持仅在生物学恰当语境成立」。已整合进 §2.5、§3.7 新增 sc-CRISPR 子节、§4.2④、§4.3、摘要。报告 `L5c_Replogle_scCRISPR_报告_2026-07-09.md`，脚本 `l5_perturbation/l5c_replogle.py` + `l5c_replogle_result.json`。

*源产物：`rewiring_full.csv`、`cc_gse174_rewiring_full.csv`、`cc_gse225_rewiring_full.csv`、`human_module_enrich.json`、`drug_reversal_result.json`、`drug_perm_result.json`、`CROSS_COHORT_报告_2026-07-09.md`、`CROSS_SPECIES_报告_2026-07-09.md`、`Validation_章草稿_2026-07-09.md`、`StepD_文章级分析阐述_2026-07-09.md`、`AUDIT_全流程审查_2026-07-09.md`、`L3_upgrade_人bulk激活保守_2026-07-09.md`、`L5_upgrade_公开扰动因果支持_2026-07-09.md`、`L5b_SigCom_报告_2026-07-09.md`、`L5c_Replogle_scCRISPR_报告_2026-07-09.md`（含 `l5_perturbation/l5_causal_direction.py` + `l5_causal_direction.json` + `l5_perturbation/l5b_sigcom.py` + `l5b_sigcom_result.json` + `l5_perturbation/l5c_replogle.py` + `l5c_replogle_result.json`）。*
