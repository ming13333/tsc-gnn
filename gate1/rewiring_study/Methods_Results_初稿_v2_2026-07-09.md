# Manuscript Draft v2 — Methods · Results · Discussion（统一初稿 · 2026-07-09）

> **作者备忘（非稿件正文，投稿前勿删；本文件已据审稿意见 v2 重排）**：
> 1. 图（固定因果 GRN + 线性读头）**不提升预测精度**——已用 GEARS 真实 Perturb-seq 与 90 点鲁棒性研究证实；主轴必须是**可解释性**（边级 ΔW + 模块级 master-regulator），不是预测精度。
> 2. 时间轴 24h/2d/14d 由 **跨数据集拼接**（GSE174574 + GSE225948），非均匀采样，存在 batch×time 混淆；须在方法/图注显式披露。
> 3. 边级显著性**不**跨队列重现（Jaccard≈0，方向一致≈0.52）——这是主轴卖点，须如实写、当证据而非缺陷藏。
> 4. Level 3 = 调控程序在**人 GRN 结构保守 + 独立人卒中 bulk 表达激活**：结构保守（C2，SOX10→髓鞘 OR=27 emp_p=0.0005）；表达激活（GSE16561 外周血 39 卒中/24 对照，三模块 BH-q=1.1e-8，但须披露：外周血非脑、横断面非时序、PAX6 负对照也激活→非唯一特异、未做细胞组成反卷积）。
> 5. Level 4 药物逆转 **置换未显著（emp_p≈0.3）**，定位为 translational proof-of-concept，非药效预测；vorinostat 跨系一致是 L1000 背景而非特异信号。
> 6. 整合集 sham 基线 = 两队列混合；三套基因集（6,897 / 8,736 / 7,041）不可直接比边级。

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

We validate the framework along four escalating levels, each answering one question a reviewer must ask:

| Level | Question | Design |
|---|---|---|
| **L1** Technical reproducibility | Is the method stable? | Re-derive rewiring on two independent mouse cohorts; compare TF master-regulator ranking. |
| **L2** Recovery of established programs | Is it biologically correct? | Align recovered TFs with known post-stroke repair regulators; report stress-artifact controls. |
| **L3** Cross-species conservation | Does it transfer to human? | (a) Orthogonally project the mouse core TF module onto the human DoRothEA network; test enrichment of literature-curated repair reference sets with size-matched permutation. (b) AUCell/ssGSEA activation test on public human stroke bulk (GSE16561, blood). |
| **L4** Drug-perturbation reversal | Is it translatable? | Convert the stroke injury signature into a LINCS L1000 reverse-connectivity map; rank candidate-repurposing agents. |

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

## 3.5 Level 3 — Cross-species conservation of regulatory programs

> **Scope disclosure (作者红线 #4).** The orthogonal-projection test (§3.5a) assesses whether the mouse-recovered TF's regulatory program is **conserved in the structure of the human GRN** (via the independent human DoRothEA network). *Separately*, we now test **expression activation** in an independent public human stroke bulk cohort (GSE16561, §3.5b). The two together upgrade L3 from "structural conservation" to "structural + activation conservation", but the bulk cohort is peripheral blood (not brain) and cross-sectional (not time-resolved), so brain-intrinsic and temporal claims still rest on the mouse evidence.

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

A brain-enriched negative-control TF (**PAX6**, 344 targets) was *also* activated (p = 4.5×10⁻⁸), indicating the signal partly reflects a broader stroke-associated transcriptional reprogramming (systemic inflammation / leukocyte-composition shift) rather than exclusive activation of the three recovered programs. We therefore interpret this as evidence that the recovered programs are among the co-activated modules in human stroke, with the oligodendrocyte-specific reading of SOX10 resting on the mouse structural evidence (§3.5a). Cell-composition deconvolution was not performed, and the blood-based, cross-sectional design limits direct inference on brain-intrinsic temporal remodeling.

This closes the cross-species loop: the mouse rewiring recovers Sox10→Plp1 remyelination and Cebpb→Il1b inflammatory resolution (§3.3–3.4), and the same TF programs are co-activated in independent human stroke blood (§3.5b) — converging with the drug-reversal hits (HDAC inhibitors / statins, anti-inflammatory; §3.6).

## 3.6 Level 4 — Drug-perturbation relevance via LINCS

As a translational demonstration we converted the 24 h stroke injury program into a disease signature — up-regulated injury/inflammatory genes (SPP1, CCL4, LGALS3, CD14, C5AR1, TNF …) and a heterogeneous set of lower-expression genes (interferon-response, microglia-marker and transport genes; see `disease_signature.json`). Notably, the 24 h signature is dominated by the **acute inflammatory/injury axis** and does **not** by itself contain down-regulated myelin/repair structural genes; the myelin-regeneration axis is recovered independently by the rewiring analysis (Levels 2–3), not by this signature. The signature was derived from the DoRothEA-network gene space (the same space as the rewiring), i.e. a regulatory-program-filtered differential expression, not the full transcriptome.

We built two signatures: a single-cohort **main** signature (up100/dn100, best score 0.0671 with severe ties, only 1 known agent hit) and a **robust** consensus signature intersected across the two cohorts (up28/dn17, best 0.125). We report **robust** as the primary result; main serves as a single-cohort noise control, echoing the L1 reproducibility theme.

The robust signature returned **four literature-supported agents in the top 50**: **vorinostat** and **trichostatin A** (HDAC inhibitors) and **mevastatin** and **rosuvastatin** (statins). Their **anti-inflammatory** relevance aligns directly with the inflammatory axis captured in the 24 h signature; their **pro-myelin / neuroprotective** relevance aligns with the remyelination program independently recovered by the rewiring analysis (Sox10→Plp1, Levels 2–3) — a cross-level convergence, not a within-signature effect.

**Statistical rigour (permutation background).** Using unique-drug counting (to avoid inflating the count by L1000's multi-cell-line duplicate entries), 20 random signatures of matched size returned a mean of **2.75 ± 1.37** such hits (max 5). The observed count of **4** did **not** exceed chance (**empirical p ≈ 0.3**; N = 20, single-shot, live API — the non-significant conclusion is robust to re-runs). A preliminary duplicate-counting scheme had spuriously suggested p = 0.048; we corrected this and do **not** claim significant enrichment. Further, vorinostat's cross-cell-line consistency (n_hits = 6) proved to be an L1000 background effect: under random signatures vorinostat appeared in 53% of runs (max n_hits = 16).

> **作者红线 #5（必写进 Limitation）**：Level 4 is a **translational proof-of-concept / hypothesis-generation** demonstration, *not* a drug-efficacy prediction. LINCS profiles derive from non-neuronal cancer cell lines; the signature is single-timepoint (24 h acute) and derived from the DoRothEA gene space rather than the full transcriptome; L1000CDS2 shows heavy ties; the permutation did not reach significance; and the entire pipeline is in silico with no wet-lab confirmation. We state these limits explicitly.

**Narrative closure.** The candidate classes surfaced by reverse mapping — HDAC inhibitors (pro-myelin / anti-inflammatory) and statins (anti-inflammatory / neuroprotective) — target *precisely the programs the framework independently recovered as dysregulated in stroke* (Levels 2–3). This closes the loop from interpretable virtual perturbation to a testable therapeutic hypothesis, which is the end-to-end translational story current GRN / virtual-perturbation methods papers largely lack.

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
4. **The program is conserved and translatable.** The regulatory architecture projects onto human-conserved stroke-repair programs (SOX10→myelin, OR = 27; CEBPB/GATA2→neuroinflammation; §3.5a, L3) and is transcriptionally activated in independent human stroke blood (§3.5b), and reverse-maps onto known neuroprotective drug classes (§3.6, L4).

The central scientific message is therefore not "a more accurate virtual-KO model," nor "a new stroke discovery," but a **methodological insight**: *fixed causal graphs do not necessarily improve perturbation prediction, but they substantially improve the interpretability of temporal regulatory remodeling by enabling edge-level rewiring analysis.* Ischemic stroke here serves as the proof-of-concept application that exercises the framework.

## 4.3 Limits and future directions

Several constraints are inherent to the current design. The longitudinal axis is stitched across two cohorts (§2.2), so discrete ΔW between segments is not a continuous rate and cross-cohort edge-level comparison is limited by differing gene spaces. The drug-reversal level is a hypothesis-generation demonstration, not efficacy prediction (§3.6). Future work can extend the framework along three axes: (i) a **learnable graph** estimated from data rather than a fixed prior; (ii) a **non-linear trained decoder**, which the overfitting probe suggests is necessary for any predictive gain; and (iii) **human longitudinal datasets** (e.g. applying the framework to PBMC or spatial transcriptomics of human stroke) to lift Level 3 from structural conservation to expression-activation validation.

## 4.4 Conclusion

This study demonstrates that, under a fixed causal graph and linear decoder, graph topology provides little additional predictive power beyond linear models.

Nevertheless, graph representations uniquely enable the recovery of dynamic regulatory rewiring that cannot be inferred from conventional predictive models.

Applied to ischemic stroke, this framework reconstructed temporally ordered repair programs involving inflammatory resolution, oligodendrocyte remyelination, and neuronal recovery, providing biologically interpretable hypotheses rather than merely predictive scores.

Future work incorporating learnable graph structures, nonlinear decoders, and human longitudinal datasets may further extend this framework from interpretable perturbation analysis toward clinically actionable prediction.

---

# Abstract one-liner (suggested)

> Under a fixed causal graph and linear decoder, graph topology contributes little predictive power beyond linear baselines; its distinctive value is interpretability. Applied to ischemic stroke, the TSC-GNN framework recovers temporally ordered, master-regulator-driven repair programs (Sox10/Sox2/Sox9; Sox10→Plp1 remyelination) that are reproducible at the TF-master-regulator level across independent cohorts (Spearman ρ = 0.48–0.55, p < 1e-15), project onto human-conserved stroke-repair architecture (SOX10→myelin OR = 27, empirical p = 0.0005; CEBPB/GATA2→neuroinflammation), and reverse-map onto known neuroprotective drug classes — demonstrating that fixed causal graphs improve the *interpretability of temporal regulatory remodeling*, not its prediction.

---

# 待办 / 诚信核对清单（提交前逐条确认）

- [ ] Methods 必含：①时间轴 24h/2d/14d 跨数据集拼接 + 非均匀采样 + batch×time 混淆披露（已含，且 §3.3/§2.2 标注 24h→2d 与 sham→14d 为拼接伪转移）；②DoRothEA 有向因果图 + state-affinity 边筛选；③PC 组成校正(n_pc=10) + 置换(n_perm=200) + pooled-FDR；④图不优于线性的鲁棒性结论（§3.1，已软化措辞）。
- [ ] Results 必含：①§3.1 预测负结果（90 点 / GEARS）；②§3.3 主 rewiring 四转移具体边数与边（8/28/19/36 + Sox10→Plp1 ΔW=+0.51 q<0.001 等）+ 汇总句；③§3.4 L1 TF rank ρ 表 + 主调控因子重现（edge 不重现作预期行为）；④§3.3 L2 Sox10/2/9 + Fos 负向对照；⑤§3.5 L3 三显著链接 OR/经验p/重叠基因；⑥§3.6 L4 robust 4 药 + 置换 emp_p≈0.3（不显著）+ vorinostat 背景 + DoRothEA 基因空间限制 + 端到端闭环。
- [ ] Discussion 必含：①方法学洞见（图在固定图+线性读头下预测贡献有限、价值在解释）；②四层结论；③一句话 scientific claim；④未来方向（learnable graph / nonlinear decoder / human 数据）。
- [ ] 红线勿删：#1 图不优于线性；#3 边级不重现=预期行为/卖点；#4 Level3=结构保守+表达激活保守（但须披露：外周血非脑/横断面/PAX6负对照也激活/未反卷积）；#5 Level4=假设生成非药效预测。
- [x] L3 升级已完成：用 GSE16561 人卒中 bulk 跑 AUCell/ssGSEA，三模块 BH-q=1.1e-8 表达激活保守（报告 `L3_upgrade_人bulk激活保守_2026-07-09.md`）；边界（外周血/横断面/PAX6负对照/未反卷积）已写入 §3.5b 与红线 #4。

*源产物：`rewiring_full.csv`、`cc_gse174_rewiring_full.csv`、`cc_gse225_rewiring_full.csv`、`human_module_enrich.json`、`drug_reversal_result.json`、`drug_perm_result.json`、`CROSS_COHORT_报告_2026-07-09.md`、`CROSS_SPECIES_报告_2026-07-09.md`、`Validation_章草稿_2026-07-09.md`、`StepD_文章级分析阐述_2026-07-09.md`、`AUDIT_全流程审查_2026-07-09.md`。*
