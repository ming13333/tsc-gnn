# Cover Letter — Submission to *Patterns* (Cell Press)

**Date:** [submission date]
**Manuscript type:** Methods article with validation
**Manuscript title:** Recovering Temporal Regulatory Rewiring: An Interpretable Graph-Based Virtual Perturbation Framework Applied to Ischemic Stroke

---

**To:** [Editor Name], Editor, *Patterns*
Cell Press, Elsevier

Dear Dr./Prof. [Editor Name],

Please find enclosed our manuscript entitled *"Recovering Temporal Regulatory Rewiring: An Interpretable Graph-Based Virtual Perturbation Framework Applied to Ischemic Stroke,"* which we submit for consideration as a **Methods article with validation** in *Patterns*.

## Why this work fits *Patterns*

*Patterns* publishes computationally rigorous, reproducible methods that turn data into usable, interpretable biology. TSC-GNN sits squarely in this remit. It is a method — not merely an application — that:

1. **Recovers edge-level gene-regulatory-network (GRN) rewiring** from multi-timepoint, multi-condition single-cell transcriptomes, with permutation significance;
2. **Couples that rewiring** to master-regulator inference, cell–cell communication remodelling (CellChat), and a LINCS L1000-based drug-repurposing map; and
3. **Validates the recovered programs** along an explicit five-level evidence ladder (L1–L5) spanning cross-cohort reproducibility, cross-species conservation, translational convergence, and independent public-perturbation causal support.

## The methodological point we most want reviewers to weigh

We report, with full transparency, that under a fixed causal graph and a linear readout, graph topology does **not** improve perturbation prediction beyond linear baselines (0 of 90 configurations). Rather than bury this negative result, we argue it *is* the contribution: it reframes what graph-based perturbation methods should be *for* — interpretability, not incremental accuracy — and we pair it with honest boundary-setting (edge-level output is presented as hypotheses to test, not settled causal truth). We believe this stance is of direct interest to the *Patterns* readership, which values methods that clarify rather than obscure mechanism.

## Novelty and scope

Ischemic stroke is used as the proof-of-concept that exercises the framework; the method itself is dataset-agnostic and applicable to any temporal single-cell remodelling question (development, neurodegeneration, cancer, immunology).

## Declarations

- The work is original, has not been published previously, and is not under consideration elsewhere.
- All data are publicly available, previously published omics datasets; no new human participants, specimens, or animal subjects were involved, so IRB/IACUC approval was not required.
- The authors declare no competing financial interests.
- In preparing the manuscript we used generative-AI assistance (WorkBuddy) solely for citation matching, reference formatting, and illustrative figures; all scientific content, data analysis, and interpretation are the authors'. The full statement is included in the manuscript.
- Code is released under an MIT licence and archived on Zenodo (DOI: 10.5281/zenodo.21289784); the GitHub repository will be made public upon acceptance.

## Suggested reviewers (optional — please confirm independence)

- [Reviewer 1 — expertise: single-cell gene-regulatory-network inference / SCENIC-family methods]
- [Reviewer 2 — expertise: graph neural networks for perturbation prediction]
- [Reviewer 3 — expertise: ischemic stroke molecular pathology / repair biology]
- [Reviewer 4 — expertise: causal inference from public perturbation atlases (Perturb-seq / LINCS)]

We would be happy to provide any additional materials needed for evaluation.

Sincerely,

**Wenwu Huang, M.D.**
Corresponding Author
Department of Neurosurgery, The Fifth Affiliated Hospital of Wenzhou Medical University (Lishui Central Hospital), Lishui, Zhejiang, China
Email: hwenwu321@gmail.com
