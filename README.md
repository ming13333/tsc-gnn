# TSC-GNN — Temporal / cell-State-Conditioned Graph Neural Network for Interpretable Regulatory Rewiring

**TSC-GNN** is a fixed-graph, linear-readout framework for recovering
*edge-level* transcriptional rewiring (ΔW with permutation significance) from
multi-timepoint, multi-condition scRNA-seq. It is designed as a **recovery
engine**, not a prediction engine: a fixed literature-curated causal GRN
(DoRothEA) is combined with a linear decoder so that the *which edge
changed and why* question stays interpretable. See the companion
manuscript (`gate1/rewiring_study/manuscript_v4.md`) for methods and
the five-level evidence ladder (L1–L5).

> **Naming note.** "GNN" here denotes the *structured causal substrate*
(the DoRothEA graph) on which a virtual perturbation is propagated and
decomposed — **not** a message-passing network with learnable
embeddings. The graph is fixed and the readout is linear by design
(see §5.10 of the manuscript). This is the source of the framework's
interpretability advantage and is distinct from end-to-end perturbation
GNNs such as GEARS.

## Install

```bash
conda env create -f environment.yml
conda activate tsc-gnn
```

The pinned environment reproduces the analyses under fixed seeds
(numpy 2.2.6 / scipy 1.15.3 / pandas 2.3.3).

## Pipeline (maps to manuscript sections)

| Step | Script(s) | Manuscript |
|------|----------|------------|
| Pre-processing | `gate1/gate1/preprocessing.py` | §5.1 |
| GRN construction (DoRothEA) | `gate1/gate1/data_acquisition.py`, `build_dorothea_grn` | §5.2 |
| Edge-level rewiring ΔW + permutation | `gate1/rewiring_study/build_cache.py`, `build_cache_single.py`, `cross_cohort_consistency.py` | §5.3 |
| Cell–cell communication | `gate1/rewiring_study/cc_diag.py` | §5.4 |
| Cross-species (L3) | `gate1/rewiring_study/build_disease_signature.py` | §5.5 |
| Drug reversal via LINCS (L4) | L1000CDS2 API client | §5.6 |
| Public TF-perturbation re-analysis (L5) | `gate1/audit_positive_control.py`, `gate1/*l5*` | §5.7 |
| Prediction benchmark (graph vs linear) | `gate1/gate1/baselines.py`, `evaluate.py` | §5.8 |

Each run emits a SHA-256 manifest (command, library versions, seeds,
cache hash) — see §5.9.

## Data

All inputs are public; no new sequencing was generated.

- Mouse stroke scRNA-seq: GSE174574, GSE225948
- Human stroke blood bulk: GSE16561
- TF knockout: GSE269122 (Sox10), GSE273163 (Cebpb)
- K562 genome-scale CRISPRi: Figshare 20029387

DoRothEA regulons are read from local TSV exports (see `gate1/data/`
or the `build_dorothea_grn` routine); the pipeline is fully offline.

## Repository layout

```
tsc-gnn-repo/
├── README.md            # this file
├── environment.yml       # pinned conda env
├── LICENSE              # MIT
├── .gitignore
└── gate1/              # analysis code (kept alongside the manuscript repo)
    ├── gate1/           # core modules: preprocessing, task_builder, baselines, evaluate, synthetic, data_acquisition
    ├── rewiring_study/  # ΔW / permutation / cross-cohort / cell-chat scripts + manuscript
    ├── gears_validation/ # GEARS comparison
    └── *.py            # L5 perturbation re-analysis, audit, positive control
```

> **Upload note.** Copy the `gate1/` tree (or symlink it) into this
> repository before `git push`. The GitHub and Zenodo URLs in the
> manuscript are placeholders to be finalised at acceptance.

## License

MIT — see `LICENSE`.
