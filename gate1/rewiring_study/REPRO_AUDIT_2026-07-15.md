# Reproducibility Audit & Humanizer Review — 2026-07-15

Skill stack: `ai-research-repro-guard` (reviewer mode) + `academic-paper-humanizer` (academic de-AI).
Scope: GSE197731 cross-cohort CellChat analysis + `manuscript_v7.md` / `manuscript_v7.docx`.

## A. Repro-guard — independent verification (not self-audit)

| Item | Check | Result |
|---|---|---|
| Download integrity | `gzip -t` on all 8 GSE197731 samples | **8/8 valid**; cell counts 4769/5344/5245/5759/4189/5184/4372/5083 — exact, **no n-drift** |
| Reproducibility of subsampling | `cellchat_py._subsample_ct` | Uses **deterministic systematic sampling** (`step=len/cap; sel=idx[floor(arange(cap)*step)]`) — **no RNG**, results perfectly reproducible without a seed |
| GSE197731 WT arm (CSV) | recomputed | 24h 59 rec / **38 sig**; 48h 49 rec / **42 sig**; arm total 108 rec / **80 sig** ✓ |
| GSE197731 KO arm (CSV) | recomputed | 24h 60 / **45**; 48h 56 / **41** ✓ |
| Cross-cohort 24h | `cross_cohort_24h_consistency.txt` | LR Jaccard = **0.000**; pathway Jaccard = **0.067**; ctpair shared = **Pericyte→Pericyte**; KO 38/45, shared 30, WT-unique 8, KO-unique 15 ✓ |
| GSE225948-only (CSV) | recomputed excluding 24h row (GSE174574) | 83 rec / **61 sig**; 2d∩14d = **22**; linked_to_MR = **51/61 (84%)**; both∩linked = **18** |

## B. Errors found and fixed (silent errors caught)

1. **GSE225948 significant-count overstatement (root cause: source CSV mixed cohorts).**
   `cellchat_rewiring_sig.csv` holds 86 rows = GSE225948 (2d+14d) **+ GSE174574 24h (3 rows)**. The manuscript had attributed the combined 63-sig to GSE225948 alone.
   - Fixed `83 LR–transition records (63 significant)` → **(61 significant)**.
   - Fixed `53 of 63 (84 %)` → **`51 of 61 (84 %)`** (recomputed `linked_to_MR` from the CSV).
   - The total "194 records (143 significant)" was already correct (61+80+2).

2. **Pathway Jaccard rounding.** Manuscript stated `0.07`; data file gives `0.067`. Changed both occurrences (§3.9 + S2 caption) to **0.067** for exact reproducibility.

3. **docx citation bug (build script, pre-existing — NOT from humanizer).** `^NN` citations nested inside a `**bold**` lead-in (e.g. `**(GSE269122 ^24).**` in §3.7a) were swallowed by the bold token and rendered as literal `^24`/`^25` in the docx. Fixed `parse_inline` in `build_docx_v7.py` to **recursively tokenize inner content** of bold/italic spans. Verified: zero literal citation `^NN` remain; `^24`/`^25` now render as superscript.

## C. Humanizer review (academic de-AI, citations/structure preserved)

Manuscript is already specific and honest (real gene names, numbers, method detail; candid "null", "not significant", "stress-artifact" framing). AI-trace scan found **no** filler vocabulary (novel/promising/crucial/delve/underscore/…). Targeted edits only:
- Varied the repetitive connective triad `This positions/aligns X as` (×3) → `Taken together… acts as`, `We therefore treat… as`, `Read this way, TSC-GNN belongs to`.
- Softened a self-congratulatory clause ("benefits the field more than another incrementally more accurate predictor" → "serves the field better than another incrementally more accurate predictor would").

**Citations untouched**: all `^NN` / `<sup>` markers preserved; 56-entry reference list intact (Kim/Redox Biol 2022 present as ^29). No numbers, tables, equations, or section structure changed.

## D. Pre-submission reproducibility table (status)

| Item | Status |
|---|---|
| All citations | ✓ DOI/PMID (GSE197731 PMID 35688114 verified) |
| All data | ✓ raw GEO public; analysis CSVs in repo |
| Code | ✓ frozen scripts + git versioned (commit 58f029a) |
| Environment | ⚠ pinned in text (numpy 2.2.6 / scipy 1.15.3 / pandas 2.3.3); `environment.yml` referenced but not shipped in this pass |
| Random seed | ✓ not required (deterministic subsampling) |
| Statistical method | ✓ MWU + permutation + BH/pooled-FDR |
| Key results | ✓ independently recomputed (this audit) |
| AI-generated text | ✓ human-reviewed + de-AI pass applied |
