# *Patterns* (Cell Press) — Submission Package Checklist

**Manuscript:** TSC-GNN (Recovering Temporal Regulatory Rewiring…) · v6 · 2026-07-11
**Package built:** 2026-07-11 · folder `Patterns_submission_package/`

This checklist maps every *Patterns* submission requirement to its status in the package. ✅ = satisfied · ⚠️ = satisfied but advise review · ❌ = blocking, must fix before submit.

---

## A. Required manuscript components (all inside `Manuscript_TSC-GNN_v6.docx`)

| # | Component | *Patterns* requirement | Status | Location / note |
|---|---|---|---|---|
| 1 | **Title** | descriptive, dataset-agnostic framing | ✅ | line 1 |
| 2 | **Author list + affiliations** | full, with equal-contribution note | ✅ | lines 3–9 |
| 3 | **Summary** | ≤150 words, single paragraph, no citations | ✅ (77 words) | line 30 |
| 4 | **Keywords** | ≤10 | ✅ (10) | line 34 |
| 5 | **Highlights** | ≤4 bullets, ≤85 chars each, no abbreviations | ✅ (4 bullets) | lines 17–22 |
| 6 | **eTOC blurb** | ≤50 words, 3rd person,大众可读 | ✅ (48 words) | line 26 |
| 7 | **Bigger picture** | 1–2 paragraphs, ≤300 words, non-expert | ✅ | lines 36–40 |
| 8 | **Terminology Ledger** | not required (optional box) | ✅ (kept as box) | lines 42–56 |
| 9 | **Main text (IMRaD)** | intro→results→discussion→methods | ✅ | full body |
| 10 | **Resource availability** | 3 subsections: Lead contact / Materials / Data+code | ✅ | lines 396–410 |
| 11 | **Author contributions** | CRediT taxonomy | ✅ | lines 414–418 |
| 12 | **Funding** | statement | ✅ | lines 422–424 |
| 13 | **Generative AI declaration** | mandatory for Cell Press | ✅ | line 432 |
| 14 | **Conflict of Interest** | statement | ✅ | lines 436–438 |
| 15 | **Ethical statement** | IRB/IACUC note | ✅ | lines 442–444 |
| 16 | **References** | sequential numbering, Cell Press style | ✅ (fixed this round) | lines 448–562; **see B below** |
| 17 | **Figure legends** | all figures | ✅ | lines 564–576 |
| 18 | **In-text citations** | superscript sequential numbers | ✅ (56 unique, aligned) | throughout |

---

## B. Open items before you hit "Submit"

### ❌ Blocking
1. **23 references still lack a DOI** (`[DOI: pending]`). Breakdown:
   - **16 are citable journal articles** (Crossref search could not disambiguate them automatically — see `Pending_DOIs_report.md`). Recommend resolving via **PubMed E-utilities** or a Zotero/EndNote batch, then a final direct-DOI verification pass.
   - **5 are non-Crossref** (ICLR/NeurIPS conference papers #15/#16/#17/#18/#38; bioRxiv preprint #40). For these, supply an **arXiv / OpenReview / PMLR / bioRxiv ID** instead of a Crossref DOI, or leave as "in press / preprint" with the identifier.
   - ⚠️ **Verify #54 Weider**: manuscript cites *Nat. Commun.* **12**:4240 (2021) but the highly-cited version is *Nat. Commun.* **9** (2018), DOI 10.1038/s41467-018-03336-3. Confirm which you mean.
2. **Editor name + suggested reviewers** in `Cover_Letter.md` are placeholders — fill before submitting.

### ⚠️ Advisory (non-blocking)
3. **Main-text length ≈ 9,800 words.** *Patterns* has no hard cap ("length aligns with the science"), and Methods articles with validation run long, but if a reviewer flags it, move Tables S1–S8 (already supplementary) and trim §4 discussion. No action required pre-submit.
4. **Graphical abstract is optional** for *Patterns* (it is on Cell Press's "optional GA" list). `Graphical_Abstract_TSC-GNN_Fig1.png` (the conceptual framework, Fig. 1) is included as a ready candidate; upload it only if you want a TOC image. Recommended — it improves discoverability.
5. **Declaration of Interests EM form**: the manuscript text declares no competing interests; the Editorial Manager system still requires the separate interest form to be completed at submission.
6. **GitHub URL placeholder**: the manuscript notes the repo goes public at acceptance; the Zenodo DOI (10.5281/zenodo.21289784) is already active and sufficient for the Data/Code availability check.

---

## C. Package contents

| File | What it is |
|---|---|
| `Manuscript_TSC-GNN_v6.docx` | Main text (title → references + figure legends), Patterns-formatted |
| `Supplementary_TSC-GNN_v6.md` | Supplementary Materials (Tables S1–S8, extended methods) |
| `Graphical_Abstract_TSC-GNN_Fig1.png` | Optional TOC/graphical abstract (conceptual framework) |
| `Cover_Letter.md` | Submission letter to the editor (fill placeholders) |
| `Pending_DOIs_report.md` | Tracker for the 23 references still missing identifiers |
| `Patterns_Submission_Checklist.md` | This file |

---

## D. Submission sequence (Editorial Manager)

1. Fill editor name + reviewers in `Cover_Letter.md`.
2. Resolve the 16 journal DOIs (item B1) — or consciously decide to submit with preprints/placeholders.
3. Upload `Manuscript_TSC-GNN_v6.docx` as the **main manuscript**; figures as separate files (Fig. 1 here doubles as GA if chosen).
4. Upload `Supplementary_TSC-GNN_v6.md` as **Supplementary Material**.
5. Paste/attach `Cover_Letter.md`.
6. Complete the **Declaration of Interests** form in the system.
7. Confirm Generative AI declaration is carried in the manuscript (it is, line 432).
