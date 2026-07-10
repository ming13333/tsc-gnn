#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
L5c positive control · In-context method validity in Replogle 2022 K562 CRISPRi
-------------------------------------------------------------------------------
Purpose: the L5c main result is a NULL for SOX10/CEBPB/GATA2 (their DoRothEA target
programs do NOT go down under their own K562 CRISPRi). To interpret that null as a
genuine CELL-TYPE CONTEXT boundary (and not a broken pipeline), we need a positive
control: TFs that ARE master regulators of K562 (an erythroleukemia line) should show
their DoRothEA regulon shifting DOWN under their own CRISPRi.

Design (identical machinery to l5c_replogle.py, no full-matrix load):
  For a curated panel of K562/erythroid-hematopoietic master TFs (positive controls)
  plus the 3 focal stroke TFs, compute:
    - self_z         : baseline-corrected effect of the TF's own locus (on-target KD check)
    - mean_target_z  : mean baseline-corrected Z of its DoRothEA target program
    - MWU one-sided  : targets shifted DOWN vs non-targets
    - Fisher OR      : targets over-represented in bottom-quartile (down) genes
    - rank / N       : rank of its own program among ALL candidate programs (1 = most down)
Interpretation:
  If in-context master TFs (GATA1/TAL1/KLF1/...) rank near the top (most-down) while the
  off-context stroke TFs do not, the method HAS power in K562 and the L5c null is a real
  context boundary, not a pipeline failure.
"""
import anndata, numpy as np, json, os, re
from scipy.stats import mannwhitneyu, fisher_exact

HERE = os.path.dirname(os.path.abspath(__file__))
H5   = os.path.join(HERE, "K562_gwps_normalized_bulk_01.h5ad")
DORO = r"C:\D 盘\科研\虚拟敲除\gate1\data\dorothea\human_dorothea_regulon.tsv"
OUT  = os.path.join(HERE, "l5c_positive_control_result.json")

# K562 = CML-derived erythroleukemia. Curated master/relevant TFs (positive controls):
POS_CTRL = ["GATA1", "TAL1", "KLF1", "MYB", "MYC", "RUNX1", "LMO2", "LDB1",
            "NFE2", "ZFPM1", "BCL11A", "ZBTB7A", "FLI1", "SPI1", "CEBPA", "E2F1", "MYC"]
FOCAL    = ["SOX10", "CEBPB", "GATA2"]
PANEL    = list(dict.fromkeys(POS_CTRL + FOCAL))  # dedup, keep order

# ---- load DoRothEA (human) ----
tgt = {}
with open(DORO, encoding="utf-8") as f:
    next(f)
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) < 2:
            continue
        tgt.setdefault(p[0], set()).add(p[1])

print("loading (backed)...", flush=True)
adata = anndata.read_h5ad(H5, backed="r")
print("shape:", adata.shape, flush=True)

obs_genes_raw = list(adata.obs.index)
def gene_of(name):
    parts = str(name).split("_")
    return parts[1].strip().upper() if len(parts) >= 2 else str(name).upper()
base_genes = np.array([gene_of(n) for n in obs_genes_raw])

ctrl_kw = ("NON-TARGETING", "INTERGENIC", "SAFE", "EMPTY", "MOCK",
           "UNTREATED", "UNKNOWN", "NEGATIVE", "NTC", "NONTARGETING")
ntc_mask = np.array([g == "NON-TARGETING" or g.startswith("NON-TARGETING")
                     or g in ctrl_kw or bool(re.match(r"^(INTERGENIC|CONTROL|SAFE|NEG|EMPTY|MOCK|UNTREATED|UNKNOWN|NTC)", g))
                     for g in base_genes])
ntc_idx = np.where(ntc_mask)[0]
print(f"NTC rows: {len(ntc_idx)}", flush=True)

# var Ensembl -> symbol
var_syms = list(adata.var.index)
if var_syms and str(var_syms[0]).startswith("ENSG"):
    import mygene
    mg = mygene.MyGeneInfo()
    q = mg.querymany(var_syms, scopes="ensembl.gene", fields="symbol", species="human", as_dataframe=True)
    symdict = q["symbol"].to_dict() if "symbol" in getattr(q, "columns", []) else {}
    sym_map = {s: (symdict[s] if (s in symdict and isinstance(symdict[s], str)) else None) for s in var_syms}
    eff_syms = [sym_map.get(s) for s in var_syms]
else:
    eff_syms = var_syms
sym2idx = {s: i for i, s in enumerate(eff_syms) if s}
print(f"var symbols mapped: {len(sym2idx)}/{len(var_syms)}", flush=True)

X = adata.X
ntc_rows = X[ntc_idx].toarray() if hasattr(X[ntc_idx], "toarray") else np.asarray(X[ntc_idx])
ntc_rows = np.nan_to_num(ntc_rows, nan=0.0, posinf=0.0, neginf=0.0)
ntc_mean = ntc_rows.mean(axis=0)

# candidate program universe (same as main analysis): perturbed genes with >=5 mapped targets
present_base = set(base_genes)
cand = [g for g in present_base if g in tgt and len([s for s in tgt[g] if s in sym2idx]) >= 5]
print(f"candidate programs: {len(cand)}", flush=True)

def analyze(tf):
    m = (base_genes == tf)
    n = int(m.sum())
    if n == 0:
        return {"present": False, "reason": "not perturbed in K562"}
    if tf not in tgt or len([s for s in tgt[tf] if s in sym2idx]) < 5:
        return {"present": True, "n_guides": n, "reason": "no DoRothEA regulon (>=5 mapped targets)"}
    rows = X[np.where(m)[0]]
    rows = rows.toarray() if hasattr(rows, "toarray") else np.asarray(rows)
    rows = np.nan_to_num(rows, nan=0.0, posinf=0.0, neginf=0.0)
    eff = rows.mean(axis=0) - ntc_mean
    self_i = sym2idx.get(tf)
    self_z = float(eff[self_i]) if self_i is not None else None

    tg = [sym2idx[s] for s in tgt[tf] if s in sym2idx]
    tgset = set(tg)
    nontg = [i for i in range(len(eff)) if i not in tgset]
    tg_z = eff[tg]; nontg_z = eff[nontg]
    U, p_mwu = mannwhitneyu(tg_z, nontg_z, alternative="less")
    thr = np.percentile(eff, 25)
    a = int((tg_z < thr).sum()); b = len(tg_z) - a
    c = int((nontg_z < thr).sum()); d = len(nontg_z) - c
    OR, p_fisher = fisher_exact([[a, b], [c, d]], alternative="greater")

    prog_means = {}
    for c_ in cand:
        ct = [sym2idx[s] for s in tgt[c_] if s in sym2idx]
        if len(ct) >= 5:
            prog_means[c_] = float(eff[ct].mean())
    all_means = np.array(list(prog_means.values()))
    self_mean = prog_means[tf]
    rank = int((all_means < self_mean).sum()) + 1
    N = len(prog_means)
    return {"present": True, "n_guides": n, "self_z": self_z,
            "n_targets_mapped": len(tg), "n_targets_total": len(tgt[tf]),
            "mean_target_z": float(tg_z.mean()),
            "mwu_p_less": float(p_mwu), "fisher_OR_down": float(OR), "fisher_p": float(p_fisher),
            "rank_among_programs": rank, "n_programs": N,
            "pct_most_down": float((N - rank) / N * 100.0)}

res = {"meta": {"dataset": "Replogle2022 K562_gwps (CRISPRi) positive control",
                "positive_controls": POS_CTRL, "focal": FOCAL,
                "n_candidate_programs": len(cand)},
       "panel": {}}
for tf in PANEL:
    r = analyze(tf)
    res["panel"][tf] = r
    if r.get("present") and "rank_among_programs" in r:
        print(f"[{tf:7s}] self_z={r['self_z'] if r['self_z'] is not None else float('nan'):.3f} "
              f"meanTgZ={r['mean_target_z']:+.4f} MWU_p={r['mwu_p_less']:.2e} "
              f"OR={r['fisher_OR_down']:.2f} rank={r['rank_among_programs']}/{r['n_programs']} "
              f"({r['pct_most_down']:.1f}% most-down)", flush=True)
    else:
        print(f"[{tf:7s}] {r.get('reason','?')}", flush=True)

with open(OUT, "w", encoding="utf-8") as f:
    json.dump(res, f, indent=2, ensure_ascii=False)
print("\nSAVED", OUT, flush=True)
