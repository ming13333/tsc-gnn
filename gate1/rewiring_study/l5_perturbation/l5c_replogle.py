#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
L5c · Public sc-CRISPR / perturb-seq re-analysis (Replogle 2022, K562 genome-scale CRISPRi)
-------------------------------------------------------------------------------------------
Single-cell-resolution directional causal support for the framework's recovered TF->target
programs. Uses Replogle et al. 2022 K562 genome-scale CRISPRi pseudo-bulk (one row per
guide/perturbation, gemgroup Z-normalized).

Data structure (confirmed):
  adata.obs.index (row names) = "<rowid>_<GENE>_<guide>_<ENSG>"  -> perturbation gene = field [1]
  adata.var.index = Ensembl IDs (ENSG...)                          -> map to symbol via mygene
  NTC rows: gene field == "non-targeting" (585 rows)

Logic (mirrors L5 bulk-KO + L5b SigCom, single-cell-derived):
  CRISPRi knockdown of an activator TF  => its direct targets go DOWN  => targets are
  enriched among the negatively-shifted genes under that TF's own perturbation.
  For each query TF q (SOX10 / CEBPB / GATA2):
    (a) self-consistency: is q's recovered DoRothEA target program the MOST down-shifted
        program among all candidate programs under q's own CRISPRi signature?
        -> rank of mean(target Z) among all programs + Mann-Whitney U + Fisher OR.
    (b) cross-specificity: is q's program more down-shifted than the other TFs' programs?
  Positioned as single-cell-resolution *causal support* (program level), NOT validation
  (K562 cancer line, not brain, not stroke; CRISPRi knockdown not knockout).

Reads only the needed rows via anndata backed mode (no full-matrix load).
"""
import anndata, numpy as np, json, os, re
from scipy.stats import mannwhitneyu, fisher_exact

HERE = os.path.dirname(os.path.abspath(__file__))
H5   = os.path.join(HERE, "K562_gwps_normalized_bulk_01.h5ad")
DORO = r"C:\D 盘\科研\虚拟敲除\gate1\data\dorothea\human_dorothea_regulon.tsv"
OUT  = os.path.join(HERE, "l5c_replogle_result.json")

TFS = ["SOX10", "CEBPB", "GATA2"]

# ---- load full DoRothEA (human) ----
tgt = {}
with open(DORO, encoding="utf-8") as f:
    next(f)
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) < 2:
            continue
        tgt.setdefault(p[0], set()).add(p[1])

# ---- open h5ad in backed mode (no full-matrix load) ----
print("loading (backed)...", flush=True)
adata = anndata.read_h5ad(H5, backed="r")
print("shape:", adata.shape, flush=True)

# perturbation gene = field[1] of obs index ("<rowid>_<GENE>_<guide>_<ENSG>")
obs_genes_raw = list(adata.obs.index)
def gene_of(name):
    parts = str(name).split("_")
    return parts[1].strip().upper() if len(parts) >= 2 else str(name).upper()
base_genes = np.array([gene_of(n) for n in obs_genes_raw])
print("n perturbations:", len(base_genes), flush=True)

# NTC detection: explicit non-targeting token + control keywords (start-anchored to
# avoid false hits on real genes like NONO / CNTN1 / KNTC1)
ctrl_kw = ("NON-TARGETING", "INTERGENIC", "SAFE", "EMPTY", "MOCK",
           "UNTREATED", "UNKNOWN", "NEGATIVE", "NTC", "NONTARGETING")
ntc_mask = np.array([g == "NON-TARGETING" or g.startswith("NON-TARGETING")
                     or g in ctrl_kw or re.match(r"^(INTERGENIC|CONTROL|SAFE|NEG|EMPTY|MOCK|UNTREATED|UNKNOWN|NTC)", g)
                     for g in base_genes])
ntc_idx = np.where(ntc_mask)[0]
print(f"NTC rows detected: {len(ntc_idx)}", flush=True)
if len(ntc_idx) == 0:
    known = set(tgt.keys())
    ntc_idx = np.where(~np.isin(base_genes, list(known)) & (np.char.str_len(base_genes.astype(str)) >= 4))[0][:50]
    print(f"fallback NTC rows: {len(ntc_idx)}", flush=True)

# var gene symbols / index (Ensembl -> symbol)
var_syms = list(adata.var.index)
if var_syms and str(var_syms[0]).startswith("ENSG"):
    import mygene
    mg = mygene.MyGeneInfo()
    q = mg.querymany(var_syms, scopes="ensembl.gene", fields="symbol", species="human", as_dataframe=True)
    symdict = q["symbol"].to_dict() if "symbol" in getattr(q, "columns", []) else {}
    sym_map = {s: (symdict[s] if (s in symdict and isinstance(symdict[s], str)) else None)
               for s in var_syms}
    eff_syms = [sym_map.get(s) for s in var_syms]
    print("mapped Ensembl->symbol via mygene", flush=True)
else:
    eff_syms = var_syms
sym2idx = {s: i for i, s in enumerate(eff_syms) if s}
print(f"var genes with symbol mapped: {len(sym2idx)}/{len(var_syms)}", flush=True)

# NTC baseline mean (read only NTC rows); sanitize non-finite (gemgroup Z can be inf
# for constant-variance genes -> division by zero) to neutral 0 (Z baseline)
X = adata.X
ntc_rows = X[ntc_idx].toarray() if hasattr(X[ntc_idx], "toarray") else np.asarray(X[ntc_idx])
ntc_rows = np.nan_to_num(ntc_rows, nan=0.0, posinf=0.0, neginf=0.0)
ntc_mean = ntc_rows.mean(axis=0)
print(f"NTC mean Z (sanity, should be ~0): min={ntc_mean.min():.3f} max={ntc_mean.max():.3f} "
      f"mean={ntc_mean.mean():.4f}", flush=True)

# candidate programs present as perturbations AND with enough targets
present_base = set(base_genes)
cand = [g for g in present_base if g in tgt and len([s for s in tgt[g] if s in sym2idx]) >= 5]
print(f"candidate programs (present + >=5 mapped targets): {len(cand)}", flush=True)

results = {"meta": {"dataset": "Replogle2022 K562_gwps normalized pseudo-bulk (CRISPRi)",
                    "n_perturbations": int(len(base_genes)),
                    "n_ntc": int(len(ntc_idx)),
                    "n_candidate_programs": len(cand)},
           "tfs": {}}

for tf in TFS:
    m = (base_genes == tf)
    n = int(m.sum())
    if n == 0:
        print(f"\n[{tf}] ABSENT in K562 (not perturbed) -> skipped", flush=True)
        results["tfs"][tf] = {"present": False, "reason": "not perturbed in K562 screen"}
        continue
    rows = X[np.where(m)[0]]
    rows = rows.toarray() if hasattr(rows, "toarray") else np.asarray(rows)
    rows = np.nan_to_num(rows, nan=0.0, posinf=0.0, neginf=0.0)
    sig = rows.mean(axis=0)               # query perturbation signature (Z, gemgroup-normalized)
    eff = sig - ntc_mean                  # baseline-corrected effect
    self_i = sym2idx.get(tf)
    self_z = float(eff[self_i]) if self_i is not None else None
    sz_disp = self_z if self_z is not None else float('nan')
    print(f"\n[{tf}] guides={n}  self-Z(CRISPRi knockdown of TF itself)={sz_disp:.3f} "
          f"(expect strongly negative)", flush=True)

    tg = [sym2idx[s] for s in tgt[tf] if s in sym2idx]
    nontg = [i for i in range(len(eff)) if i not in set(tg)]
    tg_z = eff[tg]
    nontg_z = eff[nontg]

    # (a) Mann-Whitney one-sided: targets shifted DOWN vs non-targets
    U, p_mwu = mannwhitneyu(tg_z, nontg_z, alternative="less")
    # (b) Fisher OR on "down" definition = bottom quartile of effect
    thr = np.percentile(eff, 25)
    a = int((tg_z < thr).sum()); b = len(tg_z) - a
    c = int((nontg_z < thr).sum()); d = len(nontg_z) - c
    OR, p_fisher = fisher_exact([[a, b], [c, d]], alternative="greater")
    mean_tg_z = float(tg_z.mean())

    # (c) rank among ALL candidate programs by mean(target Z) under this query signature
    prog_means = {}
    for c_ in cand:
        ct = [sym2idx[s] for s in tgt[c_] if s in sym2idx]
        if len(ct) < 5:
            continue
        prog_means[c_] = float(eff[ct].mean())
    all_means = np.array(list(prog_means.values()))
    self_mean = prog_means[tf]
    rank = int((all_means < self_mean).sum()) + 1      # 1 = most negative (most down)
    N = len(prog_means)
    pct = (N - rank) / N * 100.0                         # top percentile (most down)

    # cross-specificity: mean(target Z) of OTHER focal TFs under this query
    cross = {}
    for o in TFS:
        if o == tf or o not in tgt:
            continue
        ot = [sym2idx[s] for s in tgt[o] if s in sym2idx]
        if len(ot) >= 5:
            cross[o] = float(eff[ot].mean())

    results["tfs"][tf] = {
        "present": True, "n_guides": n, "self_z": self_z,
        "n_targets_mapped": len(tg), "n_targets_total": len(tgt[tf]),
        "mean_target_z": mean_tg_z,
        "mwu_U": float(U), "mwu_p_less": float(p_mwu),
        "down_thr_q25": float(thr),
        "fisher_OR_down": float(OR), "fisher_p": float(p_fisher),
        "down_targets": a, "n_targets": len(tg_z),
        "rank_among_programs": rank, "n_programs": N,
        "pct_most_down": float(pct),
        "cross_other_tf_meanZ": cross,
    }
    print(f"  n_targets={len(tg)} mean_target_Z={mean_tg_z:.4f} "
          f"MWU_p={p_mwu:.2e} Fisher_OR={OR:.2f} p={p_fisher:.2e} "
          f"rank={rank}/{N} ({pct:.2f}% most-down)", flush=True)
    print(f"  cross meanZ: {cross}", flush=True)

with open(OUT, "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2, ensure_ascii=False)
print("\nSAVED", OUT, flush=True)
