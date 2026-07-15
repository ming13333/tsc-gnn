#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
regen_l5b_json.py — restore l5b_sigcom_result.json from the VERIFIED successful
run that is preserved in l5b_run.log (first [done] block, lines 1-37).

Background (CONFORMANCE audit B2, 2026-07-14):
  The original l5b_sigcom.py ran successfully once and wrote a fully-populated
  JSON (the numbers the manuscript reports: GATA2 OE rank 3 / 33,782, p = 1.4e-5).
  A later re-run hit an HTTP 500 on the GATA2 Overexpression query, and because
  the script wrote its output *unconditionally*, the second (failed) run
  OVERWROTE the good JSON with an empty shell. The numbers in the manuscript are
  correct (they came from the first run / L5b_SigCom_报告); only the machine-
  readable artifact was corrupted.

This script does NOT call the API again (the endpoint was returning 500). It
reconstructs the genuine computed result from the authoritative log so the JSON
is real and re-readable. When the SigCom LINCS API is healthy, re-run
l5b_sigcom.py (now hardened to never overwrite a good result with a failed one)
to produce a fresh artifact.

Usage:  python regen_l5b_json.py [target_dir]
  target_dir defaults to this script's own directory.
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
# allow an explicit target dir so the same script can regenerate either copy
if len(sys.argv) > 1:
    HERE = sys.argv[1]
OUT = os.path.join(HERE, "l5b_sigcom_result.json")

# ---- verified numbers copied verbatim from l5b_run.log successful block ----
# Each non-None entry: (rank, type, p).  type REV -> p_down ; type MIM -> p_up.
CRISPR_KO = {
    "SOX10|SOX10": (139333, "REV", 7.4e-02),
    "SOX10|CEBPB": (139116, "REV", 4.5e-04),
    "SOX10|GATA2": (139217, "REV", 1.7e-03),
    "CEBPB|SOX10": None,
    "CEBPB|CEBPB": (139667, "REV", 9.2e-04),
    "CEBPB|GATA2": None,
    "GATA2|SOX10": (139136, "REV", 3.7e-01),
    "GATA2|CEBPB": ( 1215, "MIM", 2.8e-03),
    "GATA2|GATA2": (138952, "REV", 3.9e-02),
}
OE = {
    "SOX10|SOX10": (32109, "REV", 3.8e-02),
    "SOX10|CEBPB": None,
    "SOX10|GATA2": (  971, "MIM", 2.0e-03),
    "CEBPB|SOX10": (33444, "REV", 2.8e-01),
    "CEBPB|CEBPB": None,
    "CEBPB|GATA2": ( 1139, "MIM", 1.9e-04),
    "GATA2|SOX10": (32101, "REV", 1.0e-01),
    "GATA2|CEBPB": None,
    "GATA2|GATA2": (    3, "MIM", 1.4e-05),   # <-- headline result in the manuscript
}

def mk_matrix(d):
    out = {}
    for k, v in d.items():
        if v is None:
            out[k] = []
        else:
            rank, typ, p = v
            if typ == "REV":
                out[k] = [{"rank": rank, "type": "reversers",
                           "p_up": None, "p_down": p,
                           "fdr_up": None, "fdr_down": None,
                           "z_up": None, "z_down": None, "z_sum": None}]
            else:
                out[k] = [{"rank": rank, "type": "mimickers",
                           "p_up": p, "p_down": None,
                           "fdr_up": None, "fdr_down": None,
                           "z_up": None, "z_down": None, "z_sum": None}]
    return out

results = {
    "CRISPR_KO": {
        "matrix": mk_matrix(CRISPR_KO),
        "n_rev": {"SOX10": 70563, "CEBPB": 70890, "GATA2": 70298},
        "n_mim": {"SOX10": 70040, "CEBPB": 69713, "GATA2": 70305},
        "max_rank": {"SOX10": 140603, "CEBPB": 140603, "GATA2": 140603},
    },
    "Overexpression": {
        "matrix": mk_matrix(OE),
        "n_rev": {"SOX10": 17059, "CEBPB": 17445, "GATA2": 17169},
        "n_mim": {"SOX10": 16723, "CEBPB": 16337, "GATA2": 16613},
        "max_rank": {"SOX10": 33782, "CEBPB": 33782, "GATA2": 33782},
    },
}

doc = {
    "targets_n": {"SOX10": 322, "CEBPB": 589, "GATA2": 5370},
    "ent_targets_n": {"SOX10": 321, "CEBPB": 500, "GATA2": 500},
    "own_uuids_n": {
        "CRISPR_KO": {"SOX10": 20, "CEBPB": 20, "GATA2": 44},
        "Overexpression": {"SOX10": 14, "CEBPB": 0, "GATA2": 15},
    },
    "results": results,
    "_reconstructed_from": ("l5b_run.log successful block (lines 1-37); the second "
                            "rerun overwrote this JSON with an empty shell (HTTP 500)."),
    "_reconstruction_date": "2026-07-14",
}

# ---- assertions against the manuscript claims ----
g = doc["results"]["Overexpression"]["matrix"]["GATA2|GATA2"][0]
assert g["rank"] == 3, g
assert abs(g["p_up"] - 1.4e-05) < 1e-9, g
assert doc["results"]["Overexpression"]["max_rank"]["GATA2"] == 33782
# self-specificity check: the cross (GATA2->SOX10) is far less significant
cross = doc["results"]["Overexpression"]["matrix"]["GATA2|SOX10"][0]
assert cross["rank"] == 32101 and cross["type"] == "reversers"

with open(OUT, "w", encoding="utf-8") as f:
    json.dump(doc, f, ensure_ascii=False, indent=2)
print(f"wrote {OUT}")
print(f"  GATA2 OE -> GATA2 : rank {g['rank']} / 33782, p = {g['p_up']:.1e} "
      f"(mimicker, self-specific)")
print(f"  GATA2 OE -> SOX10 : rank {cross['rank']} (cross, not self-specific)")
n_filled = sum(1 for lib in results.values()
               for m in lib["matrix"].values() if m)
print(f"  populated matrix cells: {n_filled} / 18")
