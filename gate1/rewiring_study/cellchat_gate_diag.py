# cellchat_gate_diag.py  (轻量: 仅 2 次 vstack, 逐对门控+安全MWU, 内存安全)
import os, json
import numpy as np
import cellchat_py as cc
from scipy.stats import mannwhitneyu
from cellchat_lrdb import LR_DB

PAIRS = sorted({(l, r) for l, r, _ in LR_DB})

M17 = {
    "Microglia": ["Cx3cr1", "Tmem119", "P2ry12", "Siglech", "Aif1", "Ctss"],
    "CAM": ["Mrc1", "Cd163", "Lyve1", "Timd4", "Fcgr3", "Mafb"],
    "MdC": ["Ly6c2", "Ccr2", "S100a8", "S100a9", "Cxcl9", "Itgam"],
    "Astrocyte": ["Gfap", "Aqp4", "Slc1a2", "Aldh1l1", "Sparcl1", "Gja1"],
    "Oligodendrocyte": ["Mbp", "Mog", "Olig2", "Plp1", "Cnp"],
    "OPC": ["Pdgfra", "Cspg4", "Pcdh15", "Bcan", "Cspg5"],
    "Neuron": ["Snap25", "Syt1", "Rbfox3", "Slc17a7", "Gabbr1"],
    "aEC": ["Cldn5", "Kdr", "Gja5", "Ephb4", "Hey1", "Dll4"],
    "vEC": ["Cldn5", "Kdr", "Nr2f2", "Vwf", "Sele", "Icam1"],
    "capEC": ["Cldn5", "Kdr", "Slc38a5", "Tspan7", "Car4", "Gpihbp1"],
    "SMC": ["Acta2", "Tagln", "Myh11", "Cnn1"],
    "Pericyte": ["Pdgfrb", "Rgs5", "Abcc9", "Kcnj8", "Des"],
    "FB": ["Col1a1", "Col1a2", "Dcn", "Lum", "Pi16", "Scara5"],
    "CPC": ["Cldn5", "Kdr", "Ttr", "Aqp1", "Mfsd2a"],
    "Ependymal": ["Foxj1", "Ccdc153", "Cd24a", "Dnai1"],
    "Tcell": ["Cd3g", "Cd3e", "Cd4", "Cd8a", "Trbc1"],
    "Bcell": ["Cd79a", "Ms4a1", "Cd19", "Cd22"],
}


def build_pools(markers):
    cc.CELLTYPE_MARKERS = markers
    samples = cc.load_gse174()
    cc._align_genes_per_study(samples)
    sham = [s for s in samples if s["condition"] == "sham"]
    mc = [s for s in samples if s["condition"] == "MCAO"]
    sham_X = np.vstack([s["X"] for s in sham]); sham_ct = np.concatenate([s["cell_type"] for s in sham])
    mc_X = np.vstack([s["X"] for s in mc]); mc_ct = np.concatenate([s["cell_type"] for s in mc])
    return sham_X, sham_ct, sham[0]["genes"], mc_X, mc_ct


def evaluate(markers, label):
    sham_X, sham_ct, sham_g, mc_X, mc_ct = build_pools(markers)
    gi = {g: i for i, g in enumerate(sham_g)}
    out = {}
    for mf in (0.05, 0.10, 0.15):
        cc.MIN_FRAC = mf
        cand = sig = 0
        uniq = set()
        for lig, rec in PAIRS:
            if lig not in gi or rec not in gi:
                continue
            li, ri = gi[lig], gi[rec]
            found = None
            for s in np.unique(sham_ct):
                for r in np.unique(sham_ct):
                    si_s = np.where(sham_ct == s)[0]; ri_s = np.where(sham_ct == r)[0]
                    si_m = np.where(mc_ct == s)[0]; ri_m = np.where(mc_ct == r)[0]
                    if min(len(si_s), len(ri_s), len(si_m), len(ri_m)) < cc.MIN_CELLS:
                        continue
                    fLs = np.mean(sham_X[si_s, li] > 0); fRs = np.mean(sham_X[ri_s, ri] > 0)
                    fLm = np.mean(mc_X[si_m, li] > 0); fRm = np.mean(mc_X[ri_m, ri] > 0)
                    if min(fLs, fRs, fLm, fRm) < mf:
                        continue
                    score_s = sham_X[si_s, li].mean() * sham_X[ri_s, ri].mean()
                    score_m = mc_X[si_m, li].mean() * mc_X[ri_m, ri].mean()
                    log2fc = np.log2((score_m + cc.EPS) / (score_s + cc.EPS))
                    try:
                        pL = mannwhitneyu(mc_X[si_m, li], sham_X[si_s, li], alternative="two-sided").pvalue
                        pR = mannwhitneyu(mc_X[ri_m, ri], sham_X[ri_s, ri], alternative="two-sided").pvalue
                    except Exception:
                        pL = pR = 1.0
                    is_sig = (min(pL, pR) < 0.05) and (abs(log2fc) >= 0.5)
                    if found is None or abs(log2fc) > abs(found[0]):
                        found = (log2fc, is_sig)
            if found is not None:
                cand += 1; uniq.add((lig, rec))
                if found[1]:
                    sig += 1
        out[mf] = {"candidate": int(cand), "significant": int(sig), "uniq": len(uniq)}
        print(f"{label:<8} MIN_FRAC={mf}: candidate={cand}, significant={sig}, uniq={len(uniq)}", flush=True)
    return out


def main():
    res = {}
    res["11-type"] = evaluate(cc.CELLTYPE_MARKERS, "11-type")
    res["17-type"] = evaluate(M17, "17-type")
    # 17-type@0.10 细胞类型分布
    cc.CELLTYPE_MARKERS = M17
    samp = cc.load_gse174()
    dist = {}
    for s in samp:
        vals, counts = np.unique(s["cell_type"], return_counts=True)
        dist[s["condition"]] = {str(k): int(v) for k, v in zip(vals.tolist(), counts.tolist())}
    res["17-type_distribution@0.10"] = dist
    json.dump(res, open("cellchat_annot_result.json", "w"), indent=2)
    print("=== DONE -> cellchat_annot_result.json ===", flush=True)


if __name__ == "__main__":
    main()
