"""Generate Supplementary Figure S2: GSE197731 24h/48h cell-cell communication
+ cross-cohort (GSE174574 vs GSE197731_WT) 24h consistency at pathway/cell-type level.
Outputs PNG (manuscript embed) and SVG (editable, fonttype=none + pcolormesh).
"""
import os, re, ast
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.normpath(os.path.join(HERE, "..", "..", "rewiring_study", "figures"))
os.makedirs(FIGDIR, exist_ok=True)

csv = os.path.join(HERE, "cellchat_gse197731_all.csv")
txt = os.path.join(HERE, "cross_cohort_24h_consistency.txt")
df = pd.read_csv(csv)

# ---- Panel A: 38 significant WT 24h LR pairs, aggregated mean log2FC ----
sub = df[(df.study == "GSE197731_WT") & (df.transition == "24h") & (df.sig)]
agg = (sub.groupby(["ligand", "receiver" if False else "receptor"])
          .agg(log2FC=("log2FC", "mean"))
          .reset_index())
agg["pair"] = agg["ligand"] + "→" + agg["receptor"]
agg = agg.sort_values("log2FC")
pairs = agg["pair"].tolist()
vals = agg["log2FC"].values.astype(float)
n = len(vals)
mat = vals.reshape(-1, 1)

# ---- Parse cross-cohort consistency txt ----
def parse_list(s):
    s = s.strip()
    # convert np.str_('X') -> 'X' precisely, then literal_eval
    s = re.sub(r"np\.str_\('([^']*)'\)", r"'\1'", s)
    return ast.literal_eval(s)

p174 = p197 = None
ct174 = ct197 = None
with open(txt) as f:
    for line in f:
        if line.startswith("GSE174574_24h_pathways="):
            p174 = parse_list(line.split("=", 1)[1])
        elif line.startswith("GSE197731_24h_pathways="):
            p197 = parse_list(line.split("=", 1)[1])
        elif line.startswith("GSE174574_24h_ctpairs="):
            ct174 = parse_list(line.split("=", 1)[1])
        elif line.startswith("GSE197731_24h_ctpairs="):
            ct197 = parse_list(line.split("=", 1)[1])

shared_path = set(p174) & set(p197)
shared_ct = set(ct174) & set(ct197)

# ---- Build figure ----
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.size"] = 9

fig = plt.figure(figsize=(11, 8.2))
# Panel A: heatmap (pcolormesh -> vector in SVG)
axA = fig.add_axes([0.34, 0.08, 0.10, 0.84])
norm = TwoSlopeNorm(vmin=min(vals.min(), -0.1), vcenter=0.0, vmax=max(vals.max(), 0.1))
cmap = plt.cm.RdBu_r
mesh = axA.pcolormesh(mat, cmap=cmap, norm=norm, edgecolors="none")
axA.set_yticks(np.arange(n) + 0.5)
axA.set_yticklabels(pairs, fontsize=6)
axA.set_xticks([])
axA.set_ylabel("")
axA.set_title("")
axA.text(0.02, 0.92, "A. 24h WT LR pairs (log₂ Ipsil/Cont)", transform=axA.transAxes,
         ha="left", va="top", fontsize=9, fontweight="bold", color="#000000")
axA.invert_yaxis()
cb = fig.colorbar(mesh, ax=axA, fraction=0.9, pad=0.02, aspect=30)
cb.set_label("log₂FC", fontsize=8)

# Panel B: pathway convergence
axB = fig.add_axes([0.52, 0.45, 0.40, 0.45])
cats = ["GSE174574 24h", "GSE197731 WT 24h"]
counts = [len(p174), len(p197)]
bars = axB.bar(cats, counts, color=["#4C72B0", "#DD8452"], width=0.5)
for b, c in zip(bars, counts):
    axB.text(b.get_x() + b.get_width()/2, c + 0.3, str(c), ha="center", fontsize=9)
axB.set_ylabel("# pathways implicated")
axB.set_title(f"B. Pathway-level 24h cross-cohort\nShared: {sorted(shared_path)}  (Jaccard={len(shared_path)/len(set(p174)|set(p197)):.2f})", fontsize=9)
axB.set_ylim(0, max(counts) + 3)
# annotate shared
axB.text(0.5, 0.95, f"shared pathway = {', '.join(sorted(shared_path))}\n(EGF unique to cohort 1)",
         transform=axB.transAxes, ha="center", va="top", fontsize=7.5, color="#333333")

# Panel C: cell-type convergence
axC = fig.add_axes([0.52, 0.08, 0.40, 0.25])
cats2 = ["GSE174574 24h", "GSE197731 WT 24h"]
counts2 = [len(ct174), len(ct197)]
bars2 = axC.bar(cats2, counts2, color=["#55A868", "#C44E52"], width=0.5)
for b, c in zip(bars2, counts2):
    axC.text(b.get_x() + b.get_width()/2, c + 0.4, str(c), ha="center", fontsize=9)
axC.set_ylabel("# sender→receiver\ncell-type pairs")
axC.set_title(f"C. Cell-type-level 24h cross-cohort\nShared: {[tuple(x) for x in shared_ct]}", fontsize=9)
axC.set_ylim(0, max(counts2) + 4)
axC.text(0.5, -0.45, "Pericyte→Pericyte shared across both 24h cohorts",
         transform=axC.transAxes, ha="center", va="top", fontsize=7.5, color="#333333")

fig.suptitle("Supplementary Fig. S2 — GSE197731 24h/48h communication & cross-cohort 24h consistency",
             fontsize=11, fontweight="bold", y=1.0)

png = os.path.join(FIGDIR, "supp_fig_s2_cellchat_24h_48h.png")
svg = os.path.join(FIGDIR, "supp_fig_s2_cellchat_24h_48h.svg")
fig.savefig(png, dpi=150, bbox_inches="tight")
fig.savefig(svg, bbox_inches="tight")
print("WROTE", png, "and", svg)
print("Panel A pairs:", n, "shared pathways:", sorted(shared_path), "shared ct:", [tuple(x) for x in shared_ct])
