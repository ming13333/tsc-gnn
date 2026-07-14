"""生成 Fig 8: 跨转变一致的 重布线 LR 对 热图 (仅显著记录)。
读 cellchat_rewiring_sig.csv, 不重算 MWU, 秒级。"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

df = pd.read_csv(os.path.join(HERE, "cellchat_rewiring_sig.csv"))
sig = df[df.sig]

# 跨转变一致: 同一 LR 对在 >=2 个 transition 显著
pc = sig.groupby(["ligand", "receptor"]).size()
cons = pc[pc >= 2].sort_values(ascending=False)
cons_set = set(cons.index)
sub = sig[sig.apply(lambda r: (r.ligand, r.receptor) in cons_set, axis=1)].copy()
sub["pair"] = sub.ligand + " → " + sub.receptor
sub["trans"] = sub.apply(lambda r: "24h" if r.study == "GSE174574" else r.transition, axis=1)

heat = sub.pivot_table(index="pair", columns="trans", values="log2FC", aggfunc="first")
for c in ["24h", "2d", "14d"]:
    if c not in heat.columns:
        heat[c] = np.nan
heat = heat[["24h", "2d", "14d"]]

fig, ax = plt.subplots(figsize=(5.5, max(5, len(heat) * 0.34)))
im = ax.imshow(heat.values, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
ax.set_xticks(range(3))
ax.set_xticklabels(heat.columns)
ax.set_yticks(range(len(heat)))
ax.set_yticklabels(heat.index, fontsize=7)
ax.set_title("Cell–cell communication (LR) rewiring\nacross stroke transitions", fontsize=10)
ax.set_ylabel("Ligand → Receptor")
fig.colorbar(im, ax=ax, label="log₂(MCAO / sham) communication score")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "cellchat_rewiring_heatmap.png"), dpi=150)
print(f"Fig 8 saved: {heat.shape[0]} consistent LR pairs × 3 transitions")
