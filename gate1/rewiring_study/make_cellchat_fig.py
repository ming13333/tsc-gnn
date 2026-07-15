"""生成 Fig 8: 跨转变一致的 重布线 LR 对 热图 (仅显著记录)。
读 cellchat_rewiring_sig.csv, 不重算 MWU, 秒级。"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# SVG 可编辑化: 文字保留 <text> 可改字; 用 pcolormesh 让单元格为矢量 <path>
# (imshow 在 SVG 里会被栅格成一张 base64 位图, 无法编辑)
plt.rcParams["svg.fonttype"] = "none"

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

# 全矢量: 主图 pcolormesh + 矢量色条 (避免 matplotlib 把 colorbar 栅格成 PNG)
ny, nx = heat.shape
gx = np.arange(nx + 1)
gy = np.arange(ny + 1)

from matplotlib.gridspec import GridSpec
fig = plt.figure(figsize=(6.2, max(5, len(heat) * 0.34)))
gs = GridSpec(1, 2, width_ratios=[1, 0.045], wspace=0.08,
              left=0.30, right=0.88, top=0.90, bottom=0.08)
ax = fig.add_subplot(gs[0])
cax = fig.add_subplot(gs[1])

mesh = ax.pcolormesh(gx, gy, heat.values, cmap="RdBu_r", vmin=-3, vmax=3,
                     edgecolors="none", shading="flat")
ax.set_xlim(0, nx)
ax.set_ylim(0, ny)
ax.set_xticks(np.arange(nx) + 0.5)
ax.set_xticklabels(heat.columns)
ax.set_yticks(np.arange(ny) + 0.5)
ax.set_yticklabels(heat.index, fontsize=7)
ax.set_title("Cell–cell communication (LR) rewiring\nacross stroke transitions", fontsize=10)
ax.set_ylabel("Ligand → Receptor")

# 矢量色条: 用 pcolormesh 画渐变, 不栅格化
grad = np.tile(np.linspace(3, -3, 256).reshape(-1, 1), (1, 2))
cax.pcolormesh([0, 1, 2], np.arange(257), grad, cmap="RdBu_r", vmin=-3, vmax=3,
               shading="flat", edgecolors="none")
cax.set_xlim(0, 1)
cax.set_ylim(0, 256)
cax.invert_yaxis()
cax.set_xticks([])
cax.set_yticks([0, 64, 128, 192, 255])
cax.set_yticklabels(["3", "1.5", "0", "−1.5", "−3"], fontsize=7)
cax.set_ylabel("log₂(MCAO / sham) communication score", fontsize=8)
fig.savefig(os.path.join(OUT, "cellchat_rewiring_heatmap.png"), dpi=150)
fig.savefig(os.path.join(OUT, "cellchat_rewiring_heatmap.svg"))
print(f"Fig 8 saved: {heat.shape[0]} consistent LR pairs × 3 transitions (PNG+SVG, fully vector)")
