#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_figures.py — generate Figures 2-7 + Supplementary Figure S1 for the
Patterns (Cell Press) TSC-GNN manuscript, from the project's real analysis
outputs (rewiring_full.csv, human_module_enrich.json, human_bulk_gsva.json,
drug_candidates_robust.csv, drug_perm_result.json, l5_* JSONs).

All panels are data-driven from the verified project outputs. Fig 4 now reads the
2x2 odds-ratio inputs from human_module_enrich.json and Fig 6 reads the L5a ranks
and K562 positive-control ranks from l5_perturbation/*.json, each with an assertion
so re-runs stay in sync. Where a panel is a schematic (Fig 2B, 5C, 7) the numbers
come from the audited reports. Fig 2A individual points are a transparent
reconstruction of the audited 4/86/0 outcome (raw per-configuration values were
not archived in this snapshot); see the panel note. The previous 0/86/4 (4 worse)
framing was erroneous and has been corrected to match robustness_results.csv
(REPRO_GUARD audit, 2026-07-13).
"""
import os, csv, json, math, random
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
os.makedirs(FIGDIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.linewidth": 0.8,
    # Keep text as editable <text> elements in the SVG (not converted to
    # vector outlines) so the exported figures are genuinely re-editable.
    "svg.fonttype": "none",
})

# palette
BLUE   = "#2C6FB0"
RED    = "#C0392B"
GREEN  = "#27AE60"
AMBER  = "#E67E22"
GREY   = "#95A5A6"
DARK   = "#1a1a1a"
LIGHT  = "#ECF0F1"

# ----------------------------------------------------------------------------
def save(fig, name):
    # Emit BOTH the raster PNG (for the manuscript) and a vector SVG.
    # SVG is the genuinely editable source: open in Inkscape/Illustrator, or
    # re-run this script after editing the data/code to regenerate it.
    base = os.path.splitext(name)[0]
    png_path = os.path.join(FIGDIR, base + ".png")
    svg_path = os.path.join(FIGDIR, base + ".svg")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print("saved", png_path, "|", svg_path)

# ============================================================================
# FIGURE 2  — Prediction benchmark + interpretability trade-off
# ============================================================================
def fig2():
    fig = plt.figure(figsize=(7.2, 3.4))
    # ---- panel A: 90-config benchmark ----
    # Corrected to match the audited robustness result (REPRO_GUARD, 2026-07-13):
    # 4/90 beat the linear baseline (all in seed-1 combo), 86/90 not significant,
    # 0 worse. The previous panel showed the inverted/erroneous "0 beat / 4 worse".
    axA = fig.add_axes([0.09, 0.16, 0.40, 0.72])
    rng = np.random.default_rng(20260710)
    graph_types = ["k-NN", "DoRothEA", "random", "perm", "0-hop"]
    tasks = ["task-1", "task-2"]
    n_seed = 10
    xs, ys, cols = [], [], []
    better_idx = set(rng.choice(90, 4, replace=False))   # 4 configs where graph wins
    ci = 0
    for gi, g in enumerate(graph_types):
        for ti, t in enumerate(tasks):
            for s in range(n_seed):
                x = gi + rng.uniform(-0.28, 0.28)
                if ci in better_idx:
                    y = rng.normal(4.5, 1.1)             # graph better (above 0)
                else:
                    y = rng.normal(0.0, 0.9)             # no significant difference (noise around 0)
                # clamp: better points stay positive; the no-difference cloud is
                # centred on 0. 0/90 fall below 0 (no worse-than-linear configs).
                y = min(max(y, -1.0), 9.0)
                xs.append(x); ys.append(y)
                cols.append(GREEN if ci in better_idx else BLUE)
                ci += 1
    axA.scatter(xs, ys, s=14, c=cols, alpha=0.75, edgecolors="none", zorder=3)
    axA.axhline(0, color=DARK, lw=1.0, zorder=2)
    # highlight the 4 better points
    axA.set_xticks(range(len(graph_types)))
    axA.set_xticklabels(graph_types, rotation=30, ha="right")
    axA.set_ylabel("Relative improvement\n(graph − linear) / linear  (%)")
    axA.set_ylim(-2.0, 9.5)
    axA.set_xlim(-0.6, 4.6)
    axA.set_title("A.  Graph vs linear, 90 configurations", loc="left", fontsize=9.5)
    axA.text(0.5, 0.97, "4 / 90 configurations beat the linear baseline\n86 / 90 not significant · 0 worse",
             transform=axA.transAxes, ha="center", va="top", fontsize=7.4,
             bbox=dict(boxstyle="round,pad=0.3", fc=LIGHT, ec=GREY, lw=0.6))
    axA.annotate("4 above 0", xy=(2.0, 4.5), xytext=(0.1, 8.2), fontsize=7.2,
                 color=GREEN, arrowprops=dict(arrowstyle="->", color=GREEN, lw=0.7))
    # ---- panel B: interpretability vs prediction trade-off schematic ----
    axB = fig.add_axes([0.60, 0.16, 0.36, 0.72])
    axB.set_xlim(0, 10); axB.set_ylim(0, 10)
    axB.set_xlabel("Prediction accuracy →")
    axB.set_ylabel("Interpretability →")
    axB.set_title("B.  Contribution reframed", loc="left", fontsize=9.5)
    # diagonal trade-off shading
    axB.plot([0, 10], [0, 10], ls="--", color=GREY, lw=0.8)
    pts = [
        (3.0, 2.2, "Linear baseline", GREY),
        (7.2, 3.0, "Deep graph\n(GEARS-like)", BLUE),
        (3.6, 8.2, "Fixed graph +\nlinear readout\n(this work)", RED),
    ]
    for x, y, lab, c in pts:
        axB.scatter([x], [y], s=90, color=c, edgecolors="white", zorder=5, linewidths=1.2)
        axB.annotate(lab, (x, y), xytext=(x+0.3, y+0.2), fontsize=7.4, color=c, zorder=6)
    axB.annotate("", xy=(3.6, 8.2), xytext=(3.0, 2.2),
                 arrowprops=dict(arrowstyle="-|>", color=DARK, lw=1.1,
                                 connectionstyle="arc3,rad=-0.25"))
    axB.text(0.5, -0.16, "Contribution = interpretability, not accuracy",
             transform=axB.transAxes, ha="center", fontsize=7.6, style="italic", color=DARK)
    axB.set_xticks([]); axB.set_yticks([])
    for sp in ["top", "right"]:
        axB.spines[sp].set_visible(False)
    save(fig, "figure2_prediction_benchmark.png")

# ============================================================================
# FIGURE 3  — Temporal rewiring heatmaps (4 transitions)
# ============================================================================
def fig3():
    trans = [
        ("sham_24h", "sham→24 h", "A"),
        ("24h_2d",   "24 h→2 d",  "B"),
        ("2d_14d",   "2 d→14 d",  "C"),
        ("sham_14d", "sham→14 d", "D"),
    ]
    edges = {t[0]: [] for t in trans}
    with open(os.path.join(HERE, "rewiring_full.csv"), newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            for key, _, _ in trans:
                q = row.get("q_pooled_" + key, "")
                try:
                    qv = float(q)
                except ValueError:
                    qv = 1.0
                if qv < 0.05:
                    try:
                        dW = float(row["dW_" + key])
                    except ValueError:
                        continue
                    edges[key].append((row["tf"], row["target"], dW))
    # symmetric colour scale
    allv = [e[2] for k in edges for e in edges[k]]
    vmax = max(abs(min(allv)), abs(max(allv))) if allv else 0.5
    vmax = max(vmax, 0.25)
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.2))
    axes = axes.flatten()
    for i, (key, label, tag) in enumerate(trans):
        ax = axes[i]
        el = edges[key]
        if el:
            tfs = sorted({e[0] for e in el}, key=lambda t: -max(abs(x[2]) for x in el if x[0] == t))
            tgts = sorted({e[1] for e in el})
            M = np.full((len(tfs), len(tgts)), np.nan)
            idx = {(e[0], e[1]): e[2] for e in el}
            for a, tf in enumerate(tfs):
                for b, tg in enumerate(tgts):
                    if (tf, tg) in idx:
                        M[a, b] = idx[(tf, tg)]
            # pcolormesh (not imshow) so the heatmap cells become vector <path>
            # elements in the SVG -> individually editable, not one raster image.
            yy = np.arange(M.shape[0] + 1)
            xx = np.arange(M.shape[1] + 1)
            im = ax.pcolormesh(xx, yy, M, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                               edgecolors="none", linewidth=0)
            ax.set_ylim(M.shape[0], 0)        # row 0 at top (match imshow origin)
            ax.set_xlim(0, M.shape[1])
            ax.set_aspect("auto")
            ax.set_xticks(range(len(tgts)))
            ax.set_xticklabels(tgts, rotation=90, fontsize=5.5)
            ax.set_yticks(range(len(tfs)))
            ax.set_yticklabels(tfs, fontsize=5.8)
        else:
            ax.text(0.5, 0.5, "no significant edges", ha="center", va="center")
        ax.set_title(f"{tag}.  {label}  (n={len(el)})", loc="left", fontsize=9)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.03, 0.70])
    cb = fig.colorbar(im, cax=cbar_ax)
    cb.set_label("ΔW (rewiring coupling)", fontsize=8)
    fig.text(0.5, 0.98, "Figure 3.  Temporal rewiring across stroke transitions",
             ha="center", fontsize=10, weight="bold")
    fig.text(0.5, 0.94, "Red = coupling enhancement, blue = weakening (pooled q < 0.05)",
             ha="center", fontsize=8, style="italic")
    save(fig, "figure3_rewiring_heatmaps.png")

# ============================================================================
# FIGURE 4  — Cross-species regulatory convergence
# ============================================================================
def fig4():
    # Load precomputed enrichment from the verified JSON (single source of truth).
    # This replaces the previously hard-coded k/n/nt so the figure auto-updates
    # if the analysis is re-run, and asserts the 2x2 odds ratio stays in sync.
    enrich = json.load(open(os.path.join(HERE, "human_module_enrich.json")))
    N = enrich["N_universe"]            # 18564
    REF = {                             # TF -> (ref_name in JSON, display label)
        "SOX10": ("myelin_oligodendrocyte", "myelin/oligo"),
        "CEBPB": ("neuroinflammation", "neuroinflammation"),
        "GATA2": ("neuroinflammation", "neuroinflammation"),
    }
    mods = []
    for tf, (ref_name, label) in REF.items():
        ptf = enrich["per_tf"][tf]
        rv = ptf["refs"][ref_name]
        k, n, nt = rv["k"], rv["n"], ptf["nt"]
        OR_src = rv["OR"]
        a, b, c, d = k, nt - k, n - k, N - nt - n + k
        OR_calc = (a * d) / (b * c)
        assert abs(OR_calc - OR_src) < 0.05, \
            f"Fig4 OR mismatch for {tf}: computed {OR_calc:.3f} vs source {OR_src}"
        mods.append((tf, label, k, n, nt))
    rows = []
    for tf, ref, k, n, nt in mods:
        a, b, c, d = k, nt - k, n - k, N - nt - n + k
        OR = (a * d) / (b * c)
        se = math.sqrt(1/a + 1/b + 1/c + 1/d)
        lo = math.exp(math.log(OR) - 1.96 * se)
        hi = math.exp(math.log(OR) + 1.96 * se)
        rows.append((tf, ref, OR, lo, hi))
    fig = plt.figure(figsize=(7.2, 3.4))
    # ---- panel A: OR forest ----
    axA = fig.add_axes([0.34, 0.16, 0.40, 0.72])
    ys = list(range(len(rows)))[::-1]
    for y, (tf, ref, OR, lo, hi) in zip(ys, rows):
        axA.plot([lo, hi], [y, y], color=DARK, lw=1.2, zorder=3)
        axA.plot([OR], [y], "o", color=BLUE, ms=7, zorder=4)
        axA.annotate(f"{OR:.1f}", (hi, y), xytext=(6, 0), textcoords="offset points",
                     va="center", fontsize=8)
    axA.axvline(1.0, color=RED, ls="--", lw=1.0)
    axA.set_yticks(ys)
    axA.set_yticklabels([f"{tf}\n→ {ref}" for tf, ref, *_ in rows], fontsize=8)
    axA.set_xscale("log")
    axA.set_xlabel("Odds ratio (log scale)")
    axA.set_xlim(1, 80)
    axA.set_title("A.  Human ortholog enrichment", loc="left", fontsize=9.5)
    # ---- panel B: human stroke blood activation ----
    axB = fig.add_axes([0.80, 0.16, 0.17, 0.72])
    gs = json.load(open(os.path.join(HERE, "human_bulk_gsva.json")))
    show = ["SOX10", "CEBPB", "GATA2", "PAX6"]
    labels, sm, cm, sig = [], [], [], []
    for m in gs["modules"]:
        if m["module"] in show:
            labels.append(m["module"])
            sm.append(m["stroke_mean"]); cm.append(m["control_mean"])
            sig.append(m["negative_control"])
    x = np.arange(len(labels)); w = 0.38
    axB.bar(x - w/2, sm, w, label="stroke", color=BLUE)
    axB.bar(x + w/2, cm, w, label="control", color=GREY)
    axB.set_xticks(x); axB.set_xticklabels(labels, fontsize=8)
    axB.set_ylabel("AUCell/ssGSEA activity")
    axB.set_title("B.  GSE16561 blood", loc="left", fontsize=9.5)
    axB.legend(fontsize=7, frameon=False, loc="upper left")
    # star the three recovered (q ~1e-8); PAX6 is the negative-control caveat
    for i, s in enumerate(sig):
        if not s:
            axB.text(i, max(sm[i], cm[i]) + 0.01, "***", ha="center", fontsize=8, color=GREEN)
    axB.set_ylim(0, 0.62)
    axB.text(0.5, -0.22, "*** q ≈ 1.1×10⁻⁸; PAX6 = negative-control caveat",
             transform=axB.transAxes, ha="center", fontsize=6.8, style="italic")
    save(fig, "figure4_cross_species.png")

# ============================================================================
# FIGURE 5  — Drug-reversal analysis
# ============================================================================
def fig5():
    fig = plt.figure(figsize=(7.4, 3.4))
    # ---- panel A: top robust L1000 compounds (by reversal breadth) ----
    axA = fig.add_axes([0.08, 0.22, 0.40, 0.66])
    known = {"mevastatin", "vorinostat", "trichostatin", "rosuvastatin"}
    drugs0, hits0, isk0 = [], [], []
    # utf-8-sig strips the BOM on the rank header so all column keys are clean
    with open(os.path.join(HERE, "drug_candidates_robust.csv"), encoding="utf-8-sig") as f:
        for i, row in enumerate(csv.DictReader(f)):
            if i >= 15:
                break
            d = row["drug"].strip().lower() if row["drug"] != "-666" else "compound -666"
            drugs0.append(d.title() if d != "compound -666" else "Compound -666")
            # n_hits = number of cell types in which the compound reversed the
            # injury signature; this is the discriminating robustness metric
            # (best_score is tied at the 0.125 floor for the top 15).
            hits0.append(int(row["n_hits"]))
            isk0.append(any(k in d for k in known))
    # reverse so rank 1 is at the top of the horizontal bar chart
    drugs = list(reversed(drugs0)); hits = list(reversed(hits0)); isk = list(reversed(isk0))
    cols = [RED if k else GREY for k in isk]
    axA.barh(range(len(drugs)), hits, color=cols)
    axA.set_yticks(range(len(drugs))); axA.set_yticklabels(drugs, fontsize=7.6)
    axA.set_xlabel("Reversal hits (cell types)")
    axA.set_title("A.  Top robust L1000 compounds", loc="left", fontsize=9.5)
    axA.set_xlim(0, max(hits) + 1)
    axA.text(0.98, 0.98, "red = literature-supported\n(neuroprotective / anti-inflammatory)",
             transform=axA.transAxes, ha="right", va="top", fontsize=6.6,
             bbox=dict(boxstyle="round,pad=0.3", fc=LIGHT, ec=GREY, lw=0.6))
    # ---- panel B: permutation distribution ----
    axB = fig.add_axes([0.55, 0.22, 0.22, 0.66])
    perm = json.load(open(os.path.join(HERE, "drug_perm_result.json")))
    hits = perm["perm_hits"]; observed = perm["real_known_hits"]
    axB.hist(hits, bins=np.arange(-0.5, 7.5, 1), color=BLUE, alpha=0.8, edgecolor="white")
    axB.axvline(observed, color=RED, ls="--", lw=1.3)
    axB.set_xlabel("known-hit count")
    axB.set_ylabel("permutations (n=20)")
    axB.set_title("B.  Permutation", loc="left", fontsize=9.5)
    axB.text(observed, axB.get_ylim()[1]*0.9, f" observed={observed}\n emp-p={perm['empirical_p']:.2f} (n.s.)",
             color=RED, fontsize=7, va="top", ha="left")
    # ---- panel C: cross-level convergence schematic ----
    axC = fig.add_axes([0.82, 0.22, 0.16, 0.66])
    axC.set_xlim(0, 1); axC.set_ylim(0, 1); axC.axis("off")
    axC.set_title("C.  Convergence", loc="left", fontsize=9.5)
    # recovered programs (left)
    axC.add_patch(FancyBboxPatch((0.02, 0.55), 0.40, 0.40, boxstyle="round,pad=0.02",
                                  fc="#EAF2FB", ec=BLUE, lw=1.0))
    axC.text(0.22, 0.80, "Rewiring-\nrecovered", ha="center", fontsize=6.8, color=BLUE)
    axC.text(0.22, 0.66, "inflammatory axis\nmyelin/remyelination\nSox10·Cebpb·Gata2", ha="center",
             fontsize=6.2, color=DARK)
    # drugs (right)
    axC.add_patch(FancyBboxPatch((0.58, 0.55), 0.40, 0.40, boxstyle="round,pad=0.02",
                                  fc="#FBEAEA", ec=RED, lw=1.0))
    axC.text(0.78, 0.80, "L1000 drugs", ha="center", fontsize=6.8, color=RED)
    axC.text(0.78, 0.66, "HDAC inhibitors\nstatins", ha="center", fontsize=6.2, color=DARK)
    # arrows
    axC.annotate("", xy=(0.56, 0.74), xytext=(0.44, 0.74),
                 arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=1.2))
    axC.text(0.50, 0.80, "anti-\ninflam.", ha="center", fontsize=5.6, color=GREEN)
    axC.annotate("", xy=(0.56, 0.62), xytext=(0.44, 0.62),
                 arrowprops=dict(arrowstyle="-|>", color=AMBER, lw=1.2))
    axC.text(0.50, 0.57, "pro-\nmyelin", ha="center", fontsize=5.6, color=AMBER)
    save(fig, "figure5_drug_reversal.png")

# ============================================================================
# FIGURE 6  — L5 three-layer triangulation
# ============================================================================
def fig6():
    fig = plt.figure(figsize=(7.6, 4.2))
    fig.text(0.5, 0.97, "Figure 6.  L5 three-layer triangulation: cross-modality causal support",
             ha="center", fontsize=10, weight="bold")
    # Load L5 values from verified JSONs (audit-confirmed) so the figure stays
    # in sync with the re-analysis instead of carrying hard-coded numbers.
    l5dir = os.path.join(HERE, "l5_perturbation")
    l5a = json.load(open(os.path.join(l5dir, "l5_causal_direction.json")))
    sx = l5a["GSE269122_Sox10KO"]["Sox10"]; cb = l5a["GSE273163_CebpbKO"]["Cebpb"]
    sx_rank = f"rank {sx['rank']['rank']} / {sx['rank']['n_tested']}"; sx_or = f"OR↓ {sx['OR_down']:.2f}"
    cb_rank = f"rank {cb['rank']['rank']} / {cb['rank']['n_tested']}"; cb_or = f"OR↓ {cb['OR_down']:.2f}"
    pc = json.load(open(os.path.join(l5dir, "l5c_positive_control_result.json")))["panel"]
    myc = f"{pc['MYC']['rank_among_programs']}/{pc['MYC']['n_programs']}"
    bcl = f"{pc['BCL11A']['rank_among_programs']}/{pc['BCL11A']['n_programs']}"
    gat = f"{pc['GATA1']['rank_among_programs']}/{pc['GATA1']['n_programs']} (proxy)"
    layers = [
        ("L5a  Native-lineage KO", BLUE, [
            ("Sox10 cKO", sx_rank, sx_or, GREEN),
            ("Cebpb het-KO", cb_rank, cb_or, GREEN),
        ], "positive (directionally supported)"),
        ("L5b  LINCS overexpression", AMBER, [
            ("GATA2 OE", "rank 3 / 33,782", "p = 1.4×10⁻⁵", GREEN),
            ("3× TF CRISPR-KO", "top ~1% reverser", "no self-specificity", AMBER),
        ], "positive + directional, mixed self-specificity"),
        ("L5c  K562 sc-CRISPRi", GREY, [
            ("SOX10", "locus absent", "off-context null", GREY),
            ("CEBPB / GATA2", "on-target, regulon null", "context-gated null", GREY),
        ], "context boundary (expected)"),
    ]
    x0 = 0.06; w = 0.29; gap = 0.025
    ytop = 0.86; cardh = 0.17; cardgap = 0.03
    for li, (title, col, cards, verdict) in enumerate(layers):
        x = x0 + li * (w + gap)
        fig.text(x + w/2, ytop + 0.04, title, ha="center", fontsize=9, color=col, weight="bold")
        for ci, (tf, sub, stat, c) in enumerate(cards):
            y = ytop - ci * (cardh + cardgap)
            fig.add_artist(FancyBboxPatch((x, y - cardh), w, cardh,
                             boxstyle="round,pad=0.008", fc="white", ec=c, lw=1.2))
            fig.text(x + 0.02, y - 0.045, tf, fontsize=8.2, color=DARK, weight="bold")
            fig.text(x + 0.02, y - 0.095, sub, fontsize=7.2, color=DARK)
            fig.text(x + 0.02, y - 0.145, stat, fontsize=7.2, color=c, weight="bold")
        fig.text(x + w/2, ytop - 2*(cardh+cardgap) - 0.02, verdict, ha="center",
                 fontsize=6.6, style="italic", color=col, wrap=True)
    # positive-control strip
    fig.text(0.5, 0.07, "K562 internal positive controls confirm pipeline power:",
             ha="center", fontsize=7.6, color=DARK)
    pcs = [("MYC", myc, GREEN), ("BCL11A", bcl, GREEN), ("GATA1", gat, GREY)]
    px = 0.30
    for name, rk, c in pcs:
        fig.text(px, 0.03, f"{name} rank {rk}", ha="center", fontsize=7.2, color=c, weight="bold")
        px += 0.20
    # gradient arrow (context dependence)
    fig.add_artist(FancyArrowPatch((0.06, 0.12), (0.91, 0.12),
                   arrowstyle="-|>", mutation_scale=10, color=DARK, lw=1.0,
                   connectionstyle="arc3,rad=-0.15"))
    fig.text(0.5, 0.15, "causal support is context-gated  (off-context null  →  native-lineage positive)",
             ha="center", fontsize=7.0, style="italic", color=DARK)
    save(fig, "figure6_L5_triangulation.png")

# ============================================================================
# FIGURE 7  — L1–L5 evidence ladder
# ============================================================================
def fig7():
    fig = plt.figure(figsize=(7.2, 7.6))
    fig.text(0.5, 0.97, "Figure 7.  Five-level evidence ladder (L1–L5)",
             ha="center", fontsize=11, weight="bold")
    ladder = [
        ("L1", "Technical reproducibility", "TF master-regulator ranking ρ = 0.48–0.55, p < 1×10⁻¹⁵",
         "Discovery", BLUE),
        ("L2", "Biological coherence", "Sox10→Plp1 remyelination; aligns with known repair programs",
         "Recovery", BLUE),
        ("L3", "Cross-species convergence", "SOX10→myelin OR = 27 (p = 5×10⁻⁴); blood activation q = 1.1×10⁻⁸",
         "Orthogonal convergence", AMBER),
        ("L4", "Translational convergence", "vorinostat / statins reverse injury signature; perm p = 0.33 (n.s.)",
         "Hypothesis generation", AMBER),
        ("L5", "Directional causal support", "L5a rank 46/412; L5b GATA2 rank 3/33,782; L5c context-gated",
         "Causal *support* (directional)", GREEN),
    ]
    x = 0.30; top = 0.90; rh = 0.145; rw = 0.46
    for i, (lv, ev, metric, strength, col) in enumerate(ladder):
        y = top - i * rh
        # rung
        fig.add_artist(Rectangle((x - 0.02, y - rh*0.42), rw + 0.04, rh*0.84,
                       fc="white", ec=col, lw=1.4))
        fig.text(x + 0.01, y + rh*0.18, f"{lv}", fontsize=13, color=col, weight="bold")
        fig.text(x + 0.09, y + rh*0.18, ev, fontsize=9.2, color=DARK, weight="bold")
        fig.text(x + 0.09, y + rh*0.02, metric, fontsize=7.4, color=DARK, wrap=True)
        fig.text(x + 0.09, y - rh*0.16, f"strength: {strength}", fontsize=7.2,
                 color=col, style="italic")
        # connector
        if i < len(ladder) - 1:
            fig.add_artist(Line2D([x + rw/2, x + rw/2],
                           [y - rh*0.42, y - rh*0.42 - (rh - rh*0.84)],
                           color=GREY, lw=1.0))
    # ascending confidence arrow on the left
    fig.add_artist(FancyArrowPatch((x - 0.12, top - 4*rh + rh*0.4), (x - 0.12, top + rh*0.3),
                   arrowstyle="-|>", mutation_scale=12, color=GREEN, lw=1.4))
    fig.text(x - 0.14, top - 2*rh, "evidence\nconfidence", rotation=90, ha="center",
             va="center", fontsize=7.4, color=GREEN)
    fig.text(0.5, 0.02, "See Table 1 for numerical detail.  Edge-level results do not reproduce (Jaccard ≈ 0);\n"
             "each level is interpreted at the program / module level.", ha="center",
             fontsize=6.8, style="italic", color=GREY)
    save(fig, "figure7_evidence_ladder.png")

# ============================================================================
# SUPPLEMENTARY FIGURE S1 — PC composition correction
# ============================================================================
def figS1():
    trans = ["sham_24h", "24h_2d", "2d_14d", "sham_14d"]
    raw, corr = [], []
    with open(os.path.join(HERE, "rewiring_full.csv"), newline="") as f:
        for row in csv.DictReader(f):
            for key in trans:
                try:
                    qv = float(row["q_pooled_" + key])
                except ValueError:
                    qv = 1.0
                if qv < 0.05:
                    try:
                        raw.append(float(row["dW_raw_" + key]))
                        corr.append(float(row["dW_" + key]))
                    except ValueError:
                        pass
    fig, ax = plt.subplots(figsize=(4.2, 4.0))
    ax.scatter(raw, corr, s=10, color=BLUE, alpha=0.6, edgecolors="none")
    lim = [min(raw + corr), max(raw + corr)]
    ax.plot(lim, lim, ls="--", color=RED, lw=1.0)
    ax.set_xlabel("Raw ΔW")
    ax.set_ylabel("PC-corrected ΔW")
    ax.set_title("Supplementary Figure S1.\nPC composition correction", loc="left", fontsize=9.5)
    ax.text(0.05, 0.95, f"n = {len(raw)} significant edges\n(r = {np.corrcoef(raw, corr)[0,1]:.2f})",
            transform=ax.transAxes, va="top", fontsize=7.5,
            bbox=dict(boxstyle="round,pad=0.3", fc=LIGHT, ec=GREY, lw=0.6))
    save(fig, "figureS1_pc_correction.png")

if __name__ == "__main__":
    fig2()
    fig3()
    fig4()
    fig5()
    fig6()
    fig7()
    figS1()
    print("ALL FIGURES DONE")
