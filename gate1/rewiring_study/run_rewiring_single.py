"""run_rewiring_single.py — 参数化版 Step A rewiring（用于 cross-cohort 单数据集）。

v2 全部修复保留（n_perm=200/p_min≈0.005、BH-FDR+pooled FDR、PC 回归组成校正、
方向标签"关联增强/减弱"、manifest）。新增参数：
  --cache PATH           单数据集 cache npz
  --time-points "a,b,c" 时间点序列（相邻对推导 transitions）
  --species mouse/human  DoRothEA 图物种
  --out-prefix STR       输出文件前缀（避免覆盖整合 rewiring 产物）

输出：{prefix}rewiring_full.csv / {prefix}REWIRING_报告.md /
      {prefix}rewiring_top_edges.png / {prefix}run_manifest.txt
"""
import os
import sys
import time
import hashlib
import argparse
import numpy as np
import pandas as pd

sys.path.insert(0, "C:/D 盘/科研/虚拟敲除/gate1")
sys.path.insert(0, "C:/D 盘/科研/虚拟敲除/gate1/tsc_gnn")
from grn import build_dorothea_grn
from rewiring import rewiring_table

OUT_DIR = "C:/D 盘/科研/虚拟敲除/gate1/rewiring_study"

TF_NOTE = {
    "SPI1": "PU.1 髓系/小胶质", "SOX10": "少突胶质/髓鞘", "CEBPE": "粒细胞/炎症",
    "STAT4": "IL-12/IFNγ", "PAX5": "B细胞", "SOX2": "神经干/重编程", "RELA": "NF-kB",
    "NFKB1": "NF-kB", "IRF8": "小胶质", "CEBPB": "炎症/急性相", "EGR1": "即刻早期",
    "RUNX3": "CD8/T", "GATA2": "造血", "ERG": "内皮", "MAF": "免疫", "NFE2": "红系",
    "SPI1B": "髓系", "SNAI2": "EMT", "NR2F2": "血管", "ERGR": "内皮",
    "CEBPA": "炎症/髓系", "IRF1": "I 型干扰素", "STAT1": "IFN", "STAT3": "IL-6/急性相",
    "PPARG": "脂肪/巨噬", "MEF2C": "心脏/神经", "TCF7L2": "Wnt", "FOXP3": "Treg",
}


def stratified_subsample(time_label, n_per_time, seed=0):
    rng = np.random.default_rng(seed)
    tl = np.asarray(time_label)
    idx = []
    for t in np.unique(tl):
        m = np.where(tl == t)[0]
        if len(m) > n_per_time:
            m = rng.choice(m, size=n_per_time, replace=False)
        idx.append(m)
    return np.concatenate(idx)


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8192), b""):
            h.update(block)
    return h.hexdigest()


def main(n_per_time=2000, n_perm=200, n_pc=10, aff_q=0.5,
         cache_path=None, time_points_str="sham,24h,2d,14d",
         species="mouse", out_prefix=""):
    t0 = time.time()
    print(f"== Step A (single-cohort) rewiring: {cache_path} ==")

    with np.load(cache_path, allow_pickle=True) as z:
        X = z["X"]
        genes = z["genes"]
        state = z["state"]
        tl = z["time_label"]
    print(f"  缓存: X={X.shape} genes={genes.shape}")

    sub = stratified_subsample(tl, n_per_time, seed=0)
    Xs = np.ascontiguousarray(X[sub], dtype=np.float32)
    tls = tl[sub]
    states = state[sub]
    for t in sorted(set(tls)):
        print(f"    {t}: {int((tls == t).sum())} 细胞")
    print(f"  抽样后: {Xs.shape[0]} 细胞, {Xs.nbytes / 1e6:.0f}MB")

    time_points = [t.strip() for t in time_points_str.split(",")]
    transitions = [(time_points[i], time_points[i + 1])
                   for i in range(len(time_points) - 1)]
    print(f"  时间点: {time_points}")
    print(f"  转移: {transitions}")

    print(f"  [grn] 构建 {species} DoRothEA(ABC) 有向图 + state-affinity ...")
    grn = build_dorothea_grn(
        list(genes), species=species, confidence_levels=("A", "B", "C"),
        X=Xs, state=states, n_sub=4000, seed=1)
    print(f"       边={grn['n_edges']} 有向={grn['is_directed']}")

    print(f"  [rewire] n_perm={n_perm}, n_pc={n_pc}, aff_q={aff_q}")
    df, transitions, time_points, manifest = rewiring_table(
        grn, list(genes), Xs, tls, state=states,
        transitions=transitions, n_perm=n_perm, seed=2,
        a_aff_threshold=aff_q, n_pc=n_pc)

    OUT_CSV = os.path.join(OUT_DIR, f"{out_prefix}rewiring_full.csv")
    OUT_MD = os.path.join(OUT_DIR, f"{out_prefix}REWIRING_报告.md")
    OUT_PNG = os.path.join(OUT_DIR, f"{out_prefix}rewiring_top_edges.png")
    OUT_MANIFEST = os.path.join(OUT_DIR, f"{out_prefix}run_manifest.txt")

    df.to_csv(OUT_CSV, index=False)
    print(f"  [csv] -> {OUT_CSV} ({df.shape})")

    manifest["command"] = (f"python run_rewiring_single.py --cache {cache_path} "
                            f"--time-points {time_points_str} --species {species} "
                            f"--out-prefix {out_prefix} "
                            f"--n-per-time {n_per_time} --n-perm {n_perm} "
                            f"--n-pc {n_pc} --aff-q {aff_q}")
    manifest["cache_sha256"] = _sha256(cache_path)
    manifest["cache_path"] = cache_path
    manifest["subsample_seed"] = 0
    manifest["grn_seed"] = 1
    manifest["n_per_time"] = n_per_time
    manifest["elapsed_sec"] = round(time.time() - t0, 1)
    with open(OUT_MANIFEST, "w") as f:
        f.write("== Rewiring Run Manifest (single-cohort) ==\n")
        f.write(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        for k, v in manifest.items():
            f.write(f"{k}: {v}\n")
    print(f"  [manifest] -> {OUT_MANIFEST}")

    _write_report(df, transitions, time_points, manifest, t0, n_pc,
                  OUT_MD, out_prefix)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        _viz(df, transitions, OUT_PNG, n_pc)
        print(f"  [png] -> {OUT_PNG}")
    except Exception as e:
        print(f"  [viz] 跳过: {e}")

    print(f"\n== 完成 ({time.time() - t0:.0f}s) ==")


def _write_report(df, transitions, time_points, manifest, t0, n_pc, OUT_MD, prefix):
    L = []
    L.append(f"# {prefix or ''}单队列 Rewiring 报告（cross-cohort 验证）")
    L.append(f"（{time.strftime('%Y-%m-%d')}）")
    L.append("")
    L.append(f"- 缓存: {manifest['cache_path']}")
    L.append(f"- 抽样细胞: {manifest['n_cells']}（每时间点 ≤{manifest.get('n_per_time','?')}）")
    L.append(f"- 测试边: {manifest['n_edges_tested']}（state-conditioned, |A_aff| 前 50%）")
    L.append(f"- 置换次数: n_perm={manifest['n_perm']} → p_min={manifest['p_min']:.4f}")
    L.append(f"- PC 回归: n_pc={n_pc}（去除细胞类型组成效应）")
    L.append(f"- 时间点: {time_points}；转移: {transitions}")
    L.append("")
    L.append("## 显著性汇总")
    L.append("")
    L.append("| 转移 | p<0.05 | BH-FDR<0.05 | BH-FDR<0.1 | "
             "Pooled q<0.05 | Pooled q<0.1 | Pooled q<0.2 |")
    L.append("|---|---|---|---|---|---|---|")
    for (t1, t2) in transitions:
        p = df[f"p_{t1}_{t2}"].values
        fdr = df[f"fdr_bh_{t1}_{t2}"].values
        q = df[f"q_pooled_{t1}_{t2}"].values
        L.append(
            f"| {t1}→{t2} | {int((p < 0.05).sum())} | "
            f"{int((fdr < 0.05).sum())} | {int((fdr < 0.1).sum())} | "
            f"{int((q < 0.05).sum())} | {int((q < 0.1).sum())} | "
            f"{int((q < 0.2).sum())} |")
    L.append("")
    L.append("## Top 15 重连边（按 |ΔW|, pooled q 标注）")
    L.append("")
    for (t1, t2) in transitions:
        L.append(f"### {t1} → {t2}")
        L.append("")
        sub = df.reindex(
            df[f"dW_{t1}_{t2}"].abs().sort_values(ascending=False).index).head(15)
        L.append("| TF | target | TF功能 | DoRothEA先验 | 方向 | "
                 "ΔW(PC) | ΔW(raw) | p | pooled q | "
                 f"r_{t1}(PC) | r_{t2}(PC) |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|")
        for _, r in sub.iterrows():
            tf = r["tf"]
            note = TF_NOTE.get(tf.upper(), "")
            d = r[f"dW_{t1}_{t2}"]
            direction = "关联增强" if d > 0 else "关联减弱"
            q_val = r[f"q_pooled_{t1}_{t2}"]
            q_mark = " **" if q_val < 0.1 else ""
            L.append(
                f"| {tf} | {r['target']} | {note} | {r['dorothea_sign']} | "
                f"{direction} | {d:+.2f} | {r[f'dW_raw_{t1}_{t2}']:+.2f} | "
                f"{r[f'p_{t1}_{t2}']:.1e} | {q_val:.2e}{q_mark} | "
                f"{r[f'r_{t1}']:+.2f} | {r[f'r_{t2}']:+.2f} |")
        L.append("")
    with open(OUT_MD, "w") as f:
        f.write("\n".join(L))


def _viz(df, transitions, out_png, n_pc):
    import matplotlib.pyplot as plt
    n = len(transitions)
    cols = min(2, n)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(14, 6 * rows))
    if n == 1:
        axes = [axes]
    else:
        axes = axes.ravel()
    for i, (t1, t2) in enumerate(transitions):
        ax = axes[i]
        sub = df.reindex(
            df[f"dW_{t1}_{t2}"].abs().sort_values(ascending=False).index).head(10)
        labels = [f"{r.tf}->{r.target}" for _, r in sub.iterrows()]
        vals = sub[f"dW_{t1}_{t2}"].values
        colors = ["#2980b9" if v > 0 else "#c0392b" for v in vals]
        ax.barh(range(len(vals))[::-1], vals, color=colors)
        ax.set_yticks(range(len(labels))[::-1])
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel(f"dW ({t1}->{t2})")
        ax.set_title(f"{t1}->{t2} (PC-corrected)")
        ax.axvline(0, color="gray", linewidth=0.5)
    for j in range(n, len(axes)):
        axes[j].axis("off")
    fig.suptitle(f"Top rewired edges (n_pc={n_pc}, pooled FDR)", fontsize=12)
    plt.tight_layout()
    plt.savefig(out_png, dpi=120)
    plt.close()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True)
    ap.add_argument("--time-points", default="sham,24h,2d,14d")
    ap.add_argument("--species", default="mouse")
    ap.add_argument("--out-prefix", default="")
    ap.add_argument("--n-per-time", type=int, default=2000)
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--n-pc", type=int, default=10)
    ap.add_argument("--aff-q", type=float, default=0.5)
    args = ap.parse_args()
    main(n_per_time=args.n_per_time, n_perm=args.n_perm,
         n_pc=args.n_pc, aff_q=args.aff_q,
         cache_path=args.cache, time_points_str=args.time_points,
         species=args.species, out_prefix=args.out_prefix)
