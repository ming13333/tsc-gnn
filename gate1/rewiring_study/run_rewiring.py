"""run_rewiring.py v2 — Step A: 卒中时间×状态 DoRothEA rewiring（5 条阻断级已修）。

v2 修复：
1. n_perm=200 (p_min≈0.005) — 前台 <3min 完成（预计算 sgemv 优化）
2. BH-FDR + pooled FDR 双重多重校正
3. PC 回归组成校正 (n_pc=10) — 去除细胞类型组成混淆
4. 方向标签改为"关联增强/关联减弱"
5. manifest.txt 落盘（命令/版本/seed/参数）

输出：
  rewiring_full.csv              所有边的 r/dW/p/FDR/q（PC-corrected + raw）
  REWIRING_报告_v3_2026-07-09.md  解读 + Top rewired edges（pooled FDR）
  rewiring_top_edges.png          Top rewired 网络图
  run_manifest.txt                可复现清单
"""
import os
import sys
import time
import hashlib
import json
import numpy as np
import pandas as pd

sys.path.insert(0, "C:/D 盘/科研/虚拟敲除/gate1")
sys.path.insert(0, "C:/D 盘/科研/虚拟敲除/gate1/tsc_gnn")
from grn import build_dorothea_grn
from rewiring import rewiring_table

CACHE = "C:/D 盘/科研/虚拟敲除/gate1/rewiring_study/rewiring_doro_cache.npz"
OUT_DIR = "C:/D 盘/科研/虚拟敲除/gate1/rewiring_study"
OUT_CSV = os.path.join(OUT_DIR, "rewiring_full.csv")
OUT_MD = os.path.join(OUT_DIR, "REWIRING_报告_v3_2026-07-09.md")
OUT_PNG = os.path.join(OUT_DIR, "rewiring_top_edges.png")
OUT_MANIFEST = os.path.join(OUT_DIR, "run_manifest.txt")

TRANSITIONS = [("sham", "24h"), ("24h", "2d"), ("2d", "14d"), ("sham", "14d")]

# TF 功能注释（鼠/人通用）
TF_NOTE = {
    "SPI1": "PU.1 髓系/小胶质", "SOX10": "少突胶质/髓鞘", "CEBPE": "粒细胞/炎症",
    "STAT4": "IL-12/IFNγ", "PAX5": "B细胞", "SOX2": "神经干/重编程", "RELA": "NF-kB",
    "NFKB1": "NF-kB", "IRF8": "小胶质", "CEBPB": "炎症/急性相", "EGR1": "即刻早期",
    "RUNX3": "CD8/T", "GATA2": "造血", "ERG": "内皮", "MAF": "免疫", "NFE2": "红系",
    "SPI1B": "髓系", "SNAI2": "EMT", "NR2F2": "血管", "ERGR": "内皮", "GATA2R": "造血",
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
    """计算文件 SHA256（用于 manifest）。"""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8192), b""):
            h.update(block)
    return h.hexdigest()


def _fmt_cmd(n_per_time, n_perm, n_pc, aff_q):
    return (f"python run_rewiring.py "
            f"--n-per-time {n_per_time} --n-perm {n_perm} "
            f"--n-pc {n_pc} --aff-q {aff_q}")


def main(n_per_time=2000, n_perm=200, n_pc=10, aff_q=0.5):
    t0 = time.time()
    print("== Step A v2: 卒中时间×状态 DoRothEA rewiring (5-fix) ==")

    # ── 加载缓存 ──
    with np.load(CACHE, allow_pickle=True) as z:
        X = z["X"]
        genes = z["genes"]
        state = z["state"]
        tl = z["time_label"]
    print(f"  缓存: X={X.shape} genes={genes.shape}")

    # ── 分层抽样（seed=0）──
    sub = stratified_subsample(tl, n_per_time, seed=0)
    Xs = np.ascontiguousarray(X[sub], dtype=np.float32)
    tls = tl[sub]
    states = state[sub]
    for t in sorted(set(tls)):
        print(f"    {t}: {int((tls == t).sum())} 细胞")
    print(f"  抽样后: {Xs.shape[0]} 细胞, {Xs.nbytes / 1e6:.0f}MB")

    # ── DoRothEA 因果图 ──
    print("  [grn] 构建 mouse DoRothEA(ABC) 有向图 + state-affinity ...")
    grn = build_dorothea_grn(
        list(genes), species="mouse", confidence_levels=("A", "B", "C"),
        X=Xs, state=states, n_sub=4000, seed=1)  # grn seed=1（不同于 subsample）
    print(f"       边={grn['n_edges']} 有向={grn['is_directed']}")

    # ── rewiring 分析 ──
    print(f"  [rewire] n_perm={n_perm}, n_pc={n_pc}, aff_q={aff_q}")
    df, transitions, time_points, manifest = rewiring_table(
        grn, list(genes), Xs, tls, state=states,
        transitions=TRANSITIONS, n_perm=n_perm, seed=2,  # perm seed=2
        a_aff_threshold=aff_q, n_pc=n_pc)

    # ── 保存 CSV ──
    df.to_csv(OUT_CSV, index=False)
    print(f"  [csv] -> {OUT_CSV} ({df.shape})")

    # ── Manifest ──
    manifest["command"] = _fmt_cmd(n_per_time, n_perm, n_pc, aff_q)
    manifest["cache_sha256"] = _sha256(CACHE)
    manifest["cache_path"] = CACHE
    manifest["subsample_seed"] = 0
    manifest["grn_seed"] = 1
    manifest["n_per_time"] = n_per_time
    manifest["elapsed_sec"] = round(time.time() - t0, 1)
    with open(OUT_MANIFEST, "w") as f:
        f.write("== Rewiring Run Manifest ==\n")
        f.write(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        for k, v in manifest.items():
            f.write(f"{k}: {v}\n")
    print(f"  [manifest] -> {OUT_MANIFEST}")

    # ── 报告 ──
    _write_report(df, transitions, time_points, manifest, t0, n_pc)
    print(f"  [md] -> {OUT_MD}")

    # ── 可视化 ──
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        _viz(df, transitions, OUT_PNG, n_pc)
        print(f"  [png] -> {OUT_PNG}")
    except Exception as e:
        print(f"  [viz] 跳过: {e}")

    print(f"\n== 完成 ({time.time() - t0:.0f}s) ==")


def _write_report(df, transitions, time_points, manifest, t0, n_pc):
    L = []
    L.append("# Step A v2 — 卒中时间×状态 DoRothEA 重布线报告")
    L.append(f"（Bug Audit 5 条阻断级已修 · {time.strftime('%Y-%m-%d')}）")
    L.append("")
    L.append("## 参数")
    L.append("")
    L.append(f"- 抽样细胞: {manifest['n_cells']}（每时间点 ≤{manifest.get('n_per_time','?')}）")
    L.append(f"- 测试边: {manifest['n_edges_tested']}（state-conditioned, |A_aff| 前 50%）")
    L.append(f"- 置换次数: n_perm={manifest['n_perm']} → p_min={manifest['p_min']:.4f}")
    L.append(f"- PC 回归: n_pc={n_pc}（去除细胞类型组成效应）")
    L.append(f"- 种子: subsample=0, grn=1, perm=2")
    L.append("")

    # ── 显著性汇总 ──
    L.append("## 显著性汇总（多重校正对比）")
    L.append("")
    L.append("| 转移 | p<0.05 (未校正) | BH-FDR<0.05 | BH-FDR<0.1 | "
             "Pooled q<0.05 | Pooled q<0.1 | Pooled q<0.2 |")
    L.append("|---|---|---|---|---|---|---|")
    for (t1, t2) in transitions:
        p = df[f"p_{t1}_{t2}"].values
        fdr = df[f"fdr_bh_{t1}_{t2}"].values
        q = df[f"q_pooled_{t1}_{t2}"].values
        dw = df[f"dW_{t1}_{t2}"].values
        L.append(
            f"| {t1}→{t2} | {int((p < 0.05).sum())} | "
            f"{int((fdr < 0.05).sum())} | {int((fdr < 0.1).sum())} | "
            f"{int((q < 0.05).sum())} | {int((q < 0.1).sum())} | "
            f"{int((q < 0.2).sum())} |")
    L.append("")
    L.append("- **BH-FDR**：逐边 p 值做 Benjamini-Hochberg 校正。"
             "因 n_perm=200 使 p_min≈0.005，在 "
             f"{manifest['n_edges_tested']} 边上 BH-FDR<0.05 需 p<"
             f"{0.05/manifest['n_edges_tested']:.1e}（不可达），故 BH-FDR 可能过保守。")
    L.append("- **Pooled q**：置换 pooled FDR（标准化 + 池化 null），"
             "不依赖 p 分辨率，更适合 n_perm 有限的场景。")
    L.append("")

    # ── 组成校正对比 ──
    L.append("## PC 校正效果（组成混淆控制）")
    L.append("")
    L.append("| 转移 | |ΔW| 中位数 (raw) | |ΔW| 中位数 (PC-corrected) | "
             "Pearson(raw, corrected) |")
    L.append("|---|---|---|---|")
    for (t1, t2) in transitions:
        dw_raw = df[f"dW_raw_{t1}_{t2}"].abs().median()
        dw_pc = df[f"dW_{t1}_{t2}"].abs().median()
        corr = np.corrcoef(df[f"dW_raw_{t1}_{t2}"],
                           df[f"dW_{t1}_{t2}"])[0, 1]
        L.append(f"| {t1}→{t2} | {dw_raw:.3f} | {dw_pc:.3f} | {corr:.2f} |")
    L.append("")
    L.append(f"- PC 回归（n_pc={n_pc}）残差化后，ΔW 量级与方向的变化反映"
             "组成校正的影响。若 raw→corrected 相关高（>0.8），"
             "说明组成非主要驱动；若低，说明组成混淆严重。")
    L.append("")

    # ── Top rewired edges per transition ──
    for (t1, t2) in transitions:
        L.append(f"## {t1} → {t2}：Top 15 重连边（按 |ΔW|, pooled q 标注）")
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
        L.append(f"_方向标签说明：'关联增强'= TF-target 共表达耦合上升（ΔW>0）；"
                 f"'关联减弱'= 耦合下降（ΔW<0）。DoRothEA 先验方向"
                 f"（activation/repression）独立报告，不与耦合变化混为一谈。_")
        L.append("")

    # ── 生物学解读 ──
    L.append("## 生物学解读")
    L.append("")
    L.append("- 方法定位：本分析量化**急性时间轴上 TF→target 调控边的"
             "共表达耦合重连**（edge rewiring），属可解释性主轴。")
    L.append(f"- PC 回归（n_pc={n_pc}）控制了细胞类型组成混淆："
             "残差化后的 ΔW 反映去除主要组成变异后的耦合变化，"
             "更接近'调控重连'而非'组成偏移'。")
    L.append("- DoRothEA 提供**有向因果方向**（TF→target），"
             "rewiring 可解释为'谁的调控耦合在卒中时间轴上被增强/减弱'，"
             "线性模型无法给出。")
    L.append("- pooled q 值不依赖 p 分辨率，"
             f"在 n_perm={manifest['n_perm']} 下仍可做多重校正。")
    L.append("")
    L.append(f"耗时 {time.time() - t0:.0f}s")

    with open(OUT_MD, "w") as f:
        f.write("\n".join(L))


def _viz(df, transitions, out_png, n_pc):
    """Top rewired edges 条形图（4 个转移 × Top 10）。"""
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
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
    fig.suptitle(
        f"Top rewired edges per transition (n_pc={n_pc}, pooled FDR)", fontsize=12)
    plt.tight_layout()
    plt.savefig(out_png, dpi=120)
    plt.close()


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-per-time", type=int, default=2000)
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--n-pc", type=int, default=10)
    ap.add_argument("--aff-q", type=float, default=0.5)
    args = ap.parse_args()
    main(n_per_time=args.n_per_time, n_perm=args.n_perm,
         n_pc=args.n_pc, aff_q=args.aff_q)
