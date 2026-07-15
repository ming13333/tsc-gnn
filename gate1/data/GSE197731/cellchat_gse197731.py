#!/usr/bin/env python
# cellchat_gse197731.py
# ---------------------------------------------------------------------------
# 把下载好的 GSE197731 (第二个 24h 队列) 接入 CellChat 管线, 并做 24h 跨队列一致性。
# 复用 cellchat_py 的注释/重布线函数, 产出与 cellchat_rewiring.csv 同 schema 的结果。
#
# 数据设计 (实测确认, 无 sham):
#   8 样本 = 24h/48h × WT/Prdx1-KO × 同侧(Ipsil,缺血)/对侧(Cont,非缺血)
# 对照框架: Cont(对侧) -> "sham", Ipsil(同侧) -> "MCAO"  (神经科学常规半球内对照)
#
# 处理策略:
#   - 加载全部 8 样本, 按基因型拆成 GSE197731_WT / GSE197731_KO 两个 study 分别跑
#     (避免 WT/KO 混合到同一 transition 造成基因型混淆)
#   - time: 24h / 48h  -> 各产生 transition "24h"(Cont vs Ipsil) 与 "48h"
#   - 核心跨队列 24h 一致性: GSE197731_WT(24h) vs GSE174574(24h, 同为 WT) -> Jaccard
#   - 次级: Prdx1-KO 基因型效应 (GSE197731_KO 24h sig 与 WT 24h sig 的重叠/差异)
#
# 用法: python cellchat_gse197731.py
# 依赖: cellchat_py.py, cellchat_lrdb.py (位于 gate1/rewiring_study)
# ---------------------------------------------------------------------------
import os, sys, glob, re
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REW = os.path.abspath(os.path.join(HERE, "..", "..", "rewiring_study"))
sys.path.insert(0, REW)
import cellchat_py as cc
from cellchat_lrdb import LR_DB, pathway_of

OUT = HERE


def _parse_filename(fname):
    """GSM5929213_h24_wt_con_counts.csv.gz -> {time, genotype, side}。"""
    m = re.search(r"_h(\d+)_(wt|ko)_(con|ip)_", fname, re.I)
    if not m:
        return None
    h, geno, side = m.groups()
    return {"time": f"{h}h", "genotype": geno.upper(), "side": side.lower()}


def _load_sample_counts(path):
    if path.endswith(".csv.gz"):
        df = pd.read_csv(path, sep=",", index_col=0, compression="gzip")
        # GSE197731 格式: 行=基因(索引), 列=细胞(barcode) -> 矩阵为 genes×cells
        # 转置成 cells×genes 供 CellChat 使用; 基因符号取自行索引
        X = df.values.T.astype(np.float32)          # cells × genes
        genes = [str(g) for g in df.index]
    else:
        from scipy.io import mmread
        M = mmread(path).T.tocsr().astype(np.float32)
        gf = glob.glob(os.path.join(os.path.dirname(path), "*genes.tsv.gz")) or \
             glob.glob(os.path.join(os.path.dirname(path), "*features.tsv.gz"))
        genes = pd.read_csv(gf[0], sep="\t", header=None, compression="gzip")[1].astype(str).tolist() if gf \
            else [f"g{i}" for i in range(M.shape[1])]
        X = (M.toarray() if M.shape[0] < 60000 else M).astype(np.float32)
        if X.shape[0] < X.shape[1]:                 # 10x mtx 常为 genes×cells, 转置为 cells×genes
            X = X.T
    if X.max() > 100:                               # 原始 counts -> log1p 归一化
        rs = X.sum(axis=1); rs[rs == 0] = 1
        X = np.log1p(X / rs[:, None] * 1e4)
    return X, genes


def _annotate(samples):
    for s in samples:
        s["cell_type"] = cc._annotate_markers(s["X"], s["genes"], cc.CELLTYPE_MARKERS)
        s["X"], s["cell_type"] = cc._subsample_ct(s["X"], s["cell_type"])
    cc._align_genes_per_study(samples)


def load_all():
    """扫描 HERE 下所有 *_counts.csv.gz, 解析文件名, Cont->sham/Ipsil->MCAO 映射。"""
    out = []
    for path in sorted(glob.glob(os.path.join(HERE, "*", "*_counts.csv.gz"))):
        fname = os.path.basename(path)
        meta = _parse_filename(fname)
        if meta is None:
            continue
        X, genes = _load_sample_counts(path)
        condition = "sham" if meta["side"] == "con" else "MCAO"
        rec = {"X": X, "genes": genes, "condition": condition, "time": meta["time"],
               "genotype": meta["genotype"], "study": "GSE197731", "src": fname}
        out.append(rec)
        print(f"  [load] {fname}: cells={X.shape[0]} {condition}/{meta['time']}/{meta['genotype']}")
    return out


def _run_study(samples, study_label, pairs):
    for s in samples:
        s["study"] = study_label
    cc._align_genes_per_study(samples)
    rows = []
    for lig, rec in pairs:
        rows.extend(cc.rewiring_for_study(samples, study_label, lig, rec))
    return pd.DataFrame(rows)


def _add_sig(df):
    if df.empty:
        return df
    df["pathway"] = [pathway_of(l, r) for l, r in zip(df.ligand, df.receptor)]
    df["sig"] = (np.minimum(df.pLig, df.pRec) < 0.05) & (df.log2FC.abs() >= 0.5)
    return df


def main():
    pairs = sorted({(l, r) for l, r, _ in LR_DB})

    print("=== 加载 GSE197731 全部 8 样本 ===")
    all_s = load_all()
    if not all_s:
        raise SystemExit("无可用样本; 请确认 *_counts.csv.gz 已下载到 GSM*/ 目录")

    wt = [s for s in all_s if s["genotype"] == "WT"]
    ko = [s for s in all_s if s["genotype"] == "KO"]
    print(f"  WT 样本={len(wt)}  KO 样本={len(ko)}")

    # ---- WT study ----
    print("\n=== 运行 GSE197731_WT CellChat ===")
    _annotate(wt)
    _annotate(ko)                      # KO 同样需注释 cell_type
    dwt = _run_study(wt, "GSE197731_WT", pairs)
    dwt = _add_sig(dwt)
    # ---- KO study ----
    print("=== 运行 GSE197731_KO CellChat ===")
    dko = _run_study(ko, "GSE197731_KO", pairs) if ko else pd.DataFrame()
    dko = _add_sig(dko)

    full = pd.concat([dwt, dko], ignore_index=True) if not dko.empty else dwt
    full.to_csv(os.path.join(OUT, "cellchat_gse197731_all.csv"), index=False)
    print(f"\n  GSE197731 全样本 LR-转变记录={len(full)}, 显著={int(full.sig.sum())}")
    for lab, d in [("GSE197731_WT", dwt), ("GSE197731_KO", dko)]:
        if d.empty:
            continue
        for t in sorted(d.transition.unique()):
            sub = d[d.transition == t]
            print(f"    {lab} {t}: 候选={len(sub)} 显著={int(sub.sig.sum())}")

    # ---- 核心: 24h 跨队列一致性 (WT only, 均为缺血vs非缺血) ----
    sig197 = set((r.ligand, r.receptor) for _, r in dwt[dwt.sig & (dwt.transition == "24h")].iterrows())
    s174 = cc.load_gse174()
    # 与 GSE197731 一致: 每 CT 封顶, 既控内存又让两队列 MWU 规模可比 (公平跨队列对比)
    for s in s174:
        s["X"], s["cell_type"] = cc._subsample_ct(s["X"], s["cell_type"])
    cc._align_genes_per_study(s174)
    d174 = _add_sig(_run_study(s174, "GSE174574", pairs))
    sig174 = set((r.ligand, r.receptor) for _, r in d174[d174.sig & (d174.transition == "24h")].iterrows())
    inter = sig197 & sig174
    union = sig197 | sig174
    jac = len(inter) / max(1, len(union))
    print(f"\n=== 24h 跨队列一致性 (GSE174574 vs GSE197731, 均为 WT, 缺血vs非缺血) ===")
    print(f"  GSE174574 24h 显著LR对: {len(sig174)}")
    print(f"  GSE197731 24h 显著LR对: {len(sig197)}")
    print(f"  交集={len(inter)}  Jaccard={jac:.2f}")
    for p in sorted(inter):
        print(f"    ✔ {p[0]}->{p[1]}")

    # ---- 通路级 / 细胞型级 跨队列一致性 ----
    # LR对Jaccard受 GSE174574 24h 低功率(仅少数显著对)限制, 补更有力的通路/细胞型视角
    dwt24 = dwt[dwt.transition == "24h"]
    lr2path = {(r.ligand, r.receptor): r.pathway for _, r in d174.iterrows()}
    lr2path.update({(r.ligand, r.receptor): r.pathway for _, r in dwt24.iterrows()})
    lr2ct = {(r.ligand, r.receptor): (r.sender, r.receiver) for _, r in d174.iterrows()}
    lr2ct.update({(r.ligand, r.receptor): (r.sender, r.receiver) for _, r in dwt24.iterrows()})
    p174 = {lr2path.get(p) for p in sig174 if lr2path.get(p)}
    p197 = {lr2path.get(p) for p in sig197 if lr2path.get(p)}
    ct174 = {lr2ct.get(p) for p in sig174 if lr2ct.get(p)}
    ct197 = {lr2ct.get(p) for p in sig197 if lr2ct.get(p)}
    pjac = len(p174 & p197) / max(1, len(p174 | p197))
    print(f"\n  通路级 24h 一致性: GSE174574={sorted(p174)}")
    print(f"    GSE197731_WT={sorted(p197)}")
    print(f"    共享通路={sorted(p174 & p197)}  (Jaccard={pjac:.2f})")
    print(f"  显著LR对的发送→接收细胞型(24h): GSE174574={len(ct174)}对 GSE197731_WT={len(ct197)}对 共享={len(ct174 & ct197)}")

    # ---- 次级: Prdx1-KO 基因型效应 (24h) ----
    ko24 = set((r.ligand, r.receptor) for _, r in dko[dko.sig & (dko.transition == "24h")].iterrows()) if not dko.empty else set()
    if ko24:
        inter_g = sig197 & ko24
        print(f"\n=== Prdx1-KO 基因型效应 (24h: WT vs KO, 同 Ipsil/Cont 框架) ===")
        print(f"  WT 24h 显著={len(sig197)}  KO 24h 显著={len(ko24)}  共享={len(inter_g)}")
        print(f"  仅 WT 特有={len(sig197 - ko24)}  仅 KO 特有={len(ko24 - sig197)}")

    with open(os.path.join(OUT, "cross_cohort_24h_consistency.txt"), "w") as f:
        f.write(f"GSE174574_24h_sig={sorted(sig174)}\n")
        f.write(f"GSE197731_WT_24h_sig={sorted(sig197)}\n")
        f.write(f"intersection={sorted(inter)}\n")
        f.write(f"jaccard={jac:.3f}\n")
        f.write(f"GSE174574_24h_pathways={sorted(p174)}\n")
        f.write(f"GSE197731_24h_pathways={sorted(p197)}\n")
        f.write(f"pathway_jaccard={pjac:.3f}\n")
        f.write(f"GSE174574_24h_ctpairs={sorted(ct174)}\n")
        f.write(f"GSE197731_24h_ctpairs={sorted(ct197)}\n")
        f.write(f"ctpair_shared={sorted(ct174 & ct197)}\n")
        f.write(f"GSE197731_KO_24h_sig={sorted(ko24)}\n")
        if ko24:
            f.write(f"ko_wt_shared={sorted(inter_g)}\n")
    print("\n=== 完成 -> cellchat_gse197731_all.csv + cross_cohort_24h_consistency.txt ===")


if __name__ == "__main__":
    main()
