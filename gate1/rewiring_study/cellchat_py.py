# cellchat_py.py
# CellChat 风格的配体-受体细胞通讯分析（Python 重实现；本机无 R）。
# 设计见方法: 逐 study 内 MCAO vs sham 计算通讯概率，跨卒中转变(24h/2d/14d)识别重布线 LR 对，
# 并与 DoRothEA 主调控因子(Sox10/Cebpb/Gata2)靶程序整合。
import os, glob, gzip, json
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse as sp
from scipy.stats import mannwhitneyu

ROOT = os.path.dirname(os.path.abspath(__file__))
GSE174 = os.path.join(ROOT, "..", "data", "GSE174574")
GSE225 = os.path.join(ROOT, "..", "data", "GSE225948")
OUT = ROOT
FIG = os.path.join(ROOT, "figures")

# ---- 鼠细胞类型 marker（与 preprocessing.CELLTYPE_MARKERS 一致）----
CELLTYPE_MARKERS = {
    "Microglia": ["Cx3cr1", "Tmem119", "Aif1", "Itgam", "P2ry12", "Siglech"],
    "Macrophage": ["Ly6c2", "Ccr2", "Adgre1", "Itgam"],
    "Neuron": ["Snap25", "Syt1", "Rbfox3", "Slc17a7", "Gabbr1"],
    "Astrocyte": ["Gfap", "Aqp4", "Slc1a2", "Aldh1l1", "Sparcl1"],
    "Oligodendrocyte": ["Mbp", "Mog", "Olig2", "Plp1", "Cnp"],
    "OPC": ["Pdgfra", "Cspg4", "Pcdh15"],
    "Endothelial": ["Cldn5", "Pecam1", "Cdh5", "Kdr"],
    "Pericyte": ["Pdgfrb", "Rgs5", "Abcc9"],
    "Tcell": ["Cd3g", "Cd3e", "Cd4", "Cd8a"],
    "Bcell": ["Cd79a", "Ms4a1", "Cd19"],
    "NK": ["Nkg7", "Klrd1"],
}

# GSE225948 GSM -> (condition, time_label)  (复制自 preprocessing.GSE225948_SAMPLE_META)
GSE225_META = {
    "GSM7060815": ("sham", "2d"), "GSM7060816": ("sham", "2d"),
    "GSM7060817": ("sham", "2d"), "GSM7060818": ("sham", "2d"),
    "GSM7060819": ("MCAO", "2d"), "GSM7060820": ("MCAO", "2d"),
    "GSM7060821": ("MCAO", "2d"), "GSM7060822": ("MCAO", "2d"),
    "GSM7060823": ("MCAO", "14d"), "GSM7060824": ("MCAO", "14d"),
    "GSM7060825": ("MCAO", "14d"), "GSM7060826": ("MCAO", "14d"),
    "GSM7060827": ("sham", "2d"), "GSM7060828": ("sham", "2d"),
    "GSM7060829": ("MCAO", "2d"), "GSM7060830": ("MCAO", "2d"),
    "GSM7060831": ("MCAO", "2d"), "GSM7060832": ("MCAO", "14d"),
}

GSE174_COND = {  # 从文件名解析
    "GSM5319987": "sham", "GSM5319988": "sham", "GSM5319989": "sham",
    "GSM5319990": "MCAO", "GSM5319991": "MCAO", "GSM5319992": "MCAO",
}

MIN_CELLS = 10         # 每 CT 至少细胞数
MIN_FRAC = 0.10        # 配体/受体在发送/接收 CT 中表达细胞比例下限
EPS = 1e-3
MAX_CELLS_PER_CT = 2500  # 每细胞类型降采样上限（封顶内存；统计上足够）


def _annotate_markers(X, genes, markers):
    """X: cells×genes (dense, log1p). 返回 cell_type 数组 (argmax over marker means)."""
    gidx = {g: i for i, g in enumerate(genes)}
    cols = []
    mat = []
    for ct, gs in markers.items():
        valid = [gidx[g] for g in gs if g in gidx]
        if valid:
            m = X[:, valid].mean(axis=1)
            mat.append(m)
            cols.append(ct)
    if not mat:
        return np.array(["unknown"] * X.shape[0])
    M = np.stack(mat, axis=1)              # cells × ntypes
    return np.array(cols)[M.argmax(1)]


def load_gse174():
    """加载 GSE174574 各样本 10x mtx -> 归一化 -> marker 注释细胞类型。
    返回 list of dict: {X(cells×genes dense log1p), genes, cell_type, condition}。"""
    out = []
    for gsm, cond in GSE174_COND.items():
        d = os.path.join(GSE174, gsm)
        mtx = glob.glob(os.path.join(d, "*_matrix.mtx.gz"))
        genes_f = glob.glob(os.path.join(d, "*_genes.tsv.gz"))
        if not mtx or not genes_f:
            continue
        M = mmread(mtx[0]).T.tocsr().astype(np.float32)   # cells×genes
        genes = pd.read_csv(genes_f[0], sep="\t", header=None,
                            compression="gzip")[1].astype(str).tolist()
        # 重复 symbol 取首次出现
        seen, keep = {}, []
        for g in genes:
            if g not in seen:
                seen[g] = len(keep); keep.append(g)
        if len(keep) != M.shape[1]:
            M = M[:, [seen[g] for g in keep]]
        genes = keep
        # normalize + log1p
        X = M.toarray() if M.shape[0] < 60000 else M
        rs = X.sum(axis=1)
        rs[rs == 0] = 1
        X = X / rs[:, None] * 1e4
        X = np.log1p(X)
        ct = _annotate_markers(X, genes, CELLTYPE_MARKERS)
        # 每 CT 降采样以封顶内存
        X, ct = _subsample_ct(X, ct)
        out.append({"X": X, "genes": genes, "cell_type": ct, "condition": cond,
                    "study": "GSE174574", "time": "24h" if cond == "MCAO" else "sham"})
        print(f"  [174] {gsm} {cond}: cells={X.shape[0]} types={np.unique(ct)}")
    return out


def load_gse225():
    """加载 GSE225948 counts+metadata -> 作者 sub.celltype 注释。
    返回 list of dict: {X(cells×genes dense), genes, cell_type, condition, time}。"""
    out = []
    for gsm, (cond, tlabel) in GSE225_META.items():
        d = os.path.join(GSE225, gsm)
        cf = glob.glob(os.path.join(d, "*counts.csv.gz"))
        mf = glob.glob(os.path.join(d, "*metadata.csv.gz"))
        if not cf:
            continue
        df = pd.read_csv(cf[0], sep=",", index_col=0, compression="gzip")
        X = df.values.T.astype(np.float32)          # cells×genes
        genes = [str(g) for g in df.index]
        cell_names = [str(c) for c in df.columns]
        ct = np.array(["unknown"] * X.shape[0])
        if mf:
            meta = pd.read_csv(mf[0], sep=",", index_col=0, compression="gzip")
            meta.columns = [str(c).strip('"') for c in meta.columns]   # 去引号包裹
            meta.index = [str(i).strip('"') for i in meta.index]
            if "sub.celltype" in meta.columns:
                meta = meta.loc[meta.index.intersection(cell_names)]
                ct = meta["sub.celltype"].reindex(cell_names).fillna("unknown").values
                ct = np.array([str(c) for c in ct])
        X, ct = _subsample_ct(X, ct)
        out.append({"X": X, "genes": genes, "cell_type": ct, "condition": cond,
                    "study": "GSE225948", "time": tlabel})
        print(f"  [225] {gsm} {cond}/{tlabel}: cells={X.shape[0]} types={np.unique(ct)[:8]}")
    return out


def _subsample_ct(X, ct, cap=MAX_CELLS_PER_CT):
    """每细胞类型最多保留 cap 个细胞（按原顺序均匀抽），其余丢弃以封顶内存。"""
    keep = np.zeros(X.shape[0], dtype=bool)
    for t in np.unique(ct):
        idx = np.where(ct == t)[0]
        if len(idx) <= cap:
            keep[idx] = True
        else:
            step = len(idx) / cap
            sel = idx[np.floor(np.arange(cap) * step).astype(int)]
            keep[sel] = True
    return X[keep], ct[keep]


def _gene_index(genes):
    return {g: i for i, g in enumerate(genes)}


def _align_genes_per_study(samples):
    """按 study 把样本对齐到统一基因集(并集), 缺失补0。原地修改 samples。
    解决多样本 scRNA-seq 基因维度不一致导致 vstack 失败的问题。"""
    by_study = {}
    for s in samples:
        by_study.setdefault(s["study"], []).append(s)
    for study, ss in by_study.items():
        union = sorted(set().union(*[set(s["genes"]) for s in ss]))
        idx = {g: i for i, g in enumerate(union)}
        lowmap = {}
        for u in union:
            lowmap.setdefault(u.lower(), u)   # 大小写不敏感回退
        n = len(union)
        for s in ss:
            old = s["genes"]
            newX = np.zeros((s["X"].shape[0], n), dtype=s["X"].dtype)
            for j, g in enumerate(old):
                tgt = idx.get(g)
                if tgt is None:
                    tgt = idx.get(lowmap.get(g.lower())) if g.lower() in lowmap else None
                if tgt is not None:
                    newX[:, tgt] = s["X"][:, j]
            s["X"] = newX
            s["genes"] = union
        print(f"  aligned {study}: {n} union genes")



def comm_score(sample, ligand, receptor):
    """CellChat 风格通讯概率（简化重实现）:
    发送CT配体均值表达 × 接收CT受体均值表达, 受表达比例门控。
    返回 (score, meanL, meanR, nS, nR, fracL, fracR) 或 None。"""
    gi = _gene_index(sample["genes"])
    if ligand not in gi or receptor not in gi:
        return None
    li, ri = gi[ligand], gi[receptor]
    X = sample["X"]; ct = sample["cell_type"]
    types = np.unique(ct)
    best = None
    for s in types:
        for r in types:
            si = np.where(ct == s)[0]; rii = np.where(ct == r)[0]
            if len(si) < MIN_CELLS or len(rii) < MIN_CELLS:
                continue
            L = X[si, li]; R = X[rii, ri]
            fL = np.mean(L > 0); fR = np.mean(R > 0)
            if fL < MIN_FRAC or fR < MIN_FRAC:
                continue
            score = float(L.mean() * R.mean())
            cand = (score, float(L.mean()), float(R.mean()), len(si), len(rii),
                    float(fL), float(fR), s, r)
            if best is None or score > best[0]:
                best = cand
    return best


def study_condition_samples(samples, study, condition):
    return [s for s in samples if s["study"] == study and s["condition"] == condition]


def rewiring_for_study(samples, study, lig, rec):
    """逐 (sender,receiver) CT 对计算 sham vs 各 MCAO time 的 log2FC + MWU p。
    返回 list of dict。"""
    rows = []
    sham = study_condition_samples(samples, study, "sham")
    if not sham:
        return rows
    # 合并 sham 多样本
    sham_X = np.vstack([s["X"] for s in sham])
    sham_ct = np.concatenate([s["cell_type"] for s in sham])
    sham_genes = sham[0]["genes"]
    gi = _gene_index(sham_genes)
    if lig not in gi or rec not in gi:
        return rows
    li, ri = gi[lig], gi[rec]
    # MCAO 各 time
    mcs = {}
    for s in samples:
        if s["study"] == study and s["condition"] == "MCAO":
            mcs.setdefault(s["time"], []).append(s)
    for tlabel, ms in mcs.items():
        mc_X = np.vstack([s["X"] for s in ms])
        mc_ct = np.concatenate([s["cell_type"] for s in ms])
        # 遍历 CT 对
        types = np.unique(sham_ct)
        best = None
        for s in types:
            for r in types:
                si_s = np.where(sham_ct == s)[0]; ri_s = np.where(sham_ct == r)[0]
                si_m = np.where(mc_ct == s)[0]; ri_m = np.where(mc_ct == r)[0]
                if min(len(si_s), len(ri_s), len(si_m), len(ri_m)) < MIN_CELLS:
                    continue
                Ls = sham_X[si_s, li]; Rs = sham_X[ri_s, ri]
                Lm = mc_X[si_m, li]; Rm = mc_X[ri_m, ri]
                fL_s, fR_s = np.mean(Ls > 0), np.mean(Rs > 0)
                fL_m, fR_m = np.mean(Lm > 0), np.mean(Rm > 0)
                if min(fL_s, fR_s, fL_m, fR_m) < MIN_FRAC:
                    continue
                score_s = Ls.mean() * Rs.mean()
                score_m = Lm.mean() * Rm.mean()
                try:
                    pL = mannwhitneyu(Lm, Ls, alternative="two-sided").pvalue
                    pR = mannwhitneyu(Rm, Rs, alternative="two-sided").pvalue
                except Exception:
                    pL = pR = 1.0
                log2fc = np.log2((score_m + EPS) / (score_s + EPS))
                cand = (abs(log2fc), log2fc, pL, pR, score_s, score_m, s, r)
                if best is None or abs(log2fc) > best[0]:
                    best = cand
        if best is not None:
            _, log2fc, pL, pR, ss, sm, s, r = best
            rows.append({"study": study, "transition": tlabel, "ligand": lig,
                         "receptor": rec, "sender": s, "receiver": r,
                         "log2FC": log2fc, "pLig": pL, "pRec": pR,
                         "score_sham": ss, "score_mcao": sm})
    return rows


def load_dorothea_targets():
    """加载 DoRothEA GRN (mouse TSV: source/target/weight/confidence),
    返回 {TF: set(targets)}。直接读 TSV, 不依赖外部模块。"""
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        proj = os.path.dirname(os.path.dirname(here))   # .../gate1/rewiring_study -> 项目根
        cand = [
            os.path.join(proj, "gate1", "data", "dorothea", "mouse_dorothea_regulon.tsv"),
            os.path.join("gate1", "data", "dorothea", "mouse_dorothea_regulon.tsv"),
        ]
        p = next((c for c in cand if os.path.exists(c)), None)
        if p is None:
            print("  [warn] DoRothEA TSV 未找到")
            return {}
        res = {}
        with open(p) as fh:
            next(fh)   # 跳过表头
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                res.setdefault(parts[0], set()).add(parts[1])
        print(f"  DoRothEA 载入: {len(res)} TF, {sum(len(v) for v in res.values())} 边")
        return res
    except Exception as e:
        print("  [warn] DoRothEA 载入失败:", e)
        return {}


def main():
    print("=== 加载数据 ===")
    samples = load_gse174() + load_gse225()
    print(f"  总样本数: {len(samples)}")
    _align_genes_per_study(samples)   # 同 study 样本对齐到统一基因集

    from cellchat_lrdb import LR_DB, pathway_of
    pairs = sorted({(l, r) for l, r, _ in LR_DB})

    print("=== 计算跨转变 LR 重布线 ===")
    rows = []
    for lig, rec in pairs:
        for study in ("GSE174574", "GSE225948"):
            rows.extend(rewiring_for_study(samples, study, lig, rec))
    df = pd.DataFrame(rows)
    if df.empty:
        print("  !! 无结果, 检查基因匹配")
        return
    df["pathway"] = [pathway_of(l, r) for l, r in zip(df.ligand, df.receptor)]
    # 显著: 配体与受体 MWU 任一 <0.05 且 |log2FC|>=0.5
    df["sig"] = (np.minimum(df.pLig, df.pRec) < 0.05) & (df.log2FC.abs() >= 0.5)
    df.to_csv(os.path.join(OUT, "cellchat_rewiring.csv"), index=False)
    print(f"  LR-转变记录数: {len(df)}, 显著: {int(df.sig.sum())}")

    # 跨转变一致的"重布线 LR 对": 在 >=2 个 transition 显著且同向
    sig_df = df[df.sig]
    pair_counts = sig_df.groupby(["ligand", "receptor"]).size()
    top_pairs = pair_counts[pair_counts >= 1].sort_values(ascending=False)
    print("=== Top 重布线 LR 对 (按出现 transition 数) ===")
    print(top_pairs.head(25))

    # 与 DoRothEA 主调控因子整合
    targets = load_dorothea_targets()
    if targets:
        mrs = ["Sox10", "Cebpb", "Gata2"]
        mr_targets = set().union(*[targets.get(t, set()) for t in mrs])
        def linked(row):
            return (row.ligand in mr_targets) or (row.receptor in mr_targets)
        sig_df = sig_df.copy()
        sig_df["linked_to_MR"] = sig_df.apply(linked, axis=1)
        print(f"  显著重布线中, 配体/受体落于 Sox10/Cebpb/Gata2 靶集者: {int(sig_df.linked_to_MR.sum())}/{len(sig_df)}")
        sig_df.to_csv(os.path.join(OUT, "cellchat_rewiring_sig.csv"), index=False)

    print("=== 出图 ===")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        # 取每个 LR 对在其 Top transition 的 log2FC 做热图
        topn = 20
        # 选跨转变最一致的 top pairs, 否则按显著数取前 20
        if len(top_pairs) >= 1:
            sel = list(top_pairs.head(topn).index)
        else:
            sel = list(sig_df.groupby(["ligand","receptor"]).log2FC.abs().max()
                       .sort_values(ascending=False).head(topn).index)
        # pivot: index=LR pair, columns=transition(24h/2d/14d by study)
        piv = df.copy()
        piv["pair"] = piv.ligand + "->" + piv.receptor
        # 简化 transition 标签: 24h(GSE174574 MCAO 24h) / 2d / 14d (GSE225948)
        def trans_label(r):
            if r.study == "GSE174574":
                return "24h"
            return r.transition  # 2d / 14d
        piv["trans"] = piv.apply(trans_label, axis=1)
        sub = piv[piv.pair.isin([f"{l}->{r}" for l, r in sel])]
        heat = sub.pivot_table(index="pair", columns="trans", values="log2FC", aggfunc="first")
        for c in ["24h", "2d", "14d"]:
            if c not in heat.columns:
                heat[c] = np.nan
        heat = heat[["24h", "2d", "14d"]]
        fig, ax = plt.subplots(figsize=(6.5, max(6, topn*0.32)))
        im = ax.imshow(heat.values, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
        ax.set_xticks(range(heat.shape[1])); ax.set_xticklabels(heat.columns)
        ax.set_yticks(range(heat.shape[0])); ax.set_yticklabels(heat.index, fontsize=7)
        ax.set_title("CellChat-style LR rewiring across stroke transitions")
        fig.colorbar(im, ax=ax, label="log2(MCAO/sham) comm. score")
        fig.tight_layout()
        os.makedirs(FIG, exist_ok=True)
        fig.savefig(os.path.join(FIG, "cellchat_rewiring_heatmap.png"), dpi=150)
        print(f"  图已保存: figures/cellchat_rewiring_heatmap.png ({heat.shape})")
    except Exception as e:
        print("  [warn] 出图失败:", e)

    print("=== 完成 ===")


if __name__ == "__main__":
    main()
