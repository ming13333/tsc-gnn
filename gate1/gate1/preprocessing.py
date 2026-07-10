"""
preprocessing.py — 真实 scRNA 数据读取与预处理（Gate 1 真实可行性验证）。

当前聚焦于 GSE174574（鼠 MCAO 24h + sham，单时间点）。由于只有 24h 一个真实时间点，
严格"时间轴"需要整合 GSE225948（2d/14d）。本模块先用：
    · time_label = '0'(sham, 基线) / '1'(MCAO 24h, 损伤后)  —— 两"状态轴"
    · state score = 连续的 DAM / 炎症 表达打分（所有细胞共有，条件化独占信息）

这能验证"Gate 1 管线在真实鼠卒中 scRNA 上能否端到端跑通、并产生真实 rel_imp 信号"。
真实 rel_imp 无论 PASS/FAIL 都有意义：
    · 若条件化显著胜线性 → 方向成立，推进 TSC-GNN。
    · 若 ≈ 线性 → 正是 Ahlmann-Eltze 2025 现象在卒中数据上的复现，本身是可写进
      Introduction 的"问题陈述"；此时需把药物逆转 / Virtual CRISPR 排序升为主卖点或改 framing。
"""
import os
import glob
import gzip
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy import sparse as sp



# ---- 鼠细胞类型 marker（symbol）----
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
# ---- 小胶活化 / DAM 连续 score 基因 ----
DAM_GENES = ["Itgax", "Cst7", "Lpl", "Lyz2", "Cd9", "Axl", "Clec7a",
             "Spp1", "Apoe", "Trem2", "Tyrobp", "Ccl3", "Ccl4", "Cd68"]
# ---- 炎症 / 损伤响应连续 score 基因 ----
INFLAM_GENES = ["Il1b", "Tnf", "Il6", "Ccl2", "Cxcl10", "Nos2", "Ptgs2", "Ifitm1"]


def _read_10x(paths):
    """读带前缀的 10x mtx 三件套（features 第二列为 symbol）。返回 AnnData(cells×genes)。"""
    mtx = mmread(paths["matrix"]).T.tocsr()  # 10x mtx 是 genes×cells，转置为 cells×genes
    genes = pd.read_csv(paths["genes"], sep="\t", header=None, compression="gzip")
    var_names = genes[1].astype(str).tolist() if genes.shape[1] >= 2 else genes[0].astype(str).tolist()
    barcodes = pd.read_csv(paths["barcodes"], sep="\t", header=None, compression="gzip")[0].astype(str).tolist()
    ad = sc.AnnData(mtx)
    ad.var_names = var_names
    ad.obs_names = [f"{b}_{i}" for i, b in enumerate(barcodes)]
    return ad


def _find_sample_paths(gsm_dir):
    files = glob.glob(os.path.join(gsm_dir, "*_matrix.mtx.gz"))
    if not files:
        return None
    base = files[0].replace("_matrix.mtx.gz", "")
    return {
        "matrix": base + "_matrix.mtx.gz",
        "barcodes": base + "_barcodes.tsv.gz",
        "genes": base + "_genes.tsv.gz",
    }


def load_gse174574(data_root="data", min_cells=3, min_genes=200, n_top_genes=2000):
    """下载好的 GSE174574 → 预处理 anndata。obs 含 sample/condition/time_label/cell_type/
    dam_score/infl_score。X 为 log-normalized。"""
    from gate1.data_acquisition import parse_series_matrix

    gse = "GSE174574"
    samples = parse_series_matrix(gse)
    ads = []
    for s in samples:
        gsm_dir = os.path.join(data_root, gse, s["gsm"])
        paths = _find_sample_paths(gsm_dir)
        if paths is None:
            print(f"  [skip] no mtx in {gsm_dir}")
            continue
        ad = _read_10x(paths)
        ad.obs["sample"] = s["gsm"]
        ad.obs["condition"] = s["condition"]
        ad.obs["time_label"] = "0" if s["condition"] == "sham" else "1"
        ads.append(ad)
        print(f"  [load] {s['gsm']} ({s['condition']}) cells={ad.n_obs}")

    if not ads:
        raise RuntimeError("未读入任何样本，请先运行 data_acquisition.download_gse174574('data')")

    # 10x features 表常见同一 symbol 对应多个 Ensembl ID → var_names 重复；
    # concat(join="outer") 要求每个输入索引唯一，故拼接前先各自去重。
    for ad in ads:
        ad.var_names_make_unique()
    adata = sc.concat(ads, join="outer")
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # QC
    adata.var["mt"] = adata.var_names.str.startswith("mt-") | adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    adata = adata[adata.obs["n_genes_by_counts"] > min_genes].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # normalize + log + HVG
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)

    # 注释 + 状态 score
    annotate_celltypes(adata, CELLTYPE_MARKERS)
    add_state_scores(adata, DAM_GENES, INFLAM_GENES)
    return adata


def annotate_celltypes(adata, markers):
    """用 marker 平均 log 表达打分为每个细胞分配 cell_type（argmax）。"""
    cols = []
    for ct, genes in markers.items():
        valid = [g for g in genes if g in adata.var_names]
        if valid:
            col = f"sc_{ct}"
            adata.obs[col] = np.asarray(adata[:, valid].X.mean(axis=1)).ravel()
            cols.append(col)
    if cols:
        M = adata.obs[cols].values
        adata.obs["cell_type"] = np.array([cols[i].replace("sc_", "") for i in M.argmax(1)])
    else:
        adata.obs["cell_type"] = "unknown"
    return adata


def add_state_scores(adata, dam, infl):
    """添加连续 DAM / 炎症 score（条件化独占信息）。"""
    dv = [g for g in dam if g in adata.var_names]
    iv = [g for g in infl if g in adata.var_names]
    adata.obs["dam_score"] = (np.asarray(adata[:, dv].X.mean(axis=1)).ravel()
                              if dv else np.zeros(adata.n_obs))
    adata.obs["infl_score"] = (np.asarray(adata[:, iv].X.mean(axis=1)).ravel()
                               if iv else np.zeros(adata.n_obs))
    print(f"  [state] DAM genes={len(dv)} / INFLAM genes={len(iv)}")
    return adata


def load_processed(gse_id="gse174574", data_root="data"):
    if gse_id.lower() == "gse174574":
        return load_gse174574(data_root)
    raise NotImplementedError(f"暂不支持 {gse_id}；请先实现对应加载器")


# ----------------------------------------------------------------------------
# GSE225948（Anrather 2024, Nature Immunol）：CSV 表达矩阵 + 细胞注释
# 格式：counts.csv.gz = 基因(行) × 细胞(列)，首列基因名，值为已归一化浮点；
#       metadata.csv.gz = 细胞注释（tissue/sex/age/treatment/sub.celltype/parent）。
# 注意：本集只取 BRAIN 样本（避開 peripheral blood 组织混杂）；时间序列为 2d/14d。
# 与 GSE174574(24h) 整合前必须做跨研究批次校正。
# ----------------------------------------------------------------------------
# GSM -> (condition, time_label) 来自 series matrix title 解析（24h 仅来自 GSE174574）
GSE225948_SAMPLE_META = {
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


def _read_gse225948_sample(gsm_dir, gsm):
    """读单样本 CSV 表达矩阵 + metadata → AnnData(cells×genes, 已归一化浮点)。"""
    counts_file = glob.glob(os.path.join(gsm_dir, "*counts.csv.gz"))
    meta_file = glob.glob(os.path.join(gsm_dir, "*metadata.csv.gz"))
    if not counts_file:
        return None
    df = pd.read_csv(counts_file[0], sep=",", index_col=0, compression="gzip")  # (genes, cells)
    X = np.asarray(df.values.T, dtype=np.float32)  # (cells, genes)
    # 关键：保持稀疏 + float32，避免后续 concat/filter 把 ~37k×27k 稠密矩阵
    # 物化成 7.6 GiB 触发 OOM（沙箱内存有限）。
    ad = sc.AnnData(sp.csr_matrix(X))
    ad.var_names = [str(g) for g in df.index]
    ad.obs_names = [str(c) for c in df.columns]
    # metadata 对齐（index = barcode）
    if meta_file:
        try:
            meta = pd.read_csv(meta_file[0], sep=",", index_col=0, compression="gzip")
            meta = meta.loc[meta.index.intersection(ad.obs_names)]
            ad.obs = ad.obs.join(meta, how="left")
        except Exception as e:
            print(f"  [warn] {gsm} metadata 解析失败: {e}")
    cond, tlabel = GSE225948_SAMPLE_META.get(gsm, ("unknown", "unknown"))
    ad.obs["sample"] = gsm
    ad.obs["condition"] = cond
    ad.obs["time_label"] = tlabel
    ad.obs["study"] = "GSE225948"
    return ad


def load_gse225948(data_root="data", tissue="brain",
                   sex=None, age=None, condition=None, time_label=None):
    """加载 GSE225948 BRAIN 样本，按可选协变量过滤，拼接为单个 anndata。

    内存安全：协变量过滤在「逐样本」阶段完成（不匹配的样本直接跳过，不进
    大矩阵），且 sc.concat 用 join='inner' 取基因交集，避免 outer 把基因数
    撑到 27k+ 后 filter .copy() 物化成 ~7.6 GiB 稠密矩阵触发 OOM。"""
    base = os.path.join(data_root, "GSE225948")
    ads = []
    for gsm in GSE225948_SAMPLE_META:
        gsm_dir = os.path.join(base, gsm)
        if not os.path.isdir(gsm_dir):
            continue
        ad = _read_gse225948_sample(gsm_dir, gsm)
        if ad is None:
            continue
        # 逐样本协变量过滤（拼接前完成，避免不匹配样本进入大矩阵）
        keep = np.ones(ad.n_obs, dtype=bool)
        for col, val in (("tissue", tissue), ("sex", sex),
                         ("age", age), ("condition", condition),
                         ("time_label", time_label)):
            if val is not None and col in ad.obs:
                if col == "age":
                    # 周龄(W8/W10) vs 月龄(M20)：只比较首字母，兼容 W8/W10 等年轻批次
                    keep &= ad.obs[col].astype(str).str[0].eq(str(val)[0]).values
                else:
                    keep &= (ad.obs[col].astype(str) == str(val)).values
        ad = ad[keep].copy()
        if ad.n_obs == 0:
            continue
        ads.append(ad)
        print(f"  [load] {gsm} ({ad.obs['condition'].iloc[0]}) "
              f"cells={ad.n_obs} genes={ad.n_vars}")
    if not ads:
        raise RuntimeError("未读入任何 GSE225948 样本（协变量过滤后为空）")
    # join='inner' 取基因交集（~19.7k），避免 outer 撑到 27k+ 引发内存爆炸
    adata = sc.concat(ads, join="inner")
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    return adata


def align_genes(adata_a, adata_b):
    """取两 adata 基因名交集，各自子集化并返回（已对齐顺序）。"""
    ga = set(adata_a.var_names)
    gb = set(adata_b.var_names)
    common = sorted(ga & gb)
    if not common:
        raise RuntimeError("两数据集无公共基因，无法对齐")
    return adata_a[:, common].copy(), adata_b[:, common].copy(), common


def load_integrated_timeseries(data_root="data", cohorts=None):
    """整合 GSE174574(24h) + GSE225948(2d/14d) 为 24h→2d→14d 时间轴。

    cohorts: 可选过滤 GSE225948 的协变量 dict，如 {'sex':'male','age':'young'}。
    返回 anndata：obs 含 study/condition/time_label(24h|2d|14d|sham)，
    var 为两集公共基因，X 为各自归一化后拼接（GSE174574 log1p；GSE225948 原归一化）。
    批次校正(ComBat, batch=study)由 prep_for_gate 执行。"""
    # GSE174574（重新加载原始 log1p，不过 HVG 子集，以便与 GSE225948 对齐基因）
    a174_raw = load_gse174574_raw(data_root)
    a174_raw.obs["study"] = "GSE174574"
    # GSE225948（按 cohorts 过滤）
    kw = cohorts or {}
    a225 = load_gse225948(data_root, tissue="brain", **kw)
    # 对齐公共基因
    a174_a, a225_a, common = align_genes(a174_raw, a225)
    print(f"  [align] 公共基因数: {len(common)}")
    # 统一 time_label：GSE174574 的 MCAO=24h, sham=sham
    a174_a.obs["time_label"] = a174_a.obs["condition"].map(
        {"MCAO": "24h", "sham": "sham"}).fillna(a174_a.obs["time_label"])
    # 拼接
    adata = sc.concat([a174_a, a225_a], join="outer")
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.obs["time_label"] = adata.obs["time_label"].astype(str)
    return adata


def load_gse174574_raw(data_root="data", min_cells=3, min_genes=200):
    """GSE174574 原始加载（log1p 归一化，不过 HVG 子集），用于跨集基因对齐。"""
    from gate1.data_acquisition import parse_series_matrix
    gse = "GSE174574"
    samples = parse_series_matrix(gse)
    ads = []
    for s in samples:
        gsm_dir = os.path.join(data_root, gse, s["gsm"])
        paths = _find_sample_paths(gsm_dir)
        if paths is None:
            continue
        ad = _read_10x(paths)
        ad.obs["sample"] = s["gsm"]
        ad.obs["condition"] = s["condition"]
        ad.obs["time_label"] = "0" if s["condition"] == "sham" else "1"
        ads.append(ad)
    for ad in ads:
        ad.var_names_make_unique()
    adata = sc.concat(ads, join="outer")
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("mt-") | adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    adata = adata[adata.obs["n_genes_by_counts"] > min_genes].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells)
    # 注意：返回 RAW counts（不做 normalize/log1p），由 prep_for_gate 统一做，
    # 否则 GSE174574(log1p) 与 GSE225948(raw counts) 尺度不一致，HVG 会溢出。
    return adata


def prep_for_gate(adata, batch_key="study", n_top_genes=2000, method="zscore"):
    """对整合 adata 做尺度统一 + 跨研究 HVG + 批次校正 + 状态 score。

    关键：两个 study（GSE174574 与 GSE225948）在送 HVG 前必须统一到 log1p
    表达尺度，否则 seurat_v3 的 expm1 对原始 counts（max~600）溢出产生 inf。
    """
    # 0) 尺度统一（若含大值视为 raw counts，做 normalize+log1p）
    Xmax = float(adata.X.max()) if sp.issparse(adata.X) else float(np.asarray(adata.X).max())
    if Xmax > 30:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    # 1) 逐研究算 HVG（seurat_v3 在 log1p 上安全），union top 后压缩到 n_top_genes
    adata.var["highly_variable"] = False
    studies = [s for s in adata.obs[batch_key].astype(str).unique()]
    for st in studies:
        m = (adata.obs[batch_key].astype(str) == st).values
        sub = adata[m].copy()
        try:
            sc.pp.highly_variable_genes(sub, n_top_genes=min(n_top_genes, sub.n_vars),
                                        flavor="cell_ranger")
            top = sub.var["highly_variable"].values
        except Exception:
            Xs = sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
            sd = Xs.std(axis=0)
            thr = np.quantile(sd, max(0.0, 1.0 - n_top_genes / sub.n_vars))
            top = sd >= thr
        adata.var.loc[sub.var_names[top], "highly_variable"] = True
    # 若 union 超过 n_top_genes，全局再压一次
    if int(adata.var["highly_variable"].sum()) > n_top_genes:
        sub = adata[:, adata.var["highly_variable"]].copy()
        sc.pp.highly_variable_genes(sub, n_top_genes=n_top_genes, flavor="cell_ranger")
        adata.var["highly_variable"] = False
        adata.var.loc[sub.var_names[sub.var["highly_variable"].values],
                      "highly_variable"] = True
    adata = adata[:, adata.var["highly_variable"]].copy()
    print(f"  [hvg] 选定 {adata.n_vars} 个高变基因 (union of {len(studies)} studies)")
    # 2) 批次校正（研究内 z-score，保留 time 信号）
    if method == "zscore":
        studies = [s for s in adata.obs[batch_key].astype(str).unique()]
        n, g = adata.X.shape
        out = np.zeros((n, g), dtype=np.float32)
        for st in studies:
            m = (adata.obs[batch_key].astype(str) == st).values
            sub = adata.X[m].toarray() if sp.issparse(adata.X) else np.asarray(adata.X)[m]
            mu = sub.mean(axis=0)
            sd = sub.std(axis=0)
            sd[sd < 1e-8] = 1.0
            out[m] = ((sub - mu) / sd).astype(np.float32)
        adata.X = out
        print(f"  [batch] z-score 校正完成，studies={studies}")
    elif method == "combat":
        sc.pp.combat(adata, key=batch_key)
    # 3) 状态 score
    annotate_celltypes(adata, CELLTYPE_MARKERS)
    add_state_scores(adata, DAM_GENES, INFLAM_GENES)
    return adata
