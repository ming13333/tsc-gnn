"""build_cache_single.py — 纯 numpy 版（无 scanpy 依赖）构建单数据集 DoRothEA-subset cache。

为何纯 numpy：bbb_gnn 环境无 scanpy，而 build 仅需加载/归一化；rewiring 本身只用
numpy/scipy/pandas。本脚本用 scipy.mmread + pandas 直接读 10x mtx / CSV，numpy 做
normalize_total+log1p + state scores，协议与 build_cache.py 完全一致（保证与整合
rewiring 表达尺度可比）。

时间轴重映射（保证 baseline 定义跨队列一致）：
  - GSE174574：文件名 `_shamN_`/`_MCAOn_` 解析 condition；sham->"sham", MCAO->"24h"
  - GSE225948：sham 样本原标记 "2d"（采集窗口），生物学是基线 -> 重映射 "sham"
"""
import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse as sp

DATA_ROOT = "C:/D 盘/科研/虚拟敲除/gate1/data"
DORO = "C:/D 盘/科研/虚拟敲除/gate1/data/dorothea/mouse_dorothea_regulon.tsv"

# 小胶活化 / DAM 连续 score 基因（复制自 preprocessing.py）
DAM_GENES = ["Itgax", "Cst7", "Lpl", "Lyz2", "Cd9", "Axl", "Clec7a",
             "Spp1", "Apoe", "Trem2", "Tyrobp", "Ccl3", "Ccl4", "Cd68"]
INFLAM_GENES = ["Il1b", "Tnf", "Il6", "Ccl2", "Cxcl10", "Nos2", "Ptgs2", "Ifitm1"]

# GSE225948 GSM -> (condition, time_label)（复制自 preprocessing.py，全 brain 样本）
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


def normalize_log1p(X):
    """X: sparse csr (cells×genes). 行归一化到 1e4 + log1p."""
    X = X.tocsr().astype(np.float32)
    sums = np.asarray(X.sum(1)).ravel()
    sums[sums == 0] = 1.0
    X = X.multiply((1e4 / sums)[:, None]).tocsr()
    X.data = np.log1p(X.data)
    return X


def add_state_scores(Xd, genes):
    """state = 0.5*z(dam) + 0.5*z(infl)（与 build_cache 同协议）。"""
    gset = {g.upper(): i for i, g in enumerate(genes)}
    dv = [gset[g.upper()] for g in DAM_GENES if g.upper() in gset]
    iv = [gset[g.upper()] for g in INFLAM_GENES if g.upper() in gset]
    dam = Xd[:, dv].mean(1) if dv else np.zeros(Xd.shape[0])
    infl = Xd[:, iv].mean(1) if iv else np.zeros(Xd.shape[0])
    dz = (dam - dam.mean()) / (dam.std() + 1e-8)
    iz = (infl - infl.mean()) / (infl.std() + 1e-8)
    return dam, infl, (0.5 * dz + 0.5 * iz)


def load_gse174574_np(data_root):
    base = os.path.join(data_root, "GSE174574")
    dirs = sorted(glob.glob(os.path.join(base, "GSM*")))
    X_list, obs_rows = [], []
    ref_genes = None
    for d in dirs:
        fnames = set(os.listdir(d))
        mtx = next(f for f in fnames if f.endswith("_matrix.mtx.gz"))
        genes_f = next(f for f in fnames if f.endswith("_genes.tsv.gz"))
        M = mmread(os.path.join(d, mtx)).T.tocsr()
        genes = pd.read_csv(os.path.join(d, genes_f), sep="\t",
                            header=None, compression="gzip")
        gsym = genes[1].astype(str).tolist() if genes.shape[1] >= 2 \
            else genes[0].astype(str).tolist()
        if ref_genes is None:
            ref_genes = gsym
        cond = "MCAO" if "_MCAO" in mtx \
            else ("sham" if "_sham" in mtx else "unknown")
        X_list.append(M)
        n = M.shape[0]
        obs_rows.append(pd.DataFrame({
            "sample": [os.path.basename(d)] * n, "condition": [cond] * n,
            "time_label": ["24h" if cond == "MCAO" else "sham"] * n}))
    X = sp.vstack(X_list).tocsr()
    # make_unique on genes（65 个重复 symbol 加后缀）
    seen, out = set(), []
    for s in ref_genes:
        if s in seen:
            k = 1
            while f"{s}-{k}" in seen:
                k += 1
            out.append(f"{s}-{k}"); seen.add(f"{s}-{k}")
        else:
            out.append(s); seen.add(s)
    var_names = np.array(out)
    obs = pd.concat(obs_rows, ignore_index=True)
    return X, var_names, obs


def load_gse225948_np(data_root):
    base = os.path.join(data_root, "GSE225948")
    X_list, var_list, obs_list = [], [], []
    for gsm, (cond, tlabel) in GSE225948_SAMPLE_META.items():
        d = os.path.join(base, gsm)
        if not os.path.isdir(d):
            continue
        counts_f = glob.glob(os.path.join(d, "*counts.csv.gz"))
        if not counts_f:
            continue
        df = pd.read_csv(counts_f[0], sep=",", index_col=0, compression="gzip")
        X = sp.csr_matrix(df.values.T.astype(np.float32))
        genes = [str(g) for g in df.index]
        # tissue 过滤（仅 brain）
        meta_f = glob.glob(os.path.join(d, "*metadata.csv.gz"))
        tissue_ok = True
        if meta_f:
            meta = pd.read_csv(meta_f[0], sep=",", index_col=0,
                               compression="gzip")
            if "tissue" in meta.columns:
                tv = meta["tissue"].astype(str).str.lower()
                tissue_ok = ((tv == "brain").all()) if len(tv) else True
        if not tissue_ok:
            continue
        X_list.append(X)
        var_list.append(genes)
        n = X.shape[0]
        obs_list.append(pd.DataFrame({
            "sample": [gsm] * n, "condition": [cond] * n,
            "time_label": [tlabel] * n}))
    # join inner on genes
    common = None
    for g in var_list:
        s = set(g)
        common = s if common is None else (common & s)
    common = sorted(common)
    Xs = []
    for X, g in zip(X_list, var_list):
        idx = [g.index(c) for c in common]
        Xs.append(X[:, idx])
    X = sp.vstack(Xs).tocsr()
    obs = pd.concat(obs_list, ignore_index=True)
    return X, np.array(common), obs


def build(dataset, out_npz):
    _d = os.path.dirname(out_npz)
    if _d:
        os.makedirs(_d, exist_ok=True)
    print(f"[1] 加载 {dataset} ...")
    if dataset == "GSE174574":
        X, genes, obs = load_gse174574_np(DATA_ROOT)
    elif dataset == "GSE225948":
        X, genes, obs = load_gse225948_np(DATA_ROOT)
        cond = obs["condition"].astype(str)
        tl = obs["time_label"].astype(str)
        obs["time_label"] = np.where(cond == "sham", "sham", tl)
    else:
        raise ValueError(dataset)
    print(f"    加载完成: X={X.shape}, "
          f"time_label={dict(pd.Series(obs['time_label'].astype(str)).value_counts())}")

    # 子集 DoRothEA(ABC)
    net = pd.read_csv(DORO, sep="\t")
    net = net[net["confidence"].isin(list("ABC"))]
    doro_genes = set(net["source"].str.upper()) | set(net["target"].str.upper())
    gu = [g.upper() for g in genes]
    keep = [i for i, g in enumerate(gu) if g in doro_genes]
    X = X[:, keep].tocsr()
    genes = genes[keep]
    print(f"    子集 DoRothEA(ABC): {X.shape[1]} 基因, {len(keep)} 列")

    # normalize + log1p
    X = normalize_log1p(X)
    Xd = np.asarray(X.toarray(), dtype=np.float32)

    # state scores
    dam, infl, state = add_state_scores(Xd, genes)
    tl = np.asarray(obs["time_label"].astype(str))

    np.savez(out_npz, X=Xd, genes=genes, dam=dam, infl=infl,
             state=state, time_label=tl)
    print(f"[done] 缓存 -> {out_npz}")
    print(f"    X={Xd.shape} time_label={dict(pd.Series(tl).value_counts())}")
    symset = set(g.upper() for g in genes)
    sub = net[(net["source"].str.upper().isin(symset)) &
              (net["target"].str.upper().isin(symset))]
    print(f"    DoRothEA(ABC) 落在数据内边: {len(sub)} "
          f"(TF={sub['source'].nunique()}, target={sub['target'].nunique()})")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=["GSE174574", "GSE225948"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    build(args.dataset, args.out)
