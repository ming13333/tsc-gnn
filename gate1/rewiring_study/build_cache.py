"""build_cache.py — 为时间×状态 rewiring 准备数据缓存。

流程（与 io_data.load_phase2_data 完全一致的归一化/状态协议，仅把
HVG 限制换成 DoRothEA 覆盖基因子集，以获得有意义的因果调控网络）：

1. load_integrated_timeseries(全基因, 24h/2d/14d/sham)
2. 子集到 DoRothEA(mouse, ABC) 覆盖的基因
3. 逐研究 normalize_total + log1p（与 phase2 缓存同协议）
4. 加 DAM/炎症 状态 score（复合 z）
5. 存 rewiring_doro_cache.npz

后台运行（首次加载 ~3-5 min）。
"""
import os
import sys
import numpy as np
import scanpy as sc
from scipy import sparse as sp
import pandas as pd

# 嵌套 gate1 包（真实 preprocessing）
sys.path.insert(0, "C:/D 盘/科研/虚拟敲除/gate1")
from gate1 import preprocessing as P

DATA_ROOT = "C:/D 盘/科研/虚拟敲除/gate1/data"
OUT = "C:/D 盘/科研/虚拟敲除/gate1/rewiring_study/rewiring_doro_cache.npz"
DORO = "C:/D 盘/科研/虚拟敲除/gate1/data/dorothea/mouse_dorothea_regulon.tsv"
COHORTS = {"sex": "male", "age": "W"}   # 与 phase2 缓存一致：young male，干净时间轴


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    print("[1] 加载整合时间序列（全基因）...")
    adata = P.load_integrated_timeseries(DATA_ROOT, cohorts=COHORTS)
    print(f"    全基因 adata: {adata.shape}, time_label={set(adata.obs['time_label'])}")

    # 子集到 DoRothEA 基因
    net = pd.read_csv(DORO, sep="\t")
    net = net[net["confidence"].isin(list("ABC"))]
    doro_genes = set(net["source"].str.upper()) | set(net["target"].str.upper())
    var_up = [str(g).upper() for g in adata.var_names]
    keep = [i for i, g in enumerate(var_up) if g in doro_genes]
    adata = adata[:, keep].copy()
    print(f"    子集到 DoRothEA(ABC) 覆盖基因: {adata.n_vars} 个基因")

    # 逐研究 normalize_total + log1p（与 io_data 同协议）
    studies = [s for s in adata.obs["study"].astype(str).unique()]
    n, g = adata.shape
    out = np.zeros((n, g), dtype=np.float32)
    for st in studies:
        m = (adata.obs["study"].astype(str) == st).values
        sub = adata[m].copy()
        sc.pp.normalize_total(sub, target_sum=1e4)
        sc.pp.log1p(sub)
        Xs = sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
        out[m] = Xs.astype(np.float32)
        print(f"    [{st}] 归一化 {int(m.sum())} 细胞")
    adata.X = sp.csr_matrix(out)

    # 状态 score（与 io_data 同协议）
    P.annotate_celltypes(adata, P.CELLTYPE_MARKERS)
    P.add_state_scores(adata, P.DAM_GENES, P.INFLAM_GENES)
    X = np.asarray(adata.X.toarray(), dtype=np.float32)
    genes = np.array([str(g) for g in adata.var_names])
    dam = np.asarray(adata.obs["dam_score"], dtype=float)
    infl = np.asarray(adata.obs["infl_score"], dtype=float)
    tl = np.asarray(adata.obs["time_label"].astype(str))
    dz = (dam - dam.mean()) / (dam.std() + 1e-8)
    iz = (infl - infl.mean()) / (infl.std() + 1e-8)
    state = 0.5 * dz + 0.5 * iz

    np.savez(OUT, X=X, genes=genes, dam=dam, infl=infl, state=state, time_label=tl)
    print(f"[done] 缓存 -> {OUT}")
    print(f"    X={X.shape} genes={genes.shape} time_label={dict(pd.Series(tl).value_counts())}")
    print(f"    state range [{state.min():.2f}, {state.max():.2f}] mean={state.mean():.2f}")
    # DoRothEA 边覆盖统计
    symset = set(g.upper() for g in genes)
    sub = net[(net["source"].str.upper().isin(symset)) & (net["target"].str.upper().isin(symset))]
    print(f"    DoRothEA(ABC) 落在数据基因内的边: {len(sub)} "
          f"(TF={sub['source'].nunique()}, target={sub['target'].nunique()})")


if __name__ == "__main__":
    main()
