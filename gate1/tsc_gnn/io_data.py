"""io_data.py — Phase 2 数据加载：log 表达 + HVG + 连续状态 score。

与 gate1.preprocessing.prep_for_gate 的区别：这里**不做 z-score 跨研究校正**
（z-score 会消除用于 GRN 相关推断的生物学变异），只用 per-study log1p。
状态 score 在 log 尺度上计算（DAM / 炎症 连续打分）。
"""
import os
import sys
import numpy as np
import scanpy as sc
from scipy import sparse as sp

# 让 tsc_gnn 能 import gate1 包（运行目录为 gate1/）
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
from gate1 import preprocessing as P


def load_phase2_data(data_root="data", cohorts=None, n_top_genes=1000):
    """加载整合时间序列 → 返回 Phase 2 就绪 dict。

    首次加载较慢（GSE174574 需网络解析 series matrix），结果缓存到
    {data_root}/phase2_cache.npz，后续秒级读盘（确定性，可复现）。

    返回：
        X      : (n_cells, G) float32, per-study log1p 表达（HVG 子集）
        genes  : (G,) 基因名
        dam    : (n_cells,) DAM 连续 score
        infl   : (n_cells,) 炎症连续 score
        state  : (n_cells,) 复合状态 z（0.5*z(dam)+0.5*z(infl)）
        time_label : (n_cells,) str，{'24h','2d','14d','sham'}
        n_cells
    """
    cohorts = cohorts or {}
    cache_key = "phase2_cache_%s_%d.npz" % (
        "_".join(f"{k}={v}" for k, v in sorted(cohorts.items())), n_top_genes)
    cache_path = os.path.join(data_root, cache_key)
    if os.path.exists(cache_path):
        print(f"  [phase2] 命中缓存 {cache_path}")
        with np.load(cache_path, allow_pickle=True) as npz:
            return dict(X=npz["X"], genes=npz["genes"], dam=npz["dam"],
                        infl=npz["infl"], state=npz["state"],
                        time_label=npz["time_label"], n_cells=int(npz["X"].shape[0]))
    adata = P.load_integrated_timeseries(data_root, cohorts=cohorts)
    # per-study log1p（不 z-score）
    studies = [s for s in adata.obs["study"].astype(str).unique()]
    out = np.zeros(adata.shape, dtype=np.float32)
    for st in studies:
        m = (adata.obs["study"].astype(str) == st).values
        sub = adata[m].copy()
        sc.pp.normalize_total(sub, target_sum=1e4)
        sc.pp.log1p(sub)
        Xs = sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
        out[m] = Xs.astype(np.float32)
    adata.X = sp.csr_matrix(out)
    # HVG（cell_ranger，无需 skmisc）
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="cell_ranger")
    adata = adata[:, adata.var["highly_variable"]].copy()
    # 状态 score（log 尺度）
    P.annotate_celltypes(adata, P.CELLTYPE_MARKERS)
    P.add_state_scores(adata, P.DAM_GENES, P.INFLAM_GENES)
    X = np.asarray(adata.X.toarray() if sp.issparse(adata.X) else adata.X, dtype=np.float32)
    genes = np.asarray(adata.var_names).astype(str)
    dam = np.asarray(adata.obs["dam_score"], dtype=float)
    infl = np.asarray(adata.obs["infl_score"], dtype=float)
    time_label = np.asarray(adata.obs["time_label"].astype(str))
    dz = (dam - dam.mean()) / (dam.std() + 1e-8)
    iz = (infl - infl.mean()) / (infl.std() + 1e-8)
    comp = 0.5 * dz + 0.5 * iz
    print(f"  [phase2] X={X.shape} HVG={genes.shape[0]} "
          f"dam genes={int((dam > 0).sum())} infl genes={int((infl > 0).sum())}")
    out = dict(X=X, genes=genes, dam=dam, infl=infl, state=comp,
               time_label=time_label, n_cells=X.shape[0])
    try:
        np.savez(cache_path, X=X, genes=genes, dam=dam, infl=infl,
                 state=comp, time_label=time_label)
        print(f"  [phase2] 已缓存 -> {cache_path}")
    except Exception as e:
        print(f"  [warn] 缓存失败（不影响结果）: {e}")
    return out
