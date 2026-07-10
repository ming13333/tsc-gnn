"""
evaluate.py — Gate 1 的度量与判定逻辑。

核心问题：在 24h→2d→14d 轴上，"粗条件化模型"是否显著优于"时间盲线性基线"？
若是 → Gate 1 PASS → 值得投入 TSC-GNN 深核（Phase 2）。
若否 → 方向证伪 → 立刻 pivot（药物逆转升主卖点 / 空间组角度）。

判定规则（保守、防假阳性）：
    PASS 当且仅当：
      (a) 相对 MSE 改善 rel_imp = (mse_lin − mse_cond)/mse_lin ≥ 0.10
      (b) bootstrap 1000 次重采样的 rel_imp 95% CI 下界 > 0（显著）
      (c) 条件化模型的平均每基因相关 corr_cond > 线性 corr_lin
"""
import numpy as np
from scipy import stats


def mse(pred, true, mask):
    p = np.asarray(pred)[mask]
    t = np.asarray(true)[mask]
    return float(np.mean((p - t) ** 2))


def per_gene_corr(pred, true, mask):
    """每个基因在带 target 的细胞上的 Pearson 相关，返回均值。"""
    p = np.asarray(pred)[mask]
    t = np.asarray(true)[mask]
    if p.shape[0] < 3:
        return np.nan
    cors = []
    for g in range(p.shape[1]):
        if np.std(p[:, g]) < 1e-9 or np.std(t[:, g]) < 1e-9:
            continue
        r, _ = stats.pearsonr(p[:, g], t[:, g])
        cors.append(r)
    return float(np.nanmean(cors)) if cors else np.nan


def _rel_imp_single(pred_lin, pred_cond, true, mask):
    ml = mse(pred_lin, true, mask)
    mc = mse(pred_cond, true, mask)
    if ml <= 1e-12:
        return 0.0
    return (ml - mc) / ml


def bootstrap_rel_improvement(cell_mse_lin, cell_mse_cond, n_boot=1000, seed=0):
    """在【每细胞 MSE】上做 bootstrap（等价于对 cells×genes 的 MSE 重采样，
    但先把基因维折叠成每细胞标量，避免每次迭代搬运 (n_cells × n_genes) 大数组）。
    数学等价：MSE = mean_cells(mean_genes(...))，bootstrap 重采样细胞后取均值一致。"""
    rng = np.random.default_rng(seed)
    n = len(cell_mse_lin)
    vals = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        s = rng.choice(n, size=n, replace=True)
        ml = cell_mse_lin[s].mean()
        mc = cell_mse_cond[s].mean()
        vals[b] = (ml - mc) / ml if ml > 1e-12 else 0.0
    return float(np.mean(vals)), float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def evaluate_gate(pred_lin, pred_cond, true, mask, n_boot=1000):
    ml = mse(pred_lin, true, mask)
    mc = mse(pred_cond, true, mask)
    rel_imp = (ml - mc) / ml if ml > 1e-12 else 0.0
    corr_lin = per_gene_corr(pred_lin, true, mask)
    corr_cond = per_gene_corr(pred_cond, true, mask)
    # 折叠成每细胞 MSE 后做向量化 bootstrap（避免大数组反复搬运）
    pl = np.asarray(pred_lin)[np.asarray(mask)]
    pc = np.asarray(pred_cond)[np.asarray(mask)]
    tt = np.asarray(true)[np.asarray(mask)]
    cell_mse_lin = np.mean((pl - tt) ** 2, axis=1)
    cell_mse_cond = np.mean((pc - tt) ** 2, axis=1)
    mean_imp, ci_lo, ci_hi = bootstrap_rel_improvement(cell_mse_lin, cell_mse_cond, n_boot)
    verdict_pass = (rel_imp >= 0.10) and (ci_lo > 0) and (corr_cond > corr_lin)
    return {
        "mse_linear": ml,
        "mse_conditional": mc,
        "rel_improvement": rel_imp,
        "boot_mean": mean_imp,
        "boot_ci_lo": ci_lo,
        "boot_ci_hi": ci_hi,
        "corr_linear": corr_lin,
        "corr_conditional": corr_cond,
        "verdict_pass": bool(verdict_pass),
    }


def format_report(res, title="Gate 1"):
    lines = []
    lines.append(f"===== {title} =====")
    lines.append(f"MSE  linear      : {res['mse_linear']:.4f}")
    lines.append(f"MSE  conditional : {res['mse_conditional']:.4f}")
    lines.append(f"Rel. improvement : {res['rel_improvement']*100:.1f}%")
    lines.append(f"Bootstrap mean  : {res['boot_mean']*100:.1f}%  "
                 f"(95% CI {res['boot_ci_lo']*100:.1f}% .. {res['boot_ci_hi']*100:.1f}%)")
    lines.append(f"Corr linear      : {res['corr_linear']:.3f}")
    lines.append(f"Corr conditional : {res['corr_conditional']:.3f}")
    verdict_str = ("PASS (条件化显著优于线性，值得推进 TSC-GNN)"
                   if res["verdict_pass"] else
                   "FAIL (条件化未显著救回线性，需 pivot)")
    lines.append("VERDICT         : " + verdict_str)
    return "\n".join(lines)
