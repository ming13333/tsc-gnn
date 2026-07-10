"""train_eval.py — 训练/评估/消融 + bootstrap 显著性。

比较对象：
  · linear_pure       : 扁平 Ridge on [x, p]                    （Ahlmann-Eltze 2025 纯线性）
  · linear_coarse     : 扁平 Ridge on [x, p, state, time]       （粗条件化，Gate 1 同款）
  · tscgnn            : TSC-GNN 全条件化（图 + 状态）            （深核）
  · tscgnn_nograph    : 去图（仅 [x,p] 传播 0 步 + 状态上下文）   （消融：图贡献）
  · tscgnn_nostate    : 去状态（gamma=0 无门控 + 无状态上下文）   （消融：状态贡献）

核心指标：rel_imp = 1 − mse_model / mse_baseline （对测试集每细胞 MSE）。
bootstrap 1000 次细胞重采样得 CI；verdict = rel_imp≥0.10 且 CI 下界>0。
"""
import numpy as np
from tsc_gnn import model as M


def bootstrap_rel_improvement(mse_a, mse_b, n_boot=1000, seed=0):
    """a=baseline, b=model。返回 (mean_rel, ci_low, ci_high)。"""
    rng = np.random.default_rng(seed)
    n = len(mse_a)
    rels = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        ma = mse_a[idx].mean()
        mb = mse_b[idx].mean()
        rels[i] = (1.0 - mb / ma) if ma > 0 else 0.0
    return float(rels.mean()), float(np.percentile(rels, 2.5)), float(np.percentile(rels, 97.5))


def _verdict(rel, ci_low):
    return (rel >= 0.10) and (ci_low > 0.0)


def evaluate_all(X, p, state, time_onehot, delta, train_mask, test_mask,
                 grn, K=2, gamma=1.0, n_boot=1000, seed=0):
    base = dict(X=X, p=p, state=state, time_onehot=time_onehot, delta=delta,
                train_mask=train_mask, test_mask=test_mask, grn=grn, K=K)
    # 各模型测试集每细胞 MSE
    mse_lin_pure = M.compute_test_mse(**base, model="linear", gamma=gamma,
                                      linear_include_state=False)
    mse_lin_coarse = M.compute_test_mse(**base, model="linear", gamma=gamma,
                                        linear_include_state=True)
    mse_full = M.compute_test_mse(**base, model="tscgnn", gamma=gamma,
                                  include_graph=True, include_state=True)
    mse_nograph = M.compute_test_mse(**base, model="tscgnn", gamma=gamma,
                                     include_graph=False, include_state=True)
    mse_nostate = M.compute_test_mse(**base, model="tscgnn", gamma=0.0,
                                     include_graph=True, include_state=False)
    mse_bridge = dict(linear_pure=mse_lin_pure, linear_coarse=mse_lin_coarse,
                      tscgnn=mse_full, tscgnn_nograph=mse_nograph,
                      tscgnn_nostate=mse_nostate)

    def rep(name, mb, ma):
        rel, lo, hi = bootstrap_rel_improvement(ma, mb, n_boot, seed)
        return dict(name=name, rel=rel, ci_low=lo, ci_high=hi,
                    mse_model=float(mb.mean()), mse_base=float(ma.mean()),
                    verdict=_verdict(rel, lo))

    results = {
        "vs_linear_pure": rep("TSC-GNN vs Ahlmann-Eltze 纯线性", mse_full, mse_lin_pure),
        "vs_linear_coarse": rep("TSC-GNN vs 粗条件化(扁平+state)", mse_full, mse_lin_coarse),
        "abl_nograph": rep("TSC-GNN 去图 vs 全条件化", mse_nograph, mse_full),
        "abl_nostate": rep("TSC-GNN 去状态 vs 全条件化", mse_nostate, mse_full),
    }
    return results, mse_bridge


def format_report(results, title=""):
    lines = [f"===== {title} ====="]
    for k, r in results.items():
        v = "PASS" if r["verdict"] else "FAIL"
        lines.append(
            f"  [{r['name']}]\n"
            f"    MSE model={r['mse_model']:.4f}  base={r['mse_base']:.4f}\n"
            f"    Rel. improvement={r['rel']*100:.1f}%  "
            f"CI=[{r['ci_low']*100:.1f}%, {r['ci_high']*100:.1f}%]  -> {v}")
    return "\n".join(lines)
