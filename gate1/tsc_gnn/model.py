"""model.py — TSC-GNN 消息传递 + Ridge 读头，以及扁平 baseline。

图消息传递（固定扩散，无 learned 图权重 → 无梯度 bug，完全可复现）：
  F0 = [x (表达), p (扰动指示)]            (n, G, 2)
  H  = [F0, A_i F0, A_i^2 F0, ...]         (n, G, 2*(K+1))   # A_i 状态条件化
节点嵌入 = H + 上下文列(state, time)      (n, G, d_total)
读头 = 每基因共享线性： Δ_hat = E @ W       (d_total × 1)

扁平 baseline（Ahlmann-Eltze / coarse conditional）：
  flat = [x, p] 或 [x, p, state, time] → Ridge → Δ
  —— 无图传播，对未见过的扰动基因无法泛化。
"""
import numpy as np
from scipy import sparse as sp


def graph_propagate(x, p, z, grn, K=2, gamma=1.0, n_bins=12):
    """状态条件化图上的消息传递（向量化按 z 分箱）。

    返回 H: (n, G, 2*(K+1))，hop h 的 [x通道, p通道] 在通道索引 2h / 2h+1。
    p 通道的 hop1 = A_i p（扰动一跳传播），hop2 = A_i^2 p（二跳）。
    """
    n, G = x.shape
    F0 = np.stack([x, p], axis=-1).astype(np.float32)  # (n, G, 2)
    zc = np.clip(z, -4.0, 4.0)
    edges = np.linspace(zc.min(), zc.max() + 1e-9, n_bins + 1)
    binid = np.clip(np.digitize(zc, edges) - 1, 0, n_bins - 1)
    H = F0.copy()
    cur = F0
    for _ in range(K):
        nxt = np.zeros_like(cur)
        for b in range(n_bins):
            m = binid == b
            nb = int(m.sum())
            if nb == 0:
                continue
            zc_b = 0.5 * (edges[b] + edges[b + 1])
            Ai = grn["build_adj"](zc_b, gamma)  # (G, G) csr
            cur_b = cur[m].transpose(1, 0, 2).reshape(G, nb * 2)  # (G, nb*2)
            out = Ai @ cur_b  # (G, nb*2)
            nxt[m] = out.reshape(G, nb, 2).transpose(1, 0, 2)
        H = np.concatenate([H, nxt], axis=-1)
        cur = nxt
    return H


def _standardize(fit, apply=None, params=None):
    if params is None:
        mu = fit.mean(0)
        sd = fit.std(0)
        sd[sd < 1e-8] = 1.0
        params = (mu, sd)
    mu, sd = params
    if apply is None:
        apply = fit
    return (apply - mu) / sd, params


def _ridge_solve(Xmat, ymat, lam=1.0):
    d = Xmat.shape[1]
    XtX = Xmat.T @ Xmat
    XtX += lam * np.eye(d)
    W = np.linalg.solve(XtX, Xmat.T @ ymat)
    return W


def _gnn_test_mse(xt, pt, st, tt, dt, xv, pv, sv, tv, dv, grn,
                  K, gamma, include_graph, include_state, lam=1.0):
    Ht = graph_propagate(xt, pt, st, grn, K, gamma)
    Hv = graph_propagate(xv, pv, sv, grn, K, gamma)
    if not include_graph:
        Ht = Ht[:, :, :2]
        Hv = Hv[:, :, :2]
    if include_state:
        ct = np.concatenate([np.broadcast_to(st[:, None, None], (Ht.shape[0], Ht.shape[1], 1)),
                             np.broadcast_to(tt[:, None, :], (Ht.shape[0], Ht.shape[1], tt.shape[1]))], axis=2)
        cv = np.concatenate([np.broadcast_to(sv[:, None, None], (Hv.shape[0], Hv.shape[1], 1)),
                             np.broadcast_to(tv[:, None, :], (Hv.shape[0], Hv.shape[1], tv.shape[1]))], axis=2)
        Et, Ev = np.concatenate([Ht, ct], axis=2), np.concatenate([Hv, cv], axis=2)
    else:
        Et, Ev = Ht, Hv
    d_in = Et.shape[2]
    Etr = Et.reshape(-1, d_in)
    Evl = Ev.reshape(-1, d_in)
    Etr_z, pp = _standardize(Etr)
    Evl_z, _ = _standardize(Evl, params=pp)
    W = _ridge_solve(Etr_z, dt.reshape(-1), lam=lam)
    dv_hat = (Evl_z @ W).reshape(dv.shape)
    return ((dv_hat - dv) ** 2).mean(axis=1)


def _linear_test_mse(xt, pt, st, tt, dt, xv, pv, sv, tv, dv,
                     include_state, lam=1.0):
    if include_state:
        ft = np.concatenate([xt, pt, st[:, None], tt], axis=1)
        fv = np.concatenate([xv, pv, sv[:, None], tv], axis=1)
    else:
        ft = np.concatenate([xt, pt], axis=1)
        fv = np.concatenate([xv, pv], axis=1)
    ft_z, pp = _standardize(ft)
    fv_z, _ = _standardize(fv, params=pp)
    W = _ridge_solve(ft_z, dt, lam=lam)  # (d_flat, G)
    dv_hat = fv_z @ W  # (n_test, G)
    return ((dv_hat - dv) ** 2).mean(axis=1)


def compute_test_mse(X, p, state, time_onehot, delta, train_mask, test_mask,
                     grn, K=2, gamma=1.0, model="tscgnn",
                     include_graph=True, include_state=True, linear_include_state=True,
                     lam=1.0):
    """返回测试集每细胞 MSE (n_test,)。"""
    xt, pt, st, tt, dt = (X[train_mask], p[train_mask], state[train_mask],
                          time_onehot[train_mask], delta[train_mask])
    xv, pv, sv, tv, dv = (X[test_mask], p[test_mask], state[test_mask],
                          time_onehot[test_mask], delta[test_mask])
    if model == "linear":
        return _linear_test_mse(xt, pt, st, tt, dt, xv, pv, sv, tv, dv,
                                include_state=linear_include_state, lam=lam)
    return _gnn_test_mse(xt, pt, st, tt, dt, xv, pv, sv, tv, dv, grn,
                         K, gamma, include_graph, include_state, lam=lam)
