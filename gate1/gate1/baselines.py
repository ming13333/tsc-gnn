"""
baselines.py — Gate 1 的两个对照模型。

(1) LinearBaseline  —— 复现 Ahlmann-Eltze 2025 (Nat Methods s41592-025-02772-6) 思路：
    用【无时间/状态上下文】的线性模型预测扰动/位移效应。这是对照锚点。
(2) CrudeConditional —— 把时间 one-hot + 状态 score 拼到特征上，仍用线性/浅 MLP：
    检验"任何形式的条件化"是否能救回线性基线的失败。这是 Gate 1 的真问题。

注意：此处用 Ridge 作为线性代表（与 Ahlmann-Eltze 的线性基线一致）。
基因数 >> 细胞数，Ridge 的 L2 正则保证稳定。
"""
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.neural_network import MLPRegressor


def _safe_ridge():
    return Ridge(alpha=1.0, random_state=0, max_iter=2000)


def fit_predict_linear(X_tr, y_tr, X_te, use_mlp=False):
    """线性基线：只用表达 x 预测位移，完全不看 time/state。"""
    if use_mlp:
        model = MLPRegressor(hidden_layer_sizes=(64, 32), max_iter=500,
                             random_state=0, early_stopping=True)
    else:
        model = _safe_ridge()
    model.fit(X_tr, y_tr)
    return model.predict(X_te)


def fit_predict_conditional(feat_tr, y_tr, feat_te, use_mlp=False):
    """粗条件化：特征 = [x | onehot(time) | state_score]。"""
    if use_mlp:
        model = MLPRegressor(hidden_layer_sizes=(128, 64), max_iter=500,
                             random_state=0, early_stopping=True)
    else:
        model = _safe_ridge()
    model.fit(feat_tr, y_tr)
    return model.predict(feat_te)
