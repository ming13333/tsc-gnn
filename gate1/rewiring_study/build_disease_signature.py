#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step D · 疾病签名构建（数据驱动，与 rewiring 框架闭环）
--------------------------------------------------------
从卒中 scRNA cache 计算 MCAO vs sham 的差异表达，得到 top up/down 基因，
作为 LINCS 反向匹配的输入签名 (disease signature)。

设计要点：
- 主签名：GSE174574 24h(MCAO) vs sham（急性期最强疾病信号）。
- 稳健性复核：GSE225948 2d vs sham（亚急性），取方向一致基因交集作为 robust 签名。
- 与 rewiring 框架闭环：报告核心 rewiring TF（Sox10 髓鞘 / Cebpb·Gata2 炎症）的
  靶 marker 落在签名哪一侧（sanity：炎症应在 UP，髓鞘应在 DOWN）。
- 输出人同源大写符号（LINCS L1000 为人细胞系）。

输出：disease_signature.json
"""
import numpy as np, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
G174 = os.path.join(HERE, "gse174_cache.npz")
G225 = os.path.join(HERE, "gse225_cache.npz")

TOPN = 100          # 每侧取的基因数（L1000CDS2 geneSet 建议 <=150）
MIN_FRAC = 0.05     # 基因须在至少 5% 细胞中表达（去噪）
MIN_ABS_LFC = 0.10  # 最小 |log fold change|（log1p 空间）

# sanity marker（人符号）：炎症应在 UP 侧、髓鞘应在 DOWN 侧
INFLAM_MARK = ["IL1B","TNF","CCL2","CCL3","CXCL10","STAT1","NFKB1","TLR2","TLR4","NOS2","PTGS2","IL6"]
MYELIN_MARK = ["PLP1","MBP","MAG","MOG","MOBP","CNP","CLDN11","MAL","MPZ","UGT8","SOX10","OLIG1","OLIG2"]


def load(npz):
    d = np.load(npz, allow_pickle=True)
    X = d["X"]                                   # (cells, genes), log1p-normalized
    genes = np.asarray([str(g) for g in d["genes"]])
    tl = np.asarray([str(t) for t in d["time_label"]])
    return X, genes, tl


def deg(X, genes, tl, grp_case, grp_ctrl):
    """Welch t + log fold change (log1p 空间). 返回 per-gene DataFrame-like dict."""
    ci = np.where(tl == grp_case)[0]
    ki = np.where(tl == grp_ctrl)[0]
    Xc, Xk = X[ci], X[ki]
    mc, mk = Xc.mean(0), Xk.mean(0)
    vc, vk = Xc.var(0, ddof=1), Xk.var(0, ddof=1)
    nc, nk = len(ci), len(ki)
    se = np.sqrt(vc / nc + vk / nk) + 1e-12
    t = (mc - mk) / se
    lfc = mc - mk
    # 表达占比（在 case 或 ctrl 任一组的表达细胞比例）
    frac = np.maximum((Xc > 0).mean(0), (Xk > 0).mean(0))
    return {"gene": genes, "t": t, "lfc": lfc, "frac": frac,
            "n_case": nc, "n_ctrl": nk}


def pick(sig, topn=TOPN):
    """按 t 统计选 top up / top down，过滤低表达与小 lfc。返回 (up_human, dn_human)."""
    g = sig["gene"]; t = sig["t"]; lfc = sig["lfc"]; frac = sig["frac"]
    keep = (frac >= MIN_FRAC) & (np.abs(lfc) >= MIN_ABS_LFC) & np.isfinite(t)
    idx = np.where(keep)[0]
    up = idx[lfc[idx] > 0]
    dn = idx[lfc[idx] < 0]
    up = up[np.argsort(-t[up])][:topn]
    dn = dn[np.argsort(t[dn])][:topn]
    up_h = [g[i].upper() for i in up]
    dn_h = [g[i].upper() for i in dn]
    return up_h, dn_h, {g[i].upper(): float(t[i]) for i in np.r_[up, dn]}


def sanity(up, dn):
    su, sd = set(up), set(dn)
    inflam_up = [m for m in INFLAM_MARK if m in su]
    inflam_dn = [m for m in INFLAM_MARK if m in sd]
    myelin_up = [m for m in MYELIN_MARK if m in su]
    myelin_dn = [m for m in MYELIN_MARK if m in sd]
    return {"inflam_in_UP": inflam_up, "inflam_in_DOWN": inflam_dn,
            "myelin_in_UP": myelin_up, "myelin_in_DOWN": myelin_dn}


def main():
    out = {}

    # ---- 主签名：GSE174574 24h vs sham ----
    X, genes, tl = load(G174)
    labs = sorted(set(tl))
    print(f"[GSE174574] labels={labs}  X={X.shape}")
    # 兼容 time_label 值（可能是 '24h'/'sham'）
    case = "24h" if "24h" in labs else [l for l in labs if l != "sham"][0]
    s174 = deg(X, genes, tl, case, "sham")
    up174, dn174, tmap174 = pick(s174)
    san174 = sanity(up174, dn174)
    print(f"  main sig: up={len(up174)} dn={len(dn174)}  case={case} (n={s174['n_case']}) vs sham (n={s174['n_ctrl']})")
    print(f"  sanity: {san174}")

    # ---- 复核签名：GSE225948 2d vs sham ----
    up225 = dn225 = None
    san225 = None
    if os.path.exists(G225):
        X2, g2, tl2 = load(G225)
        labs2 = sorted(set(tl2))
        print(f"[GSE225948] labels={labs2}  X={X2.shape}")
        if "2d" in labs2 and "sham" in labs2:
            s225 = deg(X2, g2, tl2, "2d", "sham")
            up225, dn225, _ = pick(s225)
            san225 = sanity(up225, dn225)
            print(f"  recheck sig: up={len(up225)} dn={len(dn225)}")
            print(f"  sanity: {san225}")

    # ---- robust 签名：两队列方向一致交集（若可用）----
    robust_up = robust_dn = None
    if up225 is not None:
        robust_up = sorted(set(up174) & set(up225))
        robust_dn = sorted(set(dn174) & set(dn225))
        print(f"  ROBUST intersection: up={len(robust_up)} dn={len(robust_dn)}")

    out = {
        "main": {
            "source": "GSE174574 24h(MCAO) vs sham",
            "up": up174, "dn": dn174,
            "n_up": len(up174), "n_dn": len(dn174),
            "sanity": san174,
        },
        "recheck": None if up225 is None else {
            "source": "GSE225948 2d vs sham",
            "up": up225, "dn": dn225,
            "n_up": len(up225), "n_dn": len(dn225),
            "sanity": san225,
        },
        "robust_intersection": None if robust_up is None else {
            "up": robust_up, "dn": robust_dn,
            "n_up": len(robust_up), "n_dn": len(robust_dn),
        },
        "params": {"TOPN": TOPN, "MIN_FRAC": MIN_FRAC, "MIN_ABS_LFC": MIN_ABS_LFC},
    }
    with open(os.path.join(HERE, "disease_signature.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("[done] -> disease_signature.json")
    print("  main UP head:", up174[:15])
    print("  main DN head:", dn174[:15])


if __name__ == "__main__":
    main()
