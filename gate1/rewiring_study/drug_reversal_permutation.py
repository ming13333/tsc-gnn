#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step D · 药物逆转背景置换对照
------------------------------------------------
问题：robust 签名反向匹配命中已知神经保护药（vorinostat/他汀），是真实特异
      还是 L1000 top50 里本就常见这些"泛逆转"药？
方法：用与 robust 签名相同基因数(up/dn)、从表达基因池随机抽的签名，重复跑
      L1000CDS2 reverse N 次，统计每次 top50 命中已知神经保护药的数量，
      与真实 robust 命中数比较得经验 p。避免确认偏误。

输出：drug_perm_result.json
"""
import requests, json, os, time, random
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SIG = os.path.join(HERE, "disease_signature.json")
G174 = os.path.join(HERE, "gse174_cache.npz")
API = "https://maayanlab.cloud/L1000CDS2/query"
N_PERM = 20
random.seed(2026)

# 与 drug_reversal.py 同一已知药先验（子串匹配）
KNOWN_KEYS = ["minocycline", "edaravone", "fasudil", "parthenolide", "rosiglitazone",
    "pioglitazone", "dexamethasone", "curcumin", "resveratrol", "simvastatin", "lovastatin",
    "atorvastatin", "mevastatin", "rosuvastatin", "statin", "vorinostat", "trichostatin",
    "isox", "celastrol", "triptolide", "sirolimus", "rapamycin", "ibuprofen", "celecoxib",
    "wortmannin"]


def n_known(topMeta):
    # 按 unique drug 去重计数（与 drug_reversal.py 的 candidates 列表口径一致），
    # 避免同一药物在多细胞系重复出现导致虚高
    seen = set()
    c = 0
    for m in topMeta:
        d = (m.get("pert_desc") or "").lower()
        if d in seen:
            continue
        seen.add(d)
        if any(k in d for k in KNOWN_KEYS):
            c += 1
    return c


def query_reverse(up, dn):
    payload = {"data": {"upGenes": up, "dnGenes": dn},
               "config": {"aggravate": False, "searchMethod": "geneSet",
                          "share": False, "combination": False, "db-version": "latest"},
               "metadata": []}
    r = requests.post(API, json=payload, timeout=90)
    r.raise_for_status()
    return r.json().get("topMeta", [])


def main():
    sig = json.load(open(SIG, encoding="utf-8"))
    ri = sig["robust_intersection"]
    up0, dn0 = ri["up"], ri["dn"]
    nu, nd = len(up0), len(dn0)

    # 真实 robust 命中
    real_tm = query_reverse(up0, dn0)
    real_hits = n_known(real_tm)
    print(f"[real robust] up={nu} dn={nd}  known-hits in top{len(real_tm)} = {real_hits}")

    # 基因池：gse174 表达基因（大写人符号），排除签名基因
    d = np.load(G174, allow_pickle=True)
    pool = sorted(set(str(g).upper() for g in d["genes"]) - set(up0) - set(dn0))
    print(f"  gene pool size = {len(pool)}")

    perm_hits = []
    for i in range(N_PERM):
        samp = random.sample(pool, nu + nd)
        rup, rdn = samp[:nu], samp[nu:]
        try:
            tm = query_reverse(rup, rdn)
            h = n_known(tm)
        except Exception as e:
            print(f"  perm {i}: FAIL {type(e).__name__}"); continue
        perm_hits.append(h)
        print(f"  perm {i:2d}: known-hits = {h}")
        time.sleep(1)

    ph = np.array(perm_hits, dtype=float)
    emp_p = (1 + int((ph >= real_hits).sum())) / (1 + len(ph))
    res = {"real_known_hits": real_hits, "n_perm": len(ph),
           "perm_mean": float(ph.mean()), "perm_std": float(ph.std()),
           "perm_max": float(ph.max()), "perm_hits": [int(x) for x in ph],
           "empirical_p": emp_p,
           "known_keys": KNOWN_KEYS, "sig_n_up": nu, "sig_n_dn": nd}
    json.dump(res, open(os.path.join(HERE, "drug_perm_result.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    print(f"\n[RESULT] real={real_hits}  perm mean={ph.mean():.2f}±{ph.std():.2f} max={ph.max():.0f}  emp_p={emp_p:.4f}")
    print("[done] -> drug_perm_result.json")


if __name__ == "__main__":
    main()
