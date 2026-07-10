#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step D · 药物逆转（LINCS L1000 反向匹配）
------------------------------------------------
输入：disease_signature.json（Step D 数据驱动卒中签名，up/down 人基因）
方法：L1000CDS2 (Ma'ayan lab) signature-reversal API，aggravate=False → 找能
      "逆转"卒中签名的化合物（药物诱导表达与疾病签名反向）。
聚合：同一化合物在多个细胞系/剂量返回多条 → 按 pert_desc 聚合，取 best(min) score
      与 mean score + 命中次数；score 越小 = 逆转越强。
对照：标注命中的已知卒中神经保护/抗炎药类别（文献先验），作为 positive control。

输出：drug_reversal_result.json、drug_candidates_main.csv、drug_candidates_robust.csv
"""
import requests, json, os, time
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
SIG = os.path.join(HERE, "disease_signature.json")
API = "https://maayanlab.cloud/L1000CDS2/query"

# 已知卒中相关神经保护/抗炎药先验（小写子串匹配 pert_desc），作 positive-control 标注
KNOWN = {
    "minocycline": "抗炎神经保护(卒中临床试验)",
    "edaravone": "自由基清除(卒中批准药)",
    "fasudil": "Rho激酶抑制剂(卒中用药)",
    "parthenolide": "NF-κB抑制剂(抗炎)",
    "rosiglitazone": "PPARγ激动剂(神经保护)",
    "pioglitazone": "PPARγ激动剂(神经保护)",
    "dexamethasone": "糖皮质激素(抗炎)",
    "curcumin": "抗炎抗氧化(神经保护)",
    "resveratrol": "抗氧化(神经保护)",
    "simvastatin": "他汀(神经保护)",
    "lovastatin": "他汀(神经保护)",
    "atorvastatin": "他汀(神经保护)",
    "vorinostat": "HDAC抑制剂(神经保护)",
    "trichostatin": "HDAC抑制剂(神经保护)",
    "isox": "HDAC6抑制剂",
    "wortmannin": "PI3K抑制剂(炎症)",
    "celastrol": "NF-κB抑制剂(抗炎)",
    "triptolide": "免疫抑制抗炎",
    "sirolimus": "mTOR抑制剂(自噬/神经保护)",
    "rapamycin": "mTOR抑制剂(自噬/神经保护)",
    "ibuprofen": "NSAID(抗炎)",
    "sc-560": "COX-1抑制剂(抗炎)",
    "celecoxib": "COX-2抑制剂(抗炎)",
    "ammonium pyrrolidinedithiocarbamate": "NF-κB抑制剂(抗炎)",
}


def known_hit(desc):
    d = (desc or "").lower()
    for k, v in KNOWN.items():
        if k in d:
            return f"{k} · {v}"
    return ""


def query_reverse(up, dn, tag):
    payload = {
        "data": {"upGenes": up, "dnGenes": dn},
        "config": {"aggravate": False, "searchMethod": "geneSet",
                   "share": False, "combination": False, "db-version": "latest"},
        "metadata": [{"key": "tag", "value": tag}],
    }
    r = requests.post(API, json=payload, timeout=90)
    r.raise_for_status()
    j = r.json()
    return j.get("topMeta", []), j.get("shareId", "")


def aggregate(topMeta):
    """按 pert_desc 聚合多细胞系/剂量记录。score 越小逆转越强。"""
    by = defaultdict(list)
    for m in topMeta:
        desc = m.get("pert_desc") or m.get("pert_id") or "?"
        by[desc].append(m)
    rows = []
    for desc, recs in by.items():
        scores = [r.get("score", 1.0) for r in recs]
        cells = sorted(set(r.get("cell_id", "?") for r in recs))
        pid = recs[0].get("pert_id", "")
        rows.append({
            "drug": desc, "pert_id": pid,
            "best_score": round(min(scores), 4),
            "mean_score": round(sum(scores) / len(scores), 4),
            "n_hits": len(recs),
            "cells": ";".join(cells[:6]),
            "known": known_hit(desc),
        })
    rows.sort(key=lambda x: (x["best_score"], x["mean_score"]))
    return rows


def write_csv(rows, path):
    import csv
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["rank", "drug", "pert_id", "best_score",
                                          "mean_score", "n_hits", "cells", "known"])
        w.writeheader()
        for i, r in enumerate(rows, 1):
            r2 = dict(r); r2["rank"] = i
            w.writerow(r2)


def main():
    sig = json.load(open(SIG, encoding="utf-8"))
    out = {}
    runs = [("main", sig["main"]["up"], sig["main"]["dn"])]
    if sig.get("robust_intersection"):
        ri = sig["robust_intersection"]
        if ri["n_up"] >= 5 and ri["n_dn"] >= 5:
            runs.append(("robust", ri["up"], ri["dn"]))

    for tag, up, dn in runs:
        print(f"[{tag}] query reverse: up={len(up)} dn={len(dn)}")
        tm, sid = query_reverse(up, dn, tag)
        rows = aggregate(tm)
        out[tag] = {"n_raw": len(tm), "shareId": sid,
                    "n_unique_drugs": len(rows), "candidates": rows}
        # 打印 top 与命中的已知药
        print(f"  raw={len(tm)} unique_drugs={len(rows)}  shareId={sid}")
        print(f"  --- top10 reversal candidates ({tag}) ---")
        for i, r in enumerate(rows[:10], 1):
            kn = f"  [KNOWN: {r['known']}]" if r["known"] else ""
            print(f"   {i:2d}. {r['drug']:36s} best={r['best_score']:.4f} n={r['n_hits']}{kn}")
        hits = [r for r in rows if r["known"]]
        print(f"  known-drug hits in list: {len(hits)}")
        for r in hits[:12]:
            print(f"     * rank {rows.index(r)+1:3d}  {r['drug']:34s} {r['known']}")
        write_csv(rows, os.path.join(HERE, f"drug_candidates_{tag}.csv"))
        time.sleep(1)

    json.dump(out, open(os.path.join(HERE, "drug_reversal_result.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    print("[done] -> drug_reversal_result.json + drug_candidates_*.csv")


if __name__ == "__main__":
    main()
