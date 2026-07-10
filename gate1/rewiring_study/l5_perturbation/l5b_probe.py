#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Probe v4: dump all pert_desc from SOX10 query; try iLINCS correct endpoint."""
import json, requests

DORO = r"C:\D 盘\科研\虚拟敲除\gate1\data\dorothea\human_dorothea_regulon.tsv"
API = "https://maayanlab.cloud/L1000CDS2/query"

tgt = {}
with open(DORO, encoding="utf-8") as f:
    next(f)
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) < 2: continue
        tgt.setdefault(p[0], set()).add(p[1])
sox10 = sorted(tgt.get("SOX10", set()))
sox10_set = set(sox10)
neutral = ['HBB','HBA1','INS','ALB','TFRC','FGB','HBA2','HPX','SERPINA1','APOA1']
dn = [g for g in neutral if g not in sox10_set][:3]

payload = {"data": {"upGenes": sox10, "dnGenes": dn},
           "config": {"aggravate": True, "searchMethod": "geneSet",
                      "share": False, "combination": False, "db-version": "latest"},
           "metadata": [{"key": "tag", "value": "p"}]}
r = requests.post(API, json=payload, timeout=120)
tm = r.json().get("topMeta", [])
print(f"topMeta len={len(tm)}; dumping pert_desc (first 50):")
for i, m in enumerate(tm):
    desc = m.get("pert_desc", "")
    pid = m.get("pert_id", "")
    mark = " <<SOX10" if "SOX10" in str(desc).upper() or "SOX10" in str(pid).upper() else ""
    print(f"  {i+1:2d}. desc={desc!r} id={pid}{mark}")

# iLINCS correct endpoint
print("\n== iLINCS /api/ilincsSearch/LINCS?query=SOX10 ==")
try:
    ri = requests.get("https://www.ilincs.org/api/ilincsSearch/LINCS?query=SOX10", timeout=60)
    print("status:", ri.status_code, "len:", len(ri.text))
    print(ri.text[:600])
except Exception as e:
    print("err:", e)
