#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Probe3: fetch SOX10 CRISPR signature gene sets (up/down entities) from SigCom."""
import requests, json

META = "https://maayanlab.cloud/sigcom-lincs/metadata-api"
LIB_ID = "96c7b8c5-1eca-5764-88e4-e4ccaee6603f"  # l1000_xpr CRISPR

# find SOX10 signatures in this library
r = requests.post(META + "/signatures/find",
                  json={"filter": {"where": {"meta.pert_name": "SOX10", "library": LIB_ID},
                                   "limit": 20}}, timeout=60)
sigs = r.json()
print("SOX10 sigs in CRISPR lib:", len(sigs))
for s in sigs[:5]:
    print("\n--- uuid:", s.get("id"), "---")
    print("meta keys:", sorted(s.get("meta", {}).keys()))
    m = s.get("meta", {})
    print("pert_name:", m.get("pert_name"), "| cell:", m.get("cell", m.get("cell_id")),
          "| pert_type:", m.get("pert_type"), "| dose:", m.get("pert_dose"))
    up = m.get("up_entities") or m.get("up") or []
    dn = m.get("down_entities") or m.get("down") or []
    print(f"up_entities len={len(up)} sample={up[:5]}; down_entities len={len(dn)} sample={dn[:5]}")
    # also check if genes stored as symbols
    ups = m.get("up_genes") or m.get("up_symbols") or []
    dns = m.get("down_genes") or m.get("down_symbols") or []
    print(f"up_genes(sym) len={len(ups)} sample={ups[:5]}; down_genes len={len(dns)} sample={dns[:5]}")
