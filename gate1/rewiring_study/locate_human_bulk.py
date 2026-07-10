#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
定位并下载公开人卒中 bulk RNA-seq（stroke vs control）
探测 GSE16561 / GSE22255 / GSE37409 等候选，解析样本特征确认分组。
"""
import urllib.request, gzip, os, sys, io

HERE = os.path.dirname(os.path.abspath(__file__))
CANDIDATES = {
    "GSE16561": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE16nnn/GSE16561/matrix/GSE16561_series_matrix.txt.gz",
    "GSE22255": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22255/matrix/GSE22255_series_matrix.txt.gz",
    "GSE37409": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE37nnn/GSE37409/matrix/GSE37409_series_matrix.txt.gz",
}

def fetch(url, timeout=120):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    return urllib.request.urlopen(req, timeout=timeout).read()

def parse_matrix(raw_bytes):
    txt = gzip.decompress(raw_bytes).decode("utf-8", "ignore")
    lines = txt.splitlines()
    meta = {}
    table_start = None
    for i, l in enumerate(lines):
        if l.startswith("!series_matrix_table_begin"):
            table_start = i
            break
        if l.startswith("!"):
            key = l[1:].split("\t")[0]
            meta[key] = l.split("\t")[1:]
    return meta, table_start, lines

def summarize(gse, url):
    print("="*70)
    print("TRY", gse, url)
    try:
        raw = fetch(url)
        print("  downloaded bytes =", len(raw))
    except Exception as e:
        print("  DOWNLOAD FAIL:", type(e).__name__, str(e)[:160])
        return None
    meta, t0, lines = parse_matrix(raw)
    n = len(meta.get("!Sample_geo_accession", []))
    print("  n_samples =", n)
    acc = meta.get("!Sample_geo_accession", [])
    chars = meta.get("!Sample_characteristics_ch1", [])
    title = meta.get("!Sample_title", [])
    tissue = meta.get("!Sample_tissue_ch1", [])
    # print all sample characteristics to decide grouping
    for i in range(n):
        c = chars[i] if i < len(chars) else ""
        t = title[i] if i < len(title) else ""
        ti = tissue[i] if i < len(tissue) else ""
        print(f"  {acc[i]:10s} | {ti:14s} | {t[:30]:30s} | {c[:90]}")
    # save parsed matrix to npz for later use
    if t0 is not None:
        import numpy as np, re
        # find header row = line after table_begin
        hdr = lines[t0+1].split("\t")
        genes = [h for h in hdr[1:]]
        mat = []
        for l in lines[t0+2:]:
            if l.startswith("!series_matrix_table_end"):
                break
            parts = l.split("\t")
            mat.append([float(x) if x not in ("", "NA") else float("nan") for x in parts[1:]])
        M = np.array(mat)  # genes x samples
        out = os.path.join(HERE, f"{gse}_bulk.npz")
        np.savez(out, expr=M, genes=np.array(genes, dtype=object),
                 acc=np.array(acc, dtype=object),
                 chars=np.array(chars, dtype=object),
                 title=np.array(title, dtype=object),
                 tissue=np.array(tissue, dtype=object))
        print("  SAVED", out, "shape", M.shape)
    return gse

if __name__ == "__main__":
    # try candidates in order until one downloads + parses
    for gse, url in CANDIDATES.items():
        r = summarize(gse, url)
        if r:
            print("\n>>> USABLE:", r)
            break
        else:
            print()
