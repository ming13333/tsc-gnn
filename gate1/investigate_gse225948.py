"""探查 GSE225948 结构：样本清单/条件/时间点/文件格式/单样本大小。"""
import os
import sys
import time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gate1 import data_acquisition as da
import urllib.request
import re
import gzip

GEO = "https://ftp.ncbi.nlm.nih.gov"

def head_size(url, timeout=30):
    try:
        req = urllib.request.Request(url, method="HEAD")
        r = urllib.request.urlopen(req, timeout=timeout)
        cl = r.headers.get("Content-Length")
        return int(cl) if cl else None
    except Exception as e:
        # 退化：用 Range 0-0 探测
        try:
            req = urllib.request.Request(url, headers={"Range": "bytes=0-0"})
            r = urllib.request.urlopen(req, timeout=timeout)
            cr = r.headers.get("Content-Range")
            if cr and "/" in cr:
                return int(cr.split("/")[-1])
        except Exception:
            pass
        return None

print("=== 解析 series matrix ===")
try:
    samples = da.parse_series_matrix("GSE225948")
except Exception as e:
    print("parse_series_matrix 失败:", e)
    samples = []

print(f"样本数: {len(samples)}")
for s in samples[:60]:
    print(f"  {s['gsm']:14s} cond={s['condition']:8s} time={s['time_label']:6s} title={s['title'][:50]}")

print("\n=== 抽查前 3 个样本 suppl 文件与大小 ===")
total = 0
checked = 0
for s in samples[:3]:
    gsm = s["gsm"]
    try:
        files = da.list_sample_files("GSE225948", gsm)
    except Exception as e:
        print(f"  {gsm}: list 失败 {e}")
        continue
    print(f"\n  {gsm} ({s['condition']}) 文件数={len(files)}")
    for f in files[:8]:
        url = da._sample_suppl_url("GSE225948", gsm) + f
        sz = head_size(url)
        if sz:
            total += sz
            checked += 1
        print(f"    {f:50s} {sz} B" if sz else f"    {f:50s} (size?)")
print(f"\n抽查 {checked} 文件合计 {total/1e6:.1f} MB（据此外推全集规模）")
