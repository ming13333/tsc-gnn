"""下载 1 个 GSE225948 brain 样本，探查 CSV 表达矩阵格式(行列/基因名)。"""
import os, sys, gzip, subprocess, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gate1 import data_acquisition as da

GSM = "GSM7060819"  # young male Stroke Day 2 brain
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7060nnn/GSM7060819/suppl/"
outdir = "data/GSE225948_probe"
os.makedirs(outdir, exist_ok=True)

files = da.list_sample_files("GSE225948", GSM)
print("样本文件:", files)
for f in files:
    url = GEO + f
    dest = os.path.join(outdir, f)
    if not os.path.exists(dest):
        subprocess.run(["curl", "-sS", "--max-time", "300", "-o", dest, url], check=True)
    print(f"下载: {f} ({os.path.getsize(dest)} B)")

# 探查 counts.csv.gz
cf = [f for f in files if "counts" in f][0]
mf = [f for f in files if "metadata" in f][0]
cpath = os.path.join(outdir, cf)
mpath = os.path.join(outdir, mf)

print("\n=== counts.csv.gz 探查 ===")
with gzip.open(cpath, "rt") as fh:
    r = csv.reader(fh)
    header = next(r)
    print("列数:", len(header), "前5列名:", header[:5])
    row0 = next(r)
    print("行1前5值:", row0[:5])
    # 数行数
    n = 1
    for _ in r:
        n += 1
    print("数据行数(基因数):", n)

print("\n=== metadata.csv.gz 探查 ===")
with gzip.open(mpath, "rt") as fh:
    r = csv.reader(fh)
    mh = next(r)
    print("metadata 列:", mh)
    for i, row in enumerate(r):
        if i < 3:
            print("  行", i, ":", row[:6])
        else:
            break
