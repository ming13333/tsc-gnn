#!/usr/bin/env python
# download_gse197731.py
# ---------------------------------------------------------------------------
# 目的: 下载 GSE197731 (Kim et al. 2022, Redox Biol; C57BL/6J tMCAO 脑 scRNA-seq,
#       24h 与 48h 双时间点) 作为 CellChat 的 *第二个 24h 队列* 来源, 补充
#       GSE174574 在 24h 的稀疏 (仅 3 条 LR 记录), 使 24h 成为跨队列验证节点。
#
# ⚠️ 数据设计要点 (实测确认):
#   该数据集 *没有 sham 对照*。8 个样本为
#     24h/48h × WT/Prdx1-KO × 同侧(Ipsil, 缺血半球)/对侧(Cont, 非缺血半球)
#   因此 CellChat 的 "对照" 用 *对侧半球(Cont)* 代替 sham, 且只取 WT 亚组避免
#   基因型污染: 24h_WT_Ipsil vs 24h_WT_Cont 即 24h 的 "缺血 vs 非缺血" 重布线,
#   与 GSE174574 的 24h MCAO-vs-sham 概念平行, 可做跨队列复现。
#   详见 cellchat_gse197731.py。
#
# 用法:
#   python download_gse197731.py                 # 默认下载到本脚本所在目录
#   python download_gse197731.py --out /path/out # 指定输出目录
#
# 产出:
#   <out>/GSE197731_RAW.tar           原始 supplementary 包
#   <out>/GSE197731_RAW/...          解压后的各样本文件
#   <out>/gse197731_meta.json        样本->(time, genotype, side, condition) 映射
# ---------------------------------------------------------------------------
import os, sys, json, argparse, urllib.request, tarfile

GSE = "GSE197731"
# GEO 文件分桶规则: 去掉 GSE 号末 3 位再补 nnn (GSE197731 -> GSE197nnn)
BUCKET = "GSE" + GSE[3:-3] + "nnn"  # "GSE197nnn"
CANDIDATE_BASES = [
    f"https://ftp.ncbi.nlm.nih.gov/geo/series/{BUCKET}/{GSE}/",
    f"ftp://ftp.ncbi.nlm.nih.gov/geo/series/{BUCKET}/{GSE}/",
]
UA = {"User-Agent": "Mozilla/5.0 (compatible; GSE-downloader/1.0)"}


def fetch(url, timeout=600):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def find_base():
    """探测哪个基址可访问 RAW.tar, 返回基址或 None。"""
    for b in CANDIDATE_BASES:
        try:
            url = f"{b}suppl/{GSE}_RAW.tar"
            data = fetch(url, timeout=60)
            if data and len(data) > 1000:
                return b
        except Exception as e:
            print(f"  [skip] {b}: {e}")
    return None


def parse_matrix_title(text):
    """从 GEO series matrix 解析 !Sample_title, 返回 {gsm: title}。"""
    titles = {}
    samples = []
    for ln in text.split("\n"):
        if not ln.startswith("!"):
            continue
        parts = ln.rstrip("\n").split("\t")
        if parts[0] == "!Sample_geo_accession":
            samples = [c.strip().strip('"') for c in parts[1:] if c.strip()]
        elif parts[0] == "!Sample_title":
            for gsm, t in zip(samples, parts[1:]):
                titles[gsm] = t.strip().strip('"')
    return titles


def infer_meta(title):
    """'24h_WT_Cont' -> {time:'24h', genotype:'WT', side:'Cont'}。"""
    t = title.lower()
    time = "24h" if "24h" in t else ("48h" if "48h" in t else "unknown")
    geno = "WT" if "_wt_" in t or t.startswith("24h_wt") or t.startswith("48h_wt") else ("KO" if "ko" in t else "unknown")
    side = "Ipsil" if "ipsil" in t else ("Cont" if "cont" in t else "unknown")
    # condition: 缺血侧(Ipsil)=treatment, 对侧(Cont)=control
    condition = "MCAO_ischemic" if side == "Ipsil" else ("control" if side == "Cont" else "unknown")
    return {"time": time, "genotype": geno, "side": side, "condition": condition, "title": title}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)))
    args = ap.parse_args()
    out = args.out
    os.makedirs(out, exist_ok=True)
    raw_tar = os.path.join(out, f"{GSE}_RAW.tar")

    print(f"[1/4] 探测 NCBI 基址 ...")
    base = find_base()
    if base is None:
        print("\n❌ 无法访问 NCBI GEO。请确认运行环境能访问 https://ftp.ncbi.nlm.nih.gov")
        sys.exit(2)

    print(f"[2/4] 下载 RAW.tar ({base}) ...")
    if os.path.exists(raw_tar) and os.path.getsize(raw_tar) > 1_000_000:
        print(f"   已存在 ({os.path.getsize(raw_tar)//1024//1024} MB), 跳过下载")
    else:
        data = fetch(f"{base}suppl/{GSE}_RAW.tar")
        with open(raw_tar, "wb") as f:
            f.write(data)
        print(f"   下载完成 ({len(data)//1024//1024} MB)")

    print(f"[3/4] 解压 RAW.tar ...")
    extract_dir = os.path.join(out, f"{GSE}_RAW")
    os.makedirs(extract_dir, exist_ok=True)
    with tarfile.open(raw_tar) as tf:
        tf.extractall(extract_dir)
    print(f"   解压到 {extract_dir}")

    print(f"[4/4] 解析样本标签 -> gse197731_meta.json ...")
    # 优先用本地已下的 matrix; 否则从 base 下载
    matrix_path = os.path.join(out, f"{GSE}_series_matrix.txt.gz")
    if not os.path.exists(matrix_path):
        mdata = fetch(f"{base}matrix/{GSE}_series_matrix.txt.gz")
        with open(matrix_path, "wb") as f:
            f.write(mdata)
    import gzip
    with gzip.open(matrix_path, "rt") as f:
        matrix_text = f.read()
    titles = parse_matrix_title(matrix_text)
    meta = {gsm: infer_meta(t) for gsm, t in titles.items()}
    with open(os.path.join(out, "gse197731_meta.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(f"   解析到 {len(meta)} 个样本:")
    for gsm, m in meta.items():
        print(f"     {gsm}: {m['title']}  -> time={m['time']} geno={m['genotype']} side={m['side']}")
    print(f"\n完成。输出目录: {out}")
    print("下一步: 运行 cellchat_gse197731.py 做 CellChat (Ipsil-vs-Cont, WT only) + 跨队列 24h 一致性。")


if __name__ == "__main__":
    main()
