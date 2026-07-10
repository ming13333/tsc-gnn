"""
data_acquisition.py — 三个核心数据集的下载与元数据抓取（HTTPS 版）。

注意：沙箱内 FTP(21端口)被阻断，但 https://ftp.ncbi.nlm.nih.gov 可达，故全部走 HTTPS。

数据集（详见研究步骤与条件_2026-07-08.md）：
  - GSE174574  (Li 2021, J Cereb Blood Flow Metab) : 3 MCAO(24h)+3 sham 鼠 scRNA，单时间点(24h)。
                 作为 24h 锚点 + 含 BEAM 伪时间轨迹。本探索中用作"sham(基线) vs MCAO(24h 损伤)"
                 两状态轴的可行性验证。
  - GSE225948  (Anrather 2024, Nature Immunology)  : 真实 2d + 14d 术后 scRNA 时间进程(脑+血)。
                 提供 2d/14d 两个真实时间点，用于与 GSE174574 整合出 24h→2d→14d 真实时间轴。
  - HRA007397  (NGDC 中国库)                       : 3 例人 AIS 取栓前 + 术后 day1/day7 的 scRNA。
                 自带 1d/7d 人时间轴 → 跨物种独立验证队列（需 NGDC 注册/申请，本机手动下）。

下载函数直接落盘到 data/<gse>/<gsm>/，由 preprocessing.py 读取。
"""
import os
import re
import gzip
import urllib.request

GEO_HTTPS = "https://ftp.ncbi.nlm.nih.gov"


# ----------------------------------------------------------------------------
# 底层：样本级 suppl 目录 listing + 下载
# ----------------------------------------------------------------------------
def _sample_suppl_url(gse_id: str, gsm: str) -> str:
    prefix = gsm[: len(gsm) - 3]  # 'GSM5319987' -> 'GSM5319'
    return f"{GEO_HTTPS}/geo/samples/{prefix}nnn/{gsm}/suppl/"


def list_sample_files(gse_id: str, gsm: str):
    """返回该样本 suppl 目录下所有 .gz 文件名列表。"""
    url = _sample_suppl_url(gse_id, gsm)
    html = urllib.request.urlopen(url, timeout=30).read().decode("utf-8", "replace")
    return [f for f in re.findall(r'href="([^"]+)"', html)
            if f.endswith(".gz") and "disclos" not in f]


def _download_file(url: str, dest: str, verbose: bool = True):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        if verbose:
            print(f"  [skip] {os.path.basename(dest)}")
        return dest
    if verbose:
        print(f"  [get ] {url.split('/')[-1]}")
    urllib.request.urlretrieve(url, dest)
    return dest


def download_sample_mtx(gse_id: str, gsm: str, dest_dir: str, verbose: bool = True):
    """下载某 GSM 样本的 10x mtx 三件套到 dest_dir，返回含 barcodes/genes/matrix 的 dict。"""
    os.makedirs(dest_dir, exist_ok=True)
    files = list_sample_files(gse_id, gsm)
    # 找 matrix 文件名前缀（如 GSM5319987_sham1）
    base = None
    for f in files:
        if f.endswith("_matrix.mtx.gz"):
            base = f[: -len("_matrix.mtx.gz")]
            break
    if base is None:
        if verbose:
            print(f"  [warn] no matrix.mtx.gz in {gsm}: {files}")
        return None
    paths = {}
    for suffix in ("barcodes.tsv.gz", "genes.tsv.gz", "matrix.mtx.gz"):
        fn = base + "_" + suffix
        url = _sample_suppl_url(gse_id, gsm) + fn
        dest = os.path.join(dest_dir, fn)
        _download_file(url, dest, verbose)
        paths[suffix.split(".")[0]] = dest
    return paths


# ----------------------------------------------------------------------------
# 解析 series matrix：样本 -> GSM / 条件(sham|MCAO) / 时间标签
# ----------------------------------------------------------------------------
def parse_series_matrix(gse_id: str):
    """下载并解析 series matrix，返回样本清单 list[dict(gsm, title, condition, time_label)]。"""
    url = (f"{GEO_HTTPS}/geo/series/{gse_id[:len(gse_id)-3]}nnn/{gse_id}"
           f"/matrix/{gse_id}_series_matrix.txt.gz")
    data = urllib.request.urlopen(url, timeout=60).read()
    text = gzip.decompress(data).decode("utf-8", "replace")
    lines = text.splitlines()

    # 收集分节
    geo_acc, titles, chars = [], [], []
    for ln in lines:
        if ln.startswith("!Sample_geo_accession"):
            geo_acc = [x.strip().strip('"').strip("'") for x in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_title"):
            titles = [x.strip().strip('"').strip("'") for x in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_characteristics_ch1"):
            chars = [x.strip().strip('"').strip("'") for x in ln.split("\t")[1:]]
        if ln.startswith("!series_matrix_table_begin"):
            break

    samples = []
    for i, gsm in enumerate(geo_acc):
        title = titles[i] if i < len(titles) else ""
        char = chars[i] if i < len(chars) else ""
        low = (title + " " + char).lower()
        if "mcko" in low or "mca" in low or "stroke" in low:
            condition = "MCAO"
        elif "sham" in low:
            condition = "sham"
        else:
            condition = "unknown"
        # 时间标签：单点 24h
        time_label = "24h" if condition == "MCAO" else "sham"
        samples.append({"gsm": gsm, "title": title, "condition": condition,
                        "time_label": time_label})
    return samples


# ----------------------------------------------------------------------------
# 高层：下载整个 GSE
# ----------------------------------------------------------------------------
def download_gse(gse_id: str, dest_root: str, verbose: bool = True):
    """下载 GSE 所有样本的 mtx 三件套，返回 (data_dir, samples_manifest)。"""
    data_dir = os.path.join(dest_root, gse_id)
    samples = parse_series_matrix(gse_id)
    if verbose:
        print(f"[GSE {gse_id}] {len(samples)} samples: "
              + ", ".join(f"{s['gsm']}({s['condition']})" for s in samples))
    manifest = []
    for s in samples:
        sd = os.path.join(data_dir, s["gsm"])
        paths = download_sample_mtx(gse_id, s["gsm"], sd, verbose)
        if paths:
            manifest.append({**s, "mtx": paths})
    return data_dir, manifest


def download_gse174574(dest_root: str):
    return download_gse("GSE174574", dest_root)


def download_gse225948(dest_root: str):
    return download_gse("GSE225948", dest_root)


def hra007397_info():
    """HRA007397 (NGDC/GSA-Human) 获取说明（需注册/申请，本机手动下）。"""
    return {
        "accession": "HRA007397",
        "portal": "https://ngdc.cncb.ac.cn/bioproject/",
        "note": ("人 AIS 取栓前 + 术后 day1/day7 的 scRNA。请于 NGDC/GSA-Human 检索 HRA007397，"
                 "按站点指引下载（可能需要注册账号并申请数据使用）。下载后放置于 data/ 目录，"
                 "由 preprocessing.py 读取。"),
        "expected_time_axis": ["pre", "day1", "day7"],
    }


if __name__ == "__main__":
    import sys
    dr = sys.argv[1] if len(sys.argv) > 1 else "data"
    d, m = download_gse174574(dr)
    print("saved:", d, "samples:", len(m))
