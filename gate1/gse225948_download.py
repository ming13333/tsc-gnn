"""
gse225948_download.py — 下载 GSE225948 的 brain 样本(counts.csv.gz + metadata.csv.gz)。
格式不同于 GSE174574(10x mtx)，这里是 CSV 表达矩阵 + 细胞注释。
用法：python gse225948_download.py <dest_root>
"""
import os
import sys
import gzip
import subprocess
import concurrent.futures as cf

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/samples"
# GSE225948 全部 18 个 BRAIN 样本 GSM（从 series matrix 解析得到，避開 peripheral blood）
BRAIN_GSMS = [
    "GSM7060815", "GSM7060816", "GSM7060817", "GSM7060818",  # sham Day2 young
    "GSM7060819", "GSM7060820", "GSM7060821", "GSM7060822",  # Stroke Day2 young
    "GSM7060823", "GSM7060824", "GSM7060825", "GSM7060826",  # Stroke Day14 young
    "GSM7060827", "GSM7060828",                                # sham Day2 aged
    "GSM7060829", "GSM7060830", "GSM7060831",                  # Stroke Day2 aged
    "GSM7060832",                                              # Stroke Day14 aged
]


def gz_ok(path: str) -> bool:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    return subprocess.run(["gzip", "-t", path], capture_output=True).returncode == 0


def list_sample_files(gse_id, gsm):
    prefix = gsm[: len(gsm) - 3]
    url = f"{GEO}/{prefix}nnn/{gsm}/suppl/"
    import urllib.request, re
    html = urllib.request.urlopen(url, timeout=30).read().decode("utf-8", "replace")
    return [f for f in re.findall(r'href="([^"]+)"', html)
            if f.endswith(".gz") and "disclos" not in f]


def fetch(url, dest, max_retry=6):
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if gz_ok(dest):
        return True
    for _ in range(max_retry):
        r = subprocess.run(
            ["curl", "-sS", "--retry", "4", "--retry-delay", "3",
             "--max-time", "600", "-o", dest, url], capture_output=True)
        if r.returncode == 0 and gz_ok(dest):
            return True
        if os.path.exists(dest) and os.path.getsize(dest) == 0:
            os.remove(dest)
    return gz_ok(dest)


def main():
    dest_root = sys.argv[1] if len(sys.argv) > 1 else "data"
    data_dir = os.path.join(dest_root, "GSE225948")
    os.makedirs(data_dir, exist_ok=True)
    jobs = []
    for gsm in BRAIN_GSMS:
        prefix = gsm[: len(gsm) - 3]
        url_dir = f"{GEO}/{prefix}nnn/{gsm}/suppl/"
        sdir = os.path.join(data_dir, gsm)
        for fn in list_sample_files("GSE225948", gsm):
            jobs.append((url_dir + fn, os.path.join(sdir, fn)))
    print(f"[gse225948] {len(jobs)} files ({len(BRAIN_GSMS)} brain samples) -> {data_dir}")
    ok_n = 0
    with cf.ThreadPoolExecutor(max_workers=6) as ex:
        for dest, ok in ex.map(lambda j: (j[1], fetch(j[0], j[1])), jobs):
            if ok:
                ok_n += 1
            print(f"  [{'OK' if ok else 'FAIL'}] {os.path.basename(dest)} ({os.path.getsize(dest)}B)")
    print(f"\n=== 完成 {ok_n}/{len(jobs)} ===")
    if ok_n != len(jobs):
        sys.exit(1)


if __name__ == "__main__":
    main()
