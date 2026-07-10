"""
clean_download.py — GSE174574 干净下载器（curl 整文件 + 并行 + gzip 校验）。
前置：用户已重新配置网络，NCBI 限速已解除，整文件下载可稳定完成。
用法：python clean_download.py <dest_root>
"""
import os
import sys
import gzip
import subprocess
import concurrent.futures as cf

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/samples"
# (gsm, base) — base 是文件名前缀（含样本标签）
SAMPLES = [
    ("GSM5319987", "GSM5319987_sham1"),
    ("GSM5319988", "GSM5319988_sham2"),
    ("GSM5319989", "GSM5319989_sham3"),
    ("GSM5319990", "GSM5319990_MCAO1"),
    ("GSM5319991", "GSM5319991_MCAO2"),
    ("GSM5319992", "GSM5319992_MCAO3"),
]
SUFFIXES = ("barcodes.tsv.gz", "genes.tsv.gz", "matrix.mtx.gz")


def gz_ok(path: str) -> bool:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    return subprocess.run(["gzip", "-t", path], capture_output=True).returncode == 0


def fetch(url: str, dest: str, max_retry: int = 6) -> bool:
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if gz_ok(dest):
        return True
    for attempt in range(1, max_retry + 1):
        r = subprocess.run(
            ["curl", "-sS", "--retry", "4", "--retry-delay", "3",
             "--max-time", "600", "-o", dest, url],
            capture_output=True)
        if r.returncode == 0 and gz_ok(dest):
            return True
        if os.path.exists(dest) and os.path.getsize(dest) == 0:
            os.remove(dest)
    return gz_ok(dest)


def main():
    dest_root = sys.argv[1] if len(sys.argv) > 1 else "data"
    data_dir = os.path.join(dest_root, "GSE174574")
    os.makedirs(data_dir, exist_ok=True)

    jobs = []
    for gsm, base in SAMPLES:
        prefix = gsm[: len(gsm) - 3]
        url_dir = f"{GEO}/{prefix}nnn/{gsm}/suppl/"
        sdir = os.path.join(data_dir, gsm)
        for suf in SUFFIXES:
            fn = f"{base}_{suf}"
            jobs.append((url_dir + fn, os.path.join(sdir, fn)))

    print(f"[clean_download] {len(jobs)} files to fetch -> {data_dir}")

    def worker(job):
        url, dest = job
        ok = fetch(url, dest)
        return dest, ok

    ok_n = 0
    with cf.ThreadPoolExecutor(max_workers=6) as ex:
        for dest, ok in ex.map(worker, jobs):
            status = "OK " if ok else "FAIL"
            if ok:
                ok_n += 1
            print(f"  [{status}] {os.path.basename(dest)} ({os.path.getsize(dest)}B)")

    print(f"\n=== 完成 {ok_n}/{len(jobs)} 文件 ===")
    if ok_n == len(jobs):
        print("ALL_GOOD")
    else:
        print("SOME_FAILED")
        sys.exit(1)


if __name__ == "__main__":
    main()
