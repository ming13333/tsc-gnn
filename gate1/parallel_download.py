"""
parallel_download.py v2 — 并行下载 + gzip 完整性门禁。
每个文件下载后用 gzip -t 校验，失败则整文件重下（不续传，避免截断拼接），
最多重试 8 次，直到全部通过或放弃。
"""
import os
import re
import subprocess
import concurrent.futures as cf
from gate1 import data_acquisition as da

GEO_HTTPS = "https://ftp.ncbi.nlm.nih.gov"


def sample_suppl_url(gsm: str) -> str:
    prefix = gsm[: len(gsm) - 3]
    return f"{GEO_HTTPS}/geo/samples/{prefix}nnn/{gsm}/suppl/"


def list_sample_files(gsm: str):
    import urllib.request
    url = sample_suppl_url(gsm)
    html = urllib.request.urlopen(url, timeout=30).read().decode("utf-8", "replace")
    return [f for f in re.findall(r'href="([^"]+)"', html)
            if f.endswith(".gz") and "disclos" not in f]


def gz_ok(path: str) -> bool:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    r = subprocess.run(["gzip", "-t", path], capture_output=True)
    return r.returncode == 0


def curl_fetch_strict(url: str, dest: str, max_retry: int = 15):
    """整文件删一次 + 循环续传(-C -) + gzip 校验，直到完整或放弃。
    关键：只删一次（避免破坏续传），之后每次 curl 从断点继续，
    直到 gzip -t 通过（整文件完整）才停。解决 NCBI 限速截断问题。
    """
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    # 启动前清掉可能损坏的残文件，保证从头/干净续传
    if os.path.exists(dest):
        os.remove(dest)
    for attempt in range(1, max_retry + 1):
        cmd = ["curl", "-sS", "-C", "-", "--retry", "3", "--retry-delay", "2",
               "--max-time", "600", "-o", dest, url]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode == 0 and gz_ok(dest):
            return dest, True, f"attempt{attempt}"
        # 否则保留断点，下次 -C - 续传
    return dest, False, "exhausted-retries"


def main(dest_root: str = "data"):
    samples = da.parse_series_matrix("GSE174574")
    print(f"[GSE174574] {len(samples)} samples")
    jobs = []
    for s in samples:
        gsm = s["gsm"]
        files = list_sample_files(gsm)
        base = None
        for f in files:
            if f.endswith("_matrix.mtx.gz"):
                base = f[: -len("_matrix.mtx.gz")]
                break
        if base is None:
            print(f"  [warn] {gsm}: no matrix in {files}")
            continue
        sd = os.path.join(dest_root, "GSE174574", gsm)
        for suffix in ("barcodes.tsv.gz", "genes.tsv.gz", "matrix.mtx.gz"):
            fn = base + "_" + suffix
            url = sample_suppl_url(gsm) + fn
            dest = os.path.join(sd, fn)
            jobs.append((gsm, s["condition"], dest, url))

    print(f"[plan] {len(jobs)} files, parallel fetch + gz verify")
    results = []
    with cf.ThreadPoolExecutor(max_workers=12) as ex:
        futs = [ex.submit(curl_fetch_strict, url, dest) for (_, _, dest, url) in jobs]
        for (gsm, cond, dest, url), fut in zip(jobs, futs):
            d, ok, msg = fut.result()
            sz = os.path.getsize(d) if os.path.exists(d) else 0
            print(f"  {'OK ' if ok else 'FAIL'} {gsm}({cond}) {os.path.basename(d)} {sz}B {msg}")
            results.append(ok)
    n_ok = sum(results)
    print(f"[done] {n_ok}/{len(jobs)} files OK")
    return n_ok == len(jobs)


if __name__ == "__main__":
    import sys
    ok = main(sys.argv[1] if len(sys.argv) > 1 else "data")
    print("ALL_OK" if ok else "PARTIAL")
