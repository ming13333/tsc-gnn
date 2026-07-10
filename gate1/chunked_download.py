"""
chunked_download.py — 分块 HTTP Range 下载，彻底解决 NCBI 限速截断。

原理：
  - 先 HEAD 取真实 Content-Length（如 36MB，远超预期！）
  - 切成 CHUNK(4MB) 小块，逐块 `curl --range start-end` 下载
  - 每块校验实际字节数 == 请求字节数，失败重试
  - 全部块到位后按序拼接成完整文件
  - gzip -t 终检，通过才算 OK

每块只有 4MB，远低于 NCBI 单次截断阈值(~8MB)，故几乎不会被切断；
即便某块失败也只重下那一块，不污染其他块。
"""
import os
import re
import subprocess
import concurrent.futures as cf
from gate1 import data_acquisition as da

GEO_HTTPS = "https://ftp.ncbi.nlm.nih.gov"
CHUNK = 400 * 1024  # 400KB 每块（NCBI 当前截断~1MB，留足余量）


def sample_suppl_url(gsm: str) -> str:
    prefix = gsm[: len(gsm) - 3]
    return f"{GEO_HTTPS}/geo/samples/{prefix}nnn/{gsm}/suppl/"


def list_sample_files(gsm: str):
    import urllib.request
    url = sample_suppl_url(gsm)
    html = urllib.request.urlopen(url, timeout=30).read().decode("utf-8", "replace")
    return [f for f in re.findall(r'href="([^"]+)"', html)
            if f.endswith(".gz") and "disclos" not in f]


def get_size(url: str) -> int:
    r = subprocess.run(["curl", "-sI", "--max-time", "60", url],
                       capture_output=True, text=True)
    for line in r.stdout.splitlines():
        if line.lower().startswith("content-length:"):
            return int(line.split(":", 1)[1].strip())
    return -1


def download_chunk(url: str, dest: str, start: int, end: int, max_retry: int = 12):
    want = end - start + 1
    # 断点续传：已存在且大小正确的块直接跳过
    if os.path.exists(dest) and os.path.getsize(dest) == want:
        return True
    for _ in range(max_retry):
        subprocess.run(
            ["curl", "-sS", "--retry", "5", "--retry-delay", "2", "--max-time", "90",
             "-r", f"{start}-{end}", "-o", dest, url],
            capture_output=True)
        if os.path.exists(dest) and os.path.getsize(dest) == want:
            return True
        if os.path.exists(dest):
            os.remove(dest)
    return False


def gz_ok(path: str) -> bool:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    r = subprocess.run(["gzip", "-t", path], capture_output=True)
    return r.returncode == 0


def fetch_file(url: str, dest: str, label: str):
    total = get_size(url)
    if total <= 0:
        return dest, False, f"{label} HEAD_FAIL"
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    # 若已是完整 gzip，跳过
    if gz_ok(dest) and os.path.getsize(dest) == total:
        return dest, True, f"{label} CACHED"
    parts = []
    ok_all = True
    import time
    for i in range(0, total, CHUNK):
        end = min(i + CHUNK - 1, total - 1)
        part = dest + f".part{i//CHUNK:04d}"
        if not download_chunk(url, part, i, end):
            ok_all = False
            break
        parts.append(part)
        time.sleep(1.5)  # 温柔间隔，避免触发 NCBI 更强制裁
    if not ok_all:
        for p in parts:
            if os.path.exists(p):
                os.remove(p)
        return dest, False, f"{label} CHUNK_FAIL"
    # 拼接
    with open(dest, "wb") as out:
        for p in parts:
            with open(p, "rb") as inp:
                out.write(inp.read())
            os.remove(p)
    if gz_ok(dest) and os.path.getsize(dest) == total:
        return dest, True, f"{label} OK({total}B)"
    return dest, False, f"{label} GZ_FAIL({os.path.getsize(dest) if os.path.exists(dest) else 0}/{total})"


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
            continue
        sd = os.path.join(dest_root, "GSE174574", gsm)
        for suffix in ("barcodes.tsv.gz", "genes.tsv.gz", "matrix.mtx.gz"):
            fn = base + "_" + suffix
            url = sample_suppl_url(gsm) + fn
            dest = os.path.join(sd, fn)
            jobs.append((gsm, s["condition"], dest, url, fn))

    print(f"[plan] {len(jobs)} files, chunked Range download (CHUNK={CHUNK//1024//1024}MB)")
    results = []
    with cf.ThreadPoolExecutor(max_workers=3) as ex:
        futs = [ex.submit(fetch_file, url, dest, fn) for (_, _, dest, url, fn) in jobs]
        for (gsm, cond, dest, url, fn), fut in zip(jobs, futs):
            d, ok, msg = fut.result()
            print(f"  {'OK ' if ok else 'FAIL'} {gsm}({cond}) {fn} {msg}")
            results.append(ok)
    n_ok = sum(results)
    print(f"[done] {n_ok}/{len(jobs)} files OK")
    return n_ok == len(jobs)


if __name__ == "__main__":
    import sys
    ok = main(sys.argv[1] if len(sys.argv) > 1 else "data")
    print("ALL_OK" if ok else "PARTIAL")
