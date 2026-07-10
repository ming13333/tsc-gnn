"""
robust_download.py — 抗 NCBI 限流的健壮下载器。

NCBI 限流时的诡异行为：
  - 有时返回 206 + 严格 Range（正常）
  - 有时忽略 Range，返回 200 + 从 0 开始的截断数据（截断点 0.4~2MB 不等）
之前的分块校验(size==want)在 200 情况下永远失败 → 死循环。

本方案：每次请求"剩余全量" (-r pos-(total-1))，根据响应码决定如何拼接：
  - 206: 响应体 = [pos, pos+R)，整段追加，pos += R
  - 200: 响应体 = [0, R)（从0截断），只追加 [pos, R) 这段新尾部，pos = R
无论哪种，都只追加"之前没有的字节"，绝重叠、绝漏字节。
循环直到 pos==total，最后 gzip -t 校验。
单文件串行、块间留间隔，尽量不触发更强制裁。
"""
import os
import re
import time
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


def get_size(url: str) -> int:
    for _ in range(5):
        r = subprocess.run(["curl", "-sI", "--max-time", "60", url],
                           capture_output=True, text=True)
        for line in r.stdout.splitlines():
            if line.lower().startswith("content-length:"):
                try:
                    return int(line.split(":", 1)[1].strip())
                except ValueError:
                    pass
        time.sleep(2)
    return -1


def gz_ok(path: str) -> bool:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    return subprocess.run(["gzip", "-t", path], capture_output=True).returncode == 0


def fetch_one(url: str, dest: str, label: str, max_retry: int = 60):
    total = get_size(url)
    if total <= 0:
        return dest, False, f"{label} HEAD_FAIL"
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if gz_ok(dest) and os.path.getsize(dest) == total:
        return dest, True, f"{label} CACHED"
    # 断点：若 dest 已存在且是合法前缀则续传
    pos = 0
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        pos = os.path.getsize(dest)
    tmp = dest + ".tmp"
    attempt = 0
    last_pos = -1
    while pos < total and attempt < max_retry:
        attempt += 1
        # 请求剩余全量
        hdr = dest + ".hdr"
        subprocess.run(
            ["curl", "-sS", "--retry", "3", "--retry-delay", "2", "--max-time", "120",
             "-D", hdr, "-r", f"{pos}-{total-1}", "-o", tmp, url],
            capture_output=True)
        if not os.path.exists(tmp) or os.path.getsize(tmp) == 0:
            time.sleep(3)
            continue
        # 解析响应码
        code = 0
        if os.path.exists(hdr):
            with open(hdr) as f:
                first = f.readline()
                if first.startswith("HTTP"):
                    try:
                        code = int(first.split()[1])
                    except (ValueError, IndexError):
                        code = 0
        rlen = os.path.getsize(tmp)
        if code == 206:
            # 响应体 = [pos, pos+rlen)，整段追加
            with open(dest, "ab") as out, open(tmp, "rb") as inp:
                out.write(inp.read())
            pos += rlen
        else:
            # 200（忽略Range，从0截断）：响应体 = [0, rlen)，只追加 [pos, rlen)
            with open(tmp, "rb") as inp:
                data = inp.read()
            if rlen > pos:
                with open(dest, "ab") as out:
                    out.write(data[pos:])
                pos = rlen
            # 若 rlen<=pos 说明无新数据，保持 pos，等退避后重试
        os.remove(tmp)
        if os.path.exists(hdr):
            os.remove(hdr)
        if pos == last_pos:
            time.sleep(5)  # 无进展，退避
        last_pos = pos
        time.sleep(1.0)  # 温柔间隔
    if pos >= total and gz_ok(dest):
        return dest, True, f"{label} OK({total}B,{attempt}次)"
    return dest, False, f"{label} FAIL(pos={pos}/{total},{attempt}次)"


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
    print(f"[plan] {len(jobs)} files, robust prefix-accumulating download")
    results = []
    with cf.ThreadPoolExecutor(max_workers=2) as ex:
        futs = [ex.submit(fetch_one, url, dest, fn) for (_, _, dest, url, fn) in jobs]
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
