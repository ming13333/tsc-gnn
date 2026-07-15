#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Targeted, verification-gated re-lookup for pending/flagged DOI refs (v2).

Strategy (robust under flaky proxy):
  T1b  Journal article with volume + first page:
       query Crossref by TITLE, then verify container-title (journal) +
       volume + first-page + year + title-similarity. High confidence.
  T2   Fallback (no vol / preprint / conference):
       title-scan; accept only if title sim>=0.92 (high) or 0.88-0.92
       (medium) AND venue/year plausible; else pending.

Every accepted DOI prints full evidence so it can be eyeballed against the
manuscript's ground-truth authors/title/venue/year.
"""
import re, json, time, os, urllib.request, urllib.parse, difflib, sys
sys.path.insert(0, r"C:/D 盘\科研\虚拟敲除\gate1\rewiring_study")
import doi_lookup as D

HERE = r"C:/D 盘\科研\虚拟敲除\gate1\rewiring_study"
RES  = HERE + r"\doi_results.json"
MAIL = "ming.luo@med.org"

def xfetch2(url, timeout=12, retries=5):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "TSC-GNN-DOI/2.0 (mailto:%s)" % MAIL})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.load(r)["message"]["items"]
        except Exception as e:
            last = e
            time.sleep(1.5 + i)
    return []

def norm(s):
    s = (s or "").lower()
    s = re.sub(r"[^a-z0-9 ]", " ", s)
    return re.sub(r"\s+", " ", s).strip()

def first_page(pages):
    return re.split(r"[–\-]", pages)[0].strip()

def t1b(ref, verbose=True):
    """Title query + verify journal/vol/page/year/title."""
    full = D.JMAP.get(ref["journal"])
    if not full or not ref.get("vol"):
        return None
    fp = first_page(ref["pages"]); y = ref["year"]; vol = ref["vol"]
    nt = norm(ref["title"])
    q = ("query.bibliographic=" + urllib.parse.quote(ref["title"]) +
         "&rows=30&mailto=" + MAIL)
    items = xfetch2("https://api.crossref.org/works?" + q)
    if not items and verbose:
        print("    [T1b] empty (proxy reset?)")
    for it in items:
        ct = (it.get("container-title") or [""])[0]
        if norm(full) not in norm(ct):
            continue
        if str(it.get("volume") or "") != vol:
            continue
        ipage = str(it.get("page") or "")
        if not (ipage.startswith(fp) or fp in ipage):
            continue
        dp = (it.get("issued") or {}).get("date-parts")
        cy = dp[0][0] if (dp and dp[0] and dp[0][0]) else None
        if cy is not None and abs(int(cy) - y) > 1:
            continue
        mt = (it.get("title") or [""])[0]
        ts = difflib.SequenceMatcher(None, nt, norm(mt)).ratio()
        if ts >= 0.6:
            au = ", ".join((a.get("family") or "") for a in it.get("author", [])[:3])
            if verbose:
                print(f"    [T1b] HIT doi={it.get('DOI')} sim={ts:.2f} '{mt}' | {au} {ct} {vol}:{ipage} yr={cy}")
            return {"doi": it.get("DOI"), "method": "T1b(title+vol+page+yr)",
                    "conf": "high", "match_title": mt, "cy": cy, "tsim": round(ts, 2)}
    return None

def t2(ref, verbose=True):
    items = xfetch2("https://api.crossref.org/works?query.title=" +
                    urllib.parse.quote(ref["title"]) + "&rows=30&mailto=" + MAIL)
    if not items:
        if verbose: print("    [T2] empty (proxy reset?)")
        return None
    nt = norm(ref["title"]); best=None; bs=-1
    for it in items:
        sc = difflib.SequenceMatcher(None, nt, norm((it.get("title") or [""])[0])).ratio()
        if sc > bs: bs=sc; best=it
    mt = (best.get("title") or [""])[0]
    dp = (best.get("issued") or {}).get("date-parts")
    cy = dp[0][0] if (dp and dp[0] and dp[0][0]) else None
    sc = difflib.SequenceMatcher(None, nt, norm(mt)).ratio()
    conf = "high" if sc >= 0.95 else ("medium" if sc >= 0.90 else "low")
    au = ", ".join((a.get("family") or "") for a in best.get("author", [])[:3])
    ct = (best.get("container-title") or [""])[0]
    if verbose:
        print(f"    [T2] best doi={best.get('DOI')} sim={sc:.3f} conf={conf} '{mt}' | {au} ({ct}) yr={cy}")
    return {"doi": best.get("DOI"), "method": "T2(title-scan)", "conf": conf,
            "match_title": mt, "cy": cy, "sim": round(sc, 3)}

def main():
    text = open(D.SRC, encoding="utf-8").read()
    lines = text.splitlines()
    start = next(i for i, l in enumerate(lines) if l.strip() == "## References")
    raws = {}
    for l in lines[start + 1:]:
        if l.startswith("## ") and l.strip() != "## References":
            break
        if re.match(r"^\d+\.", l.strip()):
            raw = re.sub(r"\s*\[DOI:[^\]]*\]\s*$", "", l.strip())
            num = int(raw.split(".", 1)[0]); raws[num] = raw
    data = json.load(open(RES, encoding="utf-8"))
    bynum = {r["num"]: r for r in data}

    # journal refs first (T1b), then no-vol/conference/preprint (T2)
    JOURNAL = [1,2,3,4,10,11,19,25,28,29,30,37,39,41,42,44,47,48,49,50,54]
    OTHER   = [15,17,38,40]

    def process(num, use_t2_only=False):
        raw = raws.get(num)
        if not raw:
            print(f"[{num:2d}] NOT FOUND in manuscript"); return
        ref = D.parse_ref(raw)
        print(f"[{num:2d}] {ref.get('kind'):11s} {ref.get('title','')[:68]}")
        res = None
        if ref["kind"] == "journal" and not use_t2_only:
            res = t1b(ref)
            if not res:
                res = t2(ref)
        else:
            res = t2(ref)
        if res and res.get("doi") and res["conf"] in ("high", "medium"):
            res.update(num=num, title=ref.get("title",""), year=ref.get("year"),
                       kind=ref.get("kind"), journal=ref.get("journal"))
            bynum[num] = res
            print(f"    => ACCEPT {res['doi']} ({res['conf']}, {res['method']})")
        else:
            if num in bynum and bynum[num].get("doi") and bynum[num].get("conf") in ("high","medium"):
                print(f"    => KEEP existing {bynum[num]['doi']}")
            else:
                bynum[num] = {"num": num, "doi": None, "conf": "low",
                              "method": (res or {}).get("method") or "none",
                              "title": ref.get("title",""), "year": ref.get("year"),
                              "kind": ref.get("kind"), "journal": ref.get("journal")}
                print("    => PENDING (low/no match)")
        json.dump(list(bynum.values()), open(RES, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        time.sleep(0.8)

    for n in JOURNAL: process(n)
    for n in OTHER:   process(n, use_t2_only=True)

    print("\n=== DONE ===")
    resolved = [k for k,v in bynum.items() if v.get("doi") and v.get("conf") in ("high","medium")]
    pending  = [k for k,v in bynum.items() if not (v.get("doi") and v.get("conf") in ("high","medium"))]
    print("resolved:", sorted(resolved))
    print("pending :", sorted(pending))

if __name__ == "__main__":
    main()
