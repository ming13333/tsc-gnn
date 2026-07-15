#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Print top Crossref candidates (full metadata) for each target ref so the
final DOI set can be curated by hand with full evidence. No auto-accept."""
import re, json, time, urllib.request, urllib.parse, difflib, sys
sys.path.insert(0, r"C:/D 盘\科研\虚拟敲除\gate1\rewiring_study")
import doi_lookup as D

HERE = r"C:/D 盘\科研\虚拟敲除\gate1\rewiring_study"
RES  = HERE + r"\doi_results.json"
MAIL = "ming.luo@med.org"

def xfetch(url, timeout=12, retries=5):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "TSC-GNN-DOI/ev (mailto:%s)" % MAIL})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.load(r)["message"]["items"]
        except Exception:
            time.sleep(1.5 + i)
    return []

def norm(s):
    s = (s or "").lower(); s = re.sub(r"[^a-z0-9 ]", " ", s)
    return re.sub(r"\s+", " ", s).strip()

def candidates(title, rows=15):
    items = xfetch("https://api.crossref.org/works?query.title=" +
                   urllib.parse.quote(title) + f"&rows={rows}&mailto=" + MAIL)
    nt = norm(title); out = []
    for it in items:
        mt = (it.get("title") or [""])[0]
        ct = (it.get("container-title") or [""])[0]
        au = ", ".join((a.get("family") or "") for a in it.get("author", [])[:4])
        dp = (it.get("issued") or {}).get("date-parts")
        cy = dp[0][0] if (dp and dp[0] and dp[0][0]) else None
        sim = difflib.SequenceMatcher(None, nt, norm(mt)).ratio()
        out.append({"doi": it.get("DOI"), "sim": round(sim, 3), "title": mt,
                    "authors": au, "journal": ct, "vol": it.get("volume"),
                    "page": it.get("page"), "year": cy})
    out.sort(key=lambda x: -x["sim"])
    return out

def main():
    text = open(D.SRC, encoding="utf-8").read()
    lines = text.splitlines()
    start = next(i for i, l in enumerate(lines) if l.strip() == "## References")
    raws = {}
    for l in lines[start + 1:]:
        if l.startswith("## ") and l.strip() != "## References": break
        if re.match(r"^\d+\.", l.strip()):
            raw = re.sub(r"\s*\[DOI:[^\]]*\]\s*$","", l.strip())
            raws[int(raw.split(".",1)[0])] = raw
    # targets: currently pending + the 4 flagged wrong/needs-check
    data = json.load(open(RES, encoding="utf-8"))
    bynum = {r["num"]: r for r in data}
    pending = [k for k,v in bynum.items() if not (v.get("doi") and v.get("conf") in ("high","medium"))]
    flagged = [16, 18, 41, 47]
    targets = sorted(set(pending) | set(flagged))
    print(f"# Evidence for {len(targets)} refs (pending + flagged)\n")
    for num in targets:
        raw = raws.get(num)
        ref = D.parse_ref(raw)
        print(f"===== [{num}] MANUSCRIPT: {ref.get('title','')[:80]}")
        print(f"      authors={ref.get('authors','')[:60]} | {ref.get('journal','')} "
              f"{ref.get('vol') or '-'}:{ref.get('pages') or '-'} ({ref.get('year')})")
        cands = candidates(ref["title"])
        if not cands:
            print("      (no Crossref items - proxy reset?)")
        for c in cands[:4]:
            print(f"   sim={c['sim']:.3f} doi={c['doi']}")
            print(f"      '{c['title'][:90]}'")
            print(f"      {c['authors'][:60]} | {c['journal']} {c['vol']}:{c['page']} ({c['year']})")
        print()
        time.sleep(0.8)
    print("# END")

if __name__ == "__main__":
    main()
