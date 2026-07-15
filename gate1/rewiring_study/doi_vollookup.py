#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Resolve journal refs by container-title + volume + first-page (most reliable).
Prints all page-matching candidates for eyeballing. No auto-accept into cache
unless page+volume+year all agree (high confidence); writes a JSON of decisions."""
import re, json, time, urllib.request, urllib.parse, sys
sys.path.insert(0, r"C:/D 盘\科研\虚拟敲除\gate1\rewiring_study")
import doi_lookup as D

HERE = r"C:/D 盘\科研\虚拟敲除\gate1\rewiring_study"
RES  = HERE + r"\doi_results.json"
MAIL = "ming.luo@med.org"

def xfetch(url, timeout=12, retries=5):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "TSC-GNN-DOI/vol (mailto:%s)" % MAIL})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.load(r)["message"]["items"]
        except Exception:
            time.sleep(1.5 + i)
    return []

def norm(s):
    s = (s or "").lower(); s = re.sub(r"[^a-z0-9 ]"," ",s); return re.sub(r"\s+"," ",s).strip()
def first_page(p): return re.split(r"[–\-]", p)[0].strip()

def vol_lookup(journal_abbr, vol, fp, year, verbose=True):
    full = D.JMAP.get(journal_abbr)
    if not full:
        return None, "no-jmap"
    q = (f"filter=container-title:{urllib.parse.quote(full)},volume:{vol}"
         f"&rows=60&mailto={MAIL}")
    items = xfetch("https://api.crossref.org/works?" + q)
    hits = []
    for it in items:
        ipage = str(it.get("page") or "")
        if ipage.startswith(fp) or fp in ipage:
            dp = (it.get("issued") or {}).get("date-parts")
            cy = dp[0][0] if (dp and dp[0] and dp[0][0]) else None
            hits.append({"doi": it.get("DOI"), "title": (it.get("title") or [""])[0],
                         "authors": ", ".join(a.get("family") or "" for a in it.get("author",[])[:3]),
                         "journal": (it.get("container-title") or [""])[0], "vol": it.get("volume"),
                         "page": ipage, "year": cy})
    return hits, full

def main():
    text = open(D.SRC, encoding="utf-8").read()
    lines = text.splitlines()
    start = next(i for i,l in enumerate(lines) if l.strip()=="## References")
    raws = {}
    for l in lines[start+1:]:
        if l.startswith("## ") and l.strip()!="## References": break
        if re.match(r"^\d+\.", l.strip()):
            raw = re.sub(r"\s*\[DOI:[^\]]*\]\s*$","",l.strip())
            raws[int(raw.split(".",1)[0])] = raw
    # journal refs with vol+page to resolve (exclude #2/#41 handled separately; #1 has page collision)
    NUMS = [3,4,19,25,28,29,30,37,39,44,48,49,50,54]
    data = json.load(open(RES, encoding="utf-8"))
    bynum = {r["num"]: r for r in data}
    decisions = {}
    for num in NUMS:
        ref = D.parse_ref(raws[num])
        fp = first_page(ref["pages"]); vol = ref["vol"]; y = ref["year"]
        print(f"===== [{num}] {ref['journal']} {vol}:{fp} ({y})  {ref['title'][:60]}")
        hits, full = vol_lookup(ref["journal"], vol, fp, y)
        if hits is None:
            print("      (no JMAP mapping)"); continue
        if not hits:
            print("      NO page-matching item found"); continue
        for h in hits[:6]:
            yok = (h["year"] is None) or (abs(int(h["year"])-y)<=1)
            print(f"   doi={h['doi']} yr={h['year']}{'*' if yok else '!!'}"
                  f" | {h['authors'][:30]} | {h['title'][:70]} | {h['journal']} {h['vol']}:{h['page']}")
        # accept: year within +-1 and title plausibly related (sim>=0.4) OR author matches first author
        best = None
        fa = ref["authors"].split(",")[0].strip().lower()
        for h in hits:
            yok = (h["year"] is None) or (abs(int(h["year"])-y)<=1)
            if not yok: continue
            ts = difflib.SequenceMatcher(None, norm(ref["title"]), norm(h["title"])).ratio()
            auth_ok = fa and fa in (h["authors"] or "").lower()
            if ts >= 0.4 or auth_ok:
                best = h; best["tsim"]=round(ts,2); break
        if best:
            decisions[num] = {"doi": best["doi"], "conf":"high", "method":"T3(vol+page+year)",
                              "match_title": best["title"], "cy": best["year"], "tsim": best.get("tsim")}
            print(f"   => ACCEPT {best['doi']} (sim={best.get('tsim')})")
        else:
            print("   => no confident match (leave pending)")
        time.sleep(0.8)
    # write decisions into cache
    for num, dec in decisions.items():
        dec.update(num=num, title=D.parse_ref(raws[num]).get("title",""),
                   year=D.parse_ref(raws[num]).get("year"),
                   kind=D.parse_ref(raws[num]).get("kind"),
                   journal=D.parse_ref(raws[num]).get("journal"))
        bynum[num] = dec
    json.dump(list(bynum.values()), open(RES,"w",encoding="utf-8"), ensure_ascii=False, indent=1)
    print("\nUpdated cache with", len(decisions), "vol-based resolutions:", sorted(decisions))

if __name__ == "__main__":
    main()
