#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Definitive safety check: for EVERY cached DOI, fetch Crossref metadata and
verify it matches the manuscript reference (title similarity / first-author /
year). Flags any mismatch so wrong DOIs are never applied. Also removes the
two known-wrong conference DOIs (#16, #18) from the cache."""
import re, json, time, urllib.request, urllib.parse, difflib, sys
sys.path.insert(0, r"C:/D 盘\科研\虚拟敲除\gate1\rewiring_study")
import doi_lookup as D

HERE = r"C:/D 盘\科研\虚拟敲除\gate1\rewiring_study"
RES  = HERE + r"\doi_results.json"
MAIL = "ming.luo@med.org"

def meta(doi):
    url="https://api.crossref.org/works/"+urllib.parse.quote(doi)
    for i in range(6):
        try:
            req=urllib.request.Request(url, headers={"User-Agent":"verify (mailto:%s)"%MAIL})
            with urllib.request.urlopen(req, timeout=12) as r:
                it=json.load(r)["message"]
                return {"doi":it.get("DOI"),"title":(it.get("title") or [""])[0],
                        "authors":[ (a.get("family") or "").lower() for a in it.get("author",[]) ],
                        "journal":(it.get("container-title") or [""])[0],
                        "vol":it.get("volume"),"page":it.get("page"),
                        "year":(it.get("issued") or {}).get("date-parts",[["?"]])[0][0]}
        except Exception:
            time.sleep(1.5+i)
    return None

def norm(s):
    s=(s or "").lower(); s=re.sub(r"[^a-z0-9 ]"," ",s); return re.sub(r"\s+"," ",s).strip()
def fa_last(raw):
    a=raw.split(",")[0].strip().lower(); return re.sub(r"[^a-z]","",a)

def main():
    text=open(D.SRC,encoding="utf-8").read()
    lines=text.splitlines()
    start=next(i for i,l in enumerate(lines) if l.strip()=="## References")
    raws={}
    for l in lines[start+1:]:
        if l.startswith("## ") and l.strip()!="## References": break
        if re.match(r"^\d+\.",l.strip()):
            raw=re.sub(r"\s*\[DOI:[^\]]*\]\s*$","",l.strip())
            raws[int(raw.split(".",1)[0])]=raw
    data=json.load(open(RES,encoding="utf-8"))
    bynum={r["num"]:r for r in data}
    problems=[]
    print("# Comprehensive DOI verification\n")
    for num in sorted(bynum):
        e=bynum[num]
        doi=e.get("doi")
        if not doi: continue
        ref=D.parse_ref(raws[num])
        m=meta(doi)
        if not m:
            print(f"[{num:2d}] {doi} -> FETCH FAILED (mark pending)")
            bynum[num]={"num":num,"doi":None,"conf":"low","method":"verify-fetch-fail",
                        "title":ref.get("title",""),"year":ref.get("year"),
                        "kind":ref.get("kind"),"journal":ref.get("journal")}
            problems.append((num,"fetch-fail")); continue
        ts=difflib.SequenceMatcher(None,norm(ref["title"]),norm(m["title"])).ratio()
        fa=fa_last(ref.get("authors",""))
        auth_ok = fa and any(fa in (a or "") for a in m["authors"])
        # for conference/preprint, authors may be sparse; rely on title+author
        yr=ref.get("year"); yrok = (m["year"]=="?" ) or abs(int(m["year"])-int(yr))<=2
        # accept if strong title match OR (author match AND year ok)
        ok = (ts>=0.80) or (auth_ok and yrok and ts>=0.55)
        flag=""
        if not yrok: flag+=" YEAR-MISMATCH(%s vs %s)"%(yr,m["year"])
        if not auth_ok and ts<0.80: flag+=" AUTHOR-CHECK"
        status="OK" if ok else "MISMATCH"
        if status=="MISMATCH": problems.append((num,flag or "low-sim"))
        print(f"[{num:2d}] {status:9s} sim={ts:.2f} auth={auth_ok} {flag}")
        print(f"      ms: {ref['title'][:48]} | {ref.get('authors','')[:22]} ({yr})")
        print(f"      cx: {m['title'][:48]} | {', '.join(m['authors'][:3])} ({m['year']}) {m['journal']} {m['vol']}:{m['page']}")
    # remove known-wrong conference DOIs
    for bad in (16,18):
        if bad in bynum and bynum[bad].get("doi"):
            print(f"\n[!] Removing known-wrong cached DOI for #{bad} ({bynum[bad]['doi']}) -> pending")
            bynum[bad]={"num":bad,"doi":None,"conf":"low","method":"removed-wrong-conference",
                        "title":D.parse_ref(raws[bad]).get("title",""),
                        "year":D.parse_ref(raws[bad]).get("year"),
                        "kind":D.parse_ref(raws[bad]).get("kind"),
                        "journal":D.parse_ref(raws[bad]).get("journal")}
            problems.append((bad,"removed-wrong"))
    json.dump(list(bynum.values()),open(RES,"w",encoding="utf-8"),ensure_ascii=False,indent=1)
    print("\n=== PROBLEMS (do NOT apply these) ===")
    for p in problems: print("  ",p)
    resolved=[k for k,v in bynum.items() if v.get("doi") and v.get("conf") in ("high","medium")]
    pending=[k for k,v in bynum.items() if not (v.get("doi") and v.get("conf") in ("high","medium"))]
    print("\nresolved:",sorted(resolved))
    print("pending :",sorted(pending))

if __name__=="__main__":
    main()
