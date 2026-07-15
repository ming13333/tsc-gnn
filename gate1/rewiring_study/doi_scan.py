#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Strict resolver: scan rows=80 of a title bibliographic query, accept ONLY an
item whose container-title matches the manuscript journal (via JMAP), volume
matches, first-page within +-3, year within +-1, and first-author lastname is
in the authors. Prints full metadata for verification. Overwrites any prior
(decoy) entry for the target nums, so it also CLEANS false accepts."""
import re, json, time, urllib.request, urllib.parse, difflib, sys
sys.path.insert(0, r"C:/D 盘\科研\虚拟敲除\gate1\rewiring_study")
import doi_lookup as D

HERE = r"C:/D 盘\科研\虚拟敲除\gate1\rewiring_study"
RES  = HERE + r"\doi_results.json"
MAIL = "ming.luo@med.org"

def xfetch(url, timeout=12, retries=5):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent":"TSC-DOI/scan2 (mailto:%s)"%MAIL})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.load(r)["message"]["items"]
        except Exception:
            time.sleep(1.5+i)
    return []

def norm(s):
    s=(s or "").lower(); s=re.sub(r"[^a-z0-9 ]"," ",s); return re.sub(r"\s+"," ",s).strip()
def fp(p):
    try: return int(re.split(r"[–\-]", str(p))[0].strip())
    except: return None
def fa_last(raw):
    a=raw.split(",")[0].strip().lower(); return re.sub(r"[^a-z]","",a)

def resolve(ref, rows=80):
    full=D.JMAP.get(ref.get("journal"))
    if not full or not ref.get("vol"): return None
    fa=fa_last(ref.get("authors","")); y=ref.get("year")
    fp_ms=fp(ref.get("pages")); vol=ref.get("vol"); nt=norm(ref.get("title",""))
    items=xfetch("https://api.crossref.org/works?query.bibliographic="+
                 urllib.parse.quote(ref["title"])+f"&rows={rows}&mailto={MAIL}")
    best=None; bscore=-1
    for it in items or []:
        ct=(it.get("container-title") or [""])[0]
        if norm(full) not in norm(ct): continue
        if str(it.get("volume") or "") != str(vol): continue
        ip=fp(it.get("page") or "")
        if ip is None or fp_ms is None or abs(ip-fp_ms)>3: continue
        dp=(it.get("issued") or {}).get("date-parts")
        cy=dp[0][0] if (dp and dp[0] and dp[0][0]) else None
        if cy is None or abs(int(cy)-y)>1: continue
        au=", ".join((a.get("family") or "").lower() for a in it.get("author",[]))
        if not (fa and fa in au): continue
        mt=(it.get("title") or [""])[0]
        ts=difflib.SequenceMatcher(None,nt,norm(mt)).ratio()
        if ts>bscore: bscore=ts; best=it
    return best

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
    NUMS=[1,3,4,10,19,25,28,29,30,37,39,40,44,47,48,49,50,54,38]
    for num in NUMS:
        ref=D.parse_ref(raws[num])
        print(f"[{num:2d}] {ref.get('journal','')} {ref.get('vol')}:{ref.get('pages')} ({ref.get('year')}) {ref.get('title','')[:52]}")
        best=resolve(ref)
        if best:
            mt=(best.get("title") or [""])[0]
            au=", ".join(a.get("family") or "" for a in best.get("author",[])[:4])
            ct=(best.get("container-title") or [""])[0]
            dp=(best.get("issued") or {}).get("date-parts"); cy=dp[0][0] if (dp and dp[0] and dp[0][0]) else None
            ts=difflib.SequenceMatcher(None,norm(ref["title"]),norm(mt)).ratio()
            bynum[num]={"num":num,"doi":best.get("DOI"),"conf":"high","method":"T5(strict:journal+vol+page+year+author)",
                        "match_title":mt,"cy":cy,"tsim":round(ts,2),
                        "title":ref.get("title",""),"year":ref.get("year"),
                        "kind":ref.get("kind"),"journal":ref.get("journal")}
            print(f"   => {best.get('DOI')} sim={ts:.2f} | {au} | {ct} {best.get('volume')}:{best.get('page')} ({cy})")
        else:
            # overwrite any prior decoy with pending
            bynum[num]={"num":num,"doi":None,"conf":"low","method":"T5-none",
                        "title":ref.get("title",""),"year":ref.get("year"),
                        "kind":ref.get("kind"),"journal":ref.get("journal")}
            print("   => PENDING (no strict match)")
        time.sleep(0.8)
    json.dump(list(bynum.values()),open(RES,"w",encoding="utf-8"),ensure_ascii=False,indent=1)
    resolved=[k for k,v in bynum.items() if v.get("doi") and v.get("conf") in ("high","medium")]
    pending=[k for k,v in bynum.items() if not (v.get("doi") and v.get("conf") in ("high","medium"))]
    print("\nresolved:",sorted(resolved))
    print("pending :",sorted(pending))

if __name__=="__main__":
    main()
