#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Re-run the ORIGINAL working T1 (bibliographic + container-title filter,
verify vol+page+title) with a resilient fetch for the still-pending journal
refs. Plus special handling for #2/#41/#47. Writes verified DOIs into cache."""
import re, json, time, urllib.request, urllib.parse, difflib, sys
sys.path.insert(0, r"C:/D 盘\科研\虚拟敲除\gate1\rewiring_study")
import doi_lookup as D

HERE = r"C:/D 盘\科研\虚拟敲除\gate1\rewiring_study"
RES  = HERE + r"\doi_results.json"
MAIL = "ming.luo@med.org"

# resilient fetch -> monkeypatch doi_lookup.xfetch
def xfetch2(url, timeout=12, retries=5):
    last=None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent":"TSC-GNN-DOI/t1b (mailto:%s)"%MAIL})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.load(r)["message"]["items"]
        except Exception as e:
            last=e; time.sleep(1.5+i)
    return []
D.xfetch = xfetch2

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
    data = json.load(open(RES, encoding="utf-8"))
    bynum = {r["num"]: r for r in data}

    # --- special, evidence-backed overrides (manuscript citation has errors) ---
    # #2 Carmichael 2016: real paper is Neurotherapeutics 13:348-359 (NOT Nat Rev Neurosci)
    bynum[2]  = {"num":2, "doi":"10.1007/s13311-015-0408-0", "conf":"high",
                 "method":"T2(title-scan, author+year verified)", "match_title":
                 "The 3 Rs of Stroke Biology: Radial, Relayed, and Regenerative",
                 "cy":2016, "tsim":0.72,
                 "title":D.parse_ref(raws[2]).get("title"), "year":2016,
                 "kind":"journal", "journal":"Neurotherapeutics",
                 "flag":"MANUSCRIPT CITATION ERROR: journal should be Neurotherapeutics 13:348-359 (title '...Radial, Relayed, and Regenerative'), not Nat Rev Neurosci 17."}
    # #41 scPerturb: real paper is Nature Methods 21:531-540 (NOT Nat Biotechnol 42)
    bynum[41] = {"num":41, "doi":"10.1038/s41592-023-02144-y", "conf":"high",
                 "method":"T2(title-scan, author verified)", "match_title":
                 "scPerturb: harmonized single-cell perturbation data",
                 "cy":2024, "tsim":1.0,
                 "title":D.parse_ref(raws[41]).get("title"), "year":2024,
                 "kind":"journal", "journal":"Nature Methods",
                 "flag":"MANUSCRIPT CITATION ERROR: journal should be Nature Methods 21:531-540, not Nat Biotechnol 42:1311-1319."}
    # #47 Bhatt Statins: real paper is Circulation 135:1707-1720 (NOT Amarenco 2008 Future Lipidology)
    bynum[47] = {"num":47, "doi":"10.1161/CIRCULATIONAHA.116.023779", "conf":"high",
                 "method":"manual-verified (Circulation 135:1707-1720, 2017)",
                 "match_title":"Statins in Stroke Prevention", "cy":2017, "tsim":1.0,
                 "title":D.parse_ref(raws[47]).get("title"), "year":2017,
                 "kind":"journal", "journal":"Circulation",
                 "flag":"MANUSCRIPT: title 'Statins in stroke prevention' also matches Amarenco 2008 (Future Lipidology); correct target is Bhatt 2017 Circulation 135:1707-1720."}

    # --- T1 for remaining journal pending refs ---
    NUMS = [3,4,19,25,28,29,30,37,39,44,48,49,50,54]
    for num in NUMS:
        ref = D.parse_ref(raws[num])
        print(f"[{num:2d}] {ref.get('kind'):8s} {ref.get('journal','')} {ref.get('vol')}:{ref.get('pages')} ({ref.get('year')})")
        print(f"      {ref.get('title','')[:70]}")
        res = D.t1_journal_volume(ref)
        if res:
            res.update(num=num, title=ref.get("title",""), year=ref.get("year"),
                       kind=ref.get("kind"), journal=ref.get("journal"))
            bynum[num] = res
            print(f"   => T1 ACCEPT {res['doi']} (sim={res.get('tsim')}, '{res.get('match_title','')[:60]}')")
        else:
            # fallback: T2 but only show evidence, do NOT auto-accept
            r2 = D.t2_title_scan(ref)
            if r2:
                print(f"   T2 best: doi={r2['doi']} sim={r2.get('sim')} conf={r2['conf']} '{r2.get('match_title','')[:60]}'")
            else:
                print("   T2: no items")
            print("   => leave PENDING (needs manual)")
        time.sleep(0.8)

    json.dump(list(bynum.values()), open(RES,"w",encoding="utf-8"), ensure_ascii=False, indent=1)
    resolved = [k for k,v in bynum.items() if v.get("doi") and v.get("conf") in ("high","medium")]
    pending  = [k for k,v in bynum.items() if not (v.get("doi") and v.get("conf") in ("high","medium"))]
    print("\nresolved:", sorted(resolved))
    print("pending :", sorted(pending))

if __name__ == "__main__":
    main()
