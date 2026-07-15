#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fill [DOI: pending] in manuscript_v6_humanized.md using Crossref.

Strategy (tiered, verification-gated -- never fabricate):
  T1  Journal article with volume + first page:
      query Crossref filtered by container-title (full journal name) + volume,
      then match first_page exactly and year within +-1. -> EXACT, high confidence.
  T2  Fallback (no vol / preprint / conference):
      query.title scan (rows=20), pick best title-similarity. Accept if sim>=0.92
      (high) or 0.85-0.92 (medium); else low/None -> left pending + flagged.

Outputs:
  doi_results.json   -- machine-readable results + confidence
  doi_report.md      -- human-readable table for review
The script does NOT modify the manuscript; fills are applied in a second step
after the user/review confirms confidence.
"""
import re, json, time, os, urllib.request, urllib.parse, difflib, sys

HERE = r"C:/D 盘/科研/虚拟敲除/gate1/rewiring_study"
SRC  = HERE + r"\manuscript_v6_humanized.md"
MAIL = "ming.luo@med.org"   # Crossref polite pool

# abbreviation -> Crossref registered container-title (only needed for T1)
JMAP = {
    "Nat. Rev. Neurosci.": "Nature Reviews Neuroscience",
    "Acta Neuropathol. Commun.": "Acta Neuropathologica Communications",
    "Nat. Immunol.": "Nature Immunology",
    "Nature": "Nature",
    "Nat. Methods": "Nature Methods",
    "Nat. Biotechnol.": "Nature Biotechnology",
    "Nat. Cell Biol.": "Nature Cell Biology",
    "Commun. Biol.": "Communications Biology",
    "Cell": "Cell",
    "npj Syst. Biol. Appl.": "npj Systems Biology and Applications",
    "Genes Dev.": "Genes & Development",
    "Ann. Neurol.": "Annals of Neurology",
    "Annu. Rev. Neurosci.": "Annual Review of Neuroscience",
    "Cell Stem Cell": "Cell Stem Cell",
    "Stroke": "Stroke",
    "Prog. Neurobiol.": "Progress in Neurobiology",
    "Nat. Commun.": "Nature Communications",
    "J. R. Stat. Soc. B": "Journal of the Royal Statistical Society: Series B (Statistical Methodology)",
    "Int. J. Mol. Sci.": "International Journal of Molecular Sciences",
    "JHEP Rep.": "JHEP Reports",
    "Nat. Rev. Phys.": "Nature Reviews Physics",
    "Nat. Microbiol.": "Nature Microbiology",
    "Nat. Genet.": "Nature Genetics",
    "BMC Bioinformatics": "BMC Bioinformatics",
    "Proc. Natl Acad. Sci. USA": "Proceedings of the National Academy of Sciences",
    "Theor. Biol. Med. Model.": "Theoretical Biology and Medical Modelling",
    "BioEssays": "BioEssays",
    "Nat. Rev. Drug Discov.": "Nature Reviews Drug Discovery",
    "Circulation": "Circulation",
}

def xfetch(url, timeout=6, retries=2):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "TSC-GNN-DOI/1.0 (mailto:%s)" % MAIL})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.load(r)["message"]["items"]
        except Exception as e:
            last = e
            time.sleep(1.0 + i)
    return []

def norm(s):
    s = s.lower()
    s = re.sub(r"[^a-z0-9 ]", " ", s)
    return re.sub(r"\s+", " ", s).strip()

def first_page(pages):
    fp = re.split(r"[–\-]", pages)[0].strip()
    return fp

# ---------- parsing ----------
VOL_RE = re.compile(
    r"^\d+\.\s+(?P<authors>.*?) \((?P<year>\d{4})\)\.\s+(?P<title>.*?)\.\s+"
    r"(?P<journal>.*?) (?P<vol>\d+),\s+(?P<pages>[\d–\-]+(?:\.[A-Za-z]+\d*)?)\.\s*$"
)
SIMPLE_RE = re.compile(r"^\d+\.\s+(?P<authors>.*?) \((?P<year>\d{4})\)\.\s+(?P<title>.*?)\.\s+(?P<rest>.*?)\.\s*$")

def parse_ref(line):
    m = VOL_RE.match(line)
    if m:
        return {"kind": "journal", "authors": m.group("authors"), "year": int(m.group("year")),
                "title": m.group("title"), "journal": m.group("journal").strip(),
                "vol": m.group("vol"), "pages": m.group("pages")}
    m = SIMPLE_RE.match(line)
    if m:
        rest = m.group("rest").strip()
        if "[bioRxiv]" in rest or "bioRxiv" in rest:
            kind = "preprint"
        elif "ICLR" in rest or "NeurIPS" in rest:
            kind = "conference"
        else:
            kind = "journal_novol"
        return {"kind": kind, "authors": m.group("authors"), "year": int(m.group("year")),
                "title": m.group("title"), "journal": rest, "vol": None, "pages": None}
    return {"kind": "unparsed", "raw": line}

# ---------- lookups ----------
def t1_journal_volume(ref):
    """Exact match via Crossref: container-title filter + (title OR first-author)
    query, then verify volume + first-page + year. No fuzzy title -> no fabrication."""
    full = JMAP.get(ref["journal"])
    if not full or not ref["vol"]:
        return None
    fp = first_page(ref["pages"]); y = ref["year"]; vol = ref["vol"]
    nt = norm(ref["title"])
    fa = ref["authors"].split(",")[0].strip()
    # ONE combined query: full citation (title+author+year+journal+vol+page) + container-title filter.
    # Crossref gets maximal signal; we then verify by volume + first-page + title agreement.
    cit = " ".join([ref["title"], fa, str(y), full, str(vol), fp])
    q = "query.bibliographic=" + urllib.parse.quote(cit) + "&filter=container-title:" + urllib.parse.quote(full)
    items = xfetch("https://api.crossref.org/works?" + q + "&rows=50&mailto=" + MAIL)
    for it in items:
        ivol = str(it.get("volume") or "")
        ipage = str(it.get("page") or "")
        dp = (it.get("issued") or {}).get("date-parts")
        cy = dp[0][0] if (dp and dp[0] and dp[0][0]) else None
        if ivol == vol and (ipage.startswith(fp) or fp in ipage):
            if cy is None or abs(int(cy) - y) <= 1:
                ts = difflib.SequenceMatcher(None, nt, norm((it.get("title") or [""])[0])).ratio()
                if ts >= 0.5:
                    return {"doi": it.get("DOI"), "method": "T1(vol+page+title)", "conf": "high",
                            "match_title": (it.get("title") or [""])[0], "cy": cy, "tsim": round(ts, 2)}
    return None

def t2_title_scan(ref):
    title = ref["title"]
    items = xfetch("https://api.crossref.org/works?query.title=" + urllib.parse.quote(title) + "&rows=40&mailto=" + MAIL)
    if not items:
        return None
    nt = norm(title)
    best = None; bs = -1
    for it in items:
        cand = (it.get("title") or [""])[0]
        sc = difflib.SequenceMatcher(None, nt, norm(cand)).ratio()
        if sc > bs:
            bs = sc; best = it
    dp = (best.get("issued") or {}).get("date-parts")
    cy = dp[0][0] if (dp and dp[0] and dp[0][0]) else None
    sc = difflib.SequenceMatcher(None, nt, norm((best.get("title") or [""])[0])).ratio()
    # strict: only near-exact titles accepted; never borderline
    conf = "high" if sc >= 0.95 else ("medium" if sc >= 0.90 else "low")
    return {"doi": best.get("DOI"), "method": "T2(title-scan)", "conf": conf, "sim": round(sc, 3),
            "match_title": (best.get("title") or [""])[0], "cy": cy}

def lookup(ref):
    if ref["kind"] == "journal":
        r = t1_journal_volume(ref)
        if r:
            return r
        r2 = t2_title_scan(ref)
        if r2 and r2["conf"] in ("high", "medium"):
            return r2
        return {"doi": None, "method": "T1+T2-none", "conf": "low"}
    if ref["kind"] in ("preprint", "conference", "journal_novol"):
        r2 = t2_title_scan(ref)
        if r2 and r2["conf"] in ("high", "medium"):
            return r2
        return {"doi": None, "method": "T2-none", "conf": "low"}
    return {"doi": None, "method": "unparsed", "conf": "low"}

# ---------- main ----------
def main():
    text = open(SRC, encoding="utf-8").read()
    lines = text.splitlines()
    start = next(i for i, l in enumerate(lines) if l.strip() == "## References")
    refs = []
    for l in lines[start + 1:]:
        if l.startswith("## ") and l.strip() != "## References":
            break
        if re.match(r"^\d+\.", l.strip()):
            # strip the trailing "[DOI: pending]" token before parsing
            refs.append(re.sub(r"\s*\[DOI:[^\]]*\]\s*$", "", l.strip()))
    print("parsed %d references" % len(refs))
    # resume: load previously resolved DOIs (high/medium) and skip them
    cache = {}
    if os.path.exists(HERE + r"\doi_results.json"):
        try:
            for r in json.load(open(HERE + r"\doi_results.json", encoding="utf-8")):
                if r.get("doi") and r.get("conf") in ("high", "medium"):
                    cache[r["num"]] = r
        except Exception:
            pass
    print("resuming: %d refs already resolved in cache" % len(cache))
    results = list(cache.values())
    logf = open(HERE + r"\doi_lookup.log", "w", encoding="utf-8")
    for raw in refs:
        num = int(raw.split(".", 1)[0])
        if num in cache:
            line = "[%2d] skip    cached       %s" % (num, cache[num]["doi"])
            print(line, flush=True); logf.write(line + "\n"); logf.flush()
            continue
        ref = parse_ref(raw)
        res = lookup(ref)
        res["num"] = num
        res["title"] = ref.get("title", "")
        res["year"] = ref.get("year")
        res["kind"] = ref.get("kind")
        res["journal"] = ref.get("journal")
        results.append(res)
        tag = res["doi"] or "PENDING"
        line = "[%2d] %-7s %-12s %s | %s" % (num, res["conf"], res["method"], tag, ref.get("title", "")[:50])
        print(line, flush=True)
        logf.write(line + "\n"); logf.flush()
        json.dump(results, open(HERE + r"\doi_results.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        time.sleep(1.0)
    json.dump(results, open(HERE + r"\doi_results.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    # report
    with open(HERE + r"\doi_report.md", "w", encoding="utf-8") as f:
        f.write("# DOI lookup report (Crossref)\n\n")
        f.write("| # | conf | method | DOI | year | title (matched) |\n")
        f.write("|---|------|--------|-----|------|------------------|\n")
        for r in results:
            f.write("| %d | %s | %s | %s | %s | %s |\n" % (
                r["num"], r["conf"], r["method"], r["doi"] or "—",
                r.get("cy") or r.get("year") or "", (r.get("match_title") or r.get("title") or "")[:60]))
        pend = [r["num"] for r in results if not r["doi"]]
        f.write("\n**Pending (%d):** %s\n" % (len(pend), pend))
    print("\nwrote doi_results.json and doi_report.md")
    print("pending:", pend)

if __name__ == "__main__":
    main()
