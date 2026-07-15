#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Apply resolved DOIs from doi_results.json into manuscript_v6_humanized.md,
replacing '[DOI: pending]' with '[DOI: <doi>]' per reference number.

Only high/medium-confidence DOIs are applied. A backup of the manuscript is
written first. Prints a summary of applied / skipped.
"""
import json, re, shutil, os

HERE = r"C:/D 盘\科研\虚拟敲除\gate1\rewiring_study"
SRC  = HERE + r"\manuscript_v6_humanized.md"
RES  = HERE + r"\doi_results.json"
BAK  = HERE + r"\manuscript_v6_humanized.md.bak_dois"

def main():
    res = json.load(open(RES, encoding="utf-8"))
    doi_by_num = {r["num"]: r["doi"] for r in res if r.get("doi") and r.get("conf") in ("high", "medium")}
    text = open(SRC, encoding="utf-8").read()
    shutil.copyfile(SRC, BAK)
    print("backup ->", BAK)

    # replace [DOI: pending] on each numbered reference line
    def repl(m):
        num = int(m.group(1))
        doi = doi_by_num.get(num)
        if doi:
            return "%s[DOI: %s]" % (m.group(2), doi)
        return m.group(0)  # leave pending
    pat = re.compile(r"^(\d+)\.\s+(.*?)\[DOI: pending\]$", re.M)
    new_text, n = pat.subn(repl, text)
    applied = n
    pending_left = len(re.findall(r"\[DOI: pending\]", new_text))
    open(SRC, "w", encoding="utf-8").write(new_text)
    print("applied %d DOIs; [DOI: pending] remaining: %d" % (applied, pending_left))
    print("resolved nums:", sorted(doi_by_num.keys()))
    # list any resolved nums whose line still shows pending (shouldn't happen)
    if pending_left:
        print("WARNING: some [DOI: pending] remain (conference preprints / unverified).")

if __name__ == "__main__":
    main()
