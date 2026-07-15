#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
L5b · SigCom LINCS gene-perturbation directional consistency
-------------------------------------------------------------
CRISPR KO an activator TF => DOWN-regulate its targets => appears as reverser
OE       an activator TF => UP-regulate its targets   => appears as mimicker

KEY FIX (2026-07-09): database param must be STRING NAME ('l1000_xpr'),
  NOT library UUID — UUID causes HTTP 500 NullPointerException.
"""
import requests, json, os, time

HERE = os.path.dirname(os.path.abspath(__file__))
DORO = r"C:\D 盘\科研\虚拟敲除\gate1\data\dorothea\human_dorothea_regulon.tsv"
META = "https://maayanlab.cloud/sigcom-lincs/metadata-api"
DATA = "https://maayanlab.cloud/sigcom-lincs/data-api/api/v1"
OUT = os.path.join(HERE, "l5b_sigcom_result.json")

# data-api uses string names; metadata-api uses UUIDs
LIBS = {
    "CRISPR_KO":     {"name": "l1000_xpr", "uuid": "96c7b8c5-1eca-5764-88e4-e4ccaee6603f"},
    "Overexpression": {"name": "l1000_oe",  "uuid": "ef9389a8-53d3-50db-90cc-57e7d150b76c"},
}
TFS = ["SOX10", "CEBPB", "GATA2"]
CAP = 500

# ---- load human DoRothEA target sets ----
tgt = {}
with open(DORO, encoding="utf-8") as f:
    next(f)
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) < 2: continue
        tgt.setdefault(p[0], set()).add(p[1])
targets = {tf: sorted(tgt.get(tf, set())) for tf in TFS}

def convert(syms):
    r = requests.post(META + "/entities/find",
        json={"filter": {"where": {"meta.symbol": {"inq": list(syms)}}, "fields": ["id", "meta.symbol"]}},
        timeout=60)
    return {e["meta"]["symbol"]: e["id"] for e in r.json() if isinstance(e, dict) and "meta" in e}

def own_uuids(tf, lib_uuid):
    r = requests.post(META + "/signatures/find",
        json={"filter": {"where": {"meta.pert_name": tf, "library": lib_uuid}, "limit": 50}},
        timeout=60)
    return [s["id"] for s in r.json()]

def ranktwosided(up_ent, dn_ent, db_name, limit=100):
    q = {"up_entities": up_ent, "down_entities": dn_ent, "database": db_name, "limit": limit}
    backoff = [2, 4, 8, 15, 30]
    for attempt in range(5):
        try:
            r = requests.post(DATA + "/enrich/ranktwosided", json=q, timeout=180)
            if r.status_code == 200:
                return r.json()
            print(f"  [warn] attempt {attempt+1}/5 status {r.status_code}: {r.text[:120]}", flush=True)
        except Exception as e:
            print(f"  [warn] attempt {attempt+1}/5 exception: {e}", flush=True)
        time.sleep(backoff[min(attempt, len(backoff)-1)])
    raise RuntimeError("ranktwosided failed after 5 retries")

def names_for(uuids):
    if not uuids: return {}
    r = requests.post(META + "/signatures/find",
        json={"filter": {"where": {"id": {"inq": uuids}}}}, timeout=60)
    return {s["id"]: s.get("meta", {}).get("pert_name", "?") for s in r.json()}

# ---- preconvert targets ----
sym2id = {}
for tf in TFS:
    sym2id.update(convert(targets[tf]))

ent_targets = {}
for tf in TFS:
    e = [sym2id[s] for s in targets[tf] if s in sym2id]
    if len(e) > CAP:
        print(f"  [cap] {tf}: {len(e)} -> {CAP}", flush=True)
        e = e[:CAP]
    ent_targets[tf] = e

# neutral down pool
NEUTRAL_POOL = ["HBB","HBA1","INS","ALB","TFRC","FGB","HPX","SERPINA1",
                "APOA1","AMBP","ORM1","C3","TF","ACTB","GAPDH",
                "RPL13A","B2M","SDHA","ENO1","TUBB","PGK1","HPRT1"]
neut_map = convert(NEUTRAL_POOL)
print("converted targets:", {tf: len(ent_targets[tf]) for tf in TFS}, flush=True)

def neutral_down_for(tf):
    excl = set(ent_targets[tf])
    return [eid for s, eid in neut_map.items() if s not in TFS and eid not in excl][:5]

# ---- own uuids per (tf, lib) ----
own = {}
for lib, info in LIBS.items():
    own[lib] = {}
    for tf in TFS:
        own[lib][tf] = own_uuids(tf, info["uuid"])
        print(f"  own[{lib}][{tf}] = {len(own[lib][tf])} sigs", flush=True)
        time.sleep(0.5)

# ---- main loop ----
results = {}
failed_any = False
n_failed = 0
for lib, info in LIBS.items():
    print(f"\n===== {lib} (db={info['name']}) =====", flush=True)
    qres = {}
    for tf in TFS:
        dn = neutral_down_for(tf)
        print(f"  querying {tf} (up={len(ent_targets[tf])}, dn={len(dn)}) ...", flush=True)
        try:
            j = ranktwosided(ent_targets[tf], dn, info["name"], limit=2000)
        except RuntimeError as e:
            print(f"  [FAIL] {tf}: {e}", flush=True)
            qres[tf] = {"raw": [], "stat": {}, "n_rev": 0, "n_mim": 0, "error": str(e)}
            failed_any = True
            n_failed += 1
            time.sleep(3)
            continue
        res = j.get("results", [])
        n_rev = j.get("reversers", 0)
        n_mim = j.get("mimickers", 0)
        max_rank = j.get("maxRank", "?")
        stat = {x["uuid"]: x for x in res}
        qres[tf] = {"raw": res, "stat": stat, "n_rev": n_rev, "n_mim": n_mim, "max_rank": max_rank}
        print(f"  -> {tf}: {len(res)} results ({n_rev} rev, {n_mim} mim, maxRank={max_rank})", flush=True)
        time.sleep(2)

    # cross matrix
    mat = {}
    for tfq in TFS:
        for tfo in TFS:
            uus = own[lib].get(tfo, [])
            hits = []
            for u in uus:
                if u in qres[tfq]["stat"]:
                    s = qres[tfq]["stat"][u]
                    hits.append({
                        "rank": s.get("rank"), "type": s.get("type"),
                        "p_up": s.get("p-up"), "p_down": s.get("p-down"),
                        "fdr_up": s.get("fdr-up"), "fdr_down": s.get("fdr-down"),
                        "z_up": s.get("z-up"), "z_down": s.get("z-down"),
                        "z_sum": s.get("z-sum"),
                    })
            mat[(tfq, tfo)] = hits
    results[lib] = {"matrix": {f"{a}|{b}": mat[(a,b)] for (a,b) in mat},
                    "n_rev": {tf: qres[tf]["n_rev"] for tf in TFS},
                    "n_mim": {tf: qres[tf]["n_mim"] for tf in TFS},
                    "max_rank": {tf: qres[tf].get("max_rank") for tf in TFS}}

    print(f"\n  {'query':8s} | " + " | ".join(f"{tfo:20s}" for tfo in TFS), flush=True)
    for tfq in TFS:
        row = []
        for tfo in TFS:
            h = mat[(tfq, tfo)]
            if h:
                s = h[0]
                t = 'REV' if s['type'] == 'reversers' else 'MIM'
                pv = s['p_down'] if s['type'] == 'reversers' else s['p_up']
                pvs = f"{pv:.1e}" if pv is not None else "N/A"
                row.append(f"rk{s['rank']:4d}/{t}/p={pvs}")
            else:
                row.append("absent")
        print(f"  {tfq:8s} | " + " | ".join(f"{x:20s}" for x in row), flush=True)

def safe_write(path, obj):
    """Atomic write: never leave a half-written / corrupt JSON on disk."""
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)
    os.replace(tmp, path)

doc = {
    "targets_n": {tf: len(targets[tf]) for tf in TFS},
    "ent_targets_n": {tf: len(ent_targets[tf]) for tf in TFS},
    "own_uuids_n": {lib: {tf: len(own[lib][tf]) for tf in TFS} for lib in LIBS},
    "results": results,
}

if failed_any:
    # HARDENING (CONFORMANCE audit B2, 2026-07-14): a failed run must NOT
    # clobber a previously good result. Write the partial output to a
    # sidecar file and leave OUT untouched.
    partial = OUT + ".partial"
    safe_write(partial, doc)
    print(f"\n[WARN] {n_failed} query(ies) failed (see log above).", flush=True)
    print(f"[WARN] The existing good result at {OUT} was NOT overwritten.", flush=True)
    print(f"[WARN] Partial/failed output written to {partial}", flush=True)
    print(f"[WARN] Re-run only after the SigCom LINCS API recovers (HTTP 500).", flush=True)
else:
    safe_write(OUT, doc)
    print(f"\n[done] -> {OUT}", flush=True)
