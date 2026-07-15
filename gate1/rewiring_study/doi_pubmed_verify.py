#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Verify & recover the 23 pending DOIs via PubMed E-utilities + Crossref direct check.
Corporate proxy MITMs TLS, so we use an unverified SSL context for these read-only public APIs.
NEVER inject a DOI that is not confirmed by matching metadata (vol+page or title+author+year).
"""
import urllib.request, urllib.parse, ssl, json, time, re, sys

ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
HDR = {'User-Agent': 'tsc-gnn-ref-verify/1.0 (mailto:hwenwu321@gmail.com)'}

def get(url):
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers=HDR)
            return urllib.request.urlopen(req, timeout=30, context=ctx).read().decode('utf-8', 'replace')
        except Exception as e:
            if attempt == 2:
                return f"__ERROR__:{e}"
            time.sleep(2)

def esearch(term):
    u = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=8&term=' + urllib.parse.quote(term)
    r = get(u)
    if r.startswith('__ERROR__'): return []
    try:
        return json.loads(r)['esearchresult'].get('idlist', [])
    except Exception:
        return []

def efetch(pmids):
    if not pmids: return ''
    u = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&id=' + ','.join(pmids)
    return get(u)

def parse_articles(xml):
    """Return list of dicts from a PubmedArticleSet XML (crude regex parse, good enough for fields)."""
    arts = re.split(r'<PubmedArticle>', xml)[1:]
    out = []
    for a in arts:
        d = {}
        m = re.search(r'<PMID[^>]*>(\d+)</PMID>', a); d['pmid'] = m.group(1) if m else ''
        m = re.search(r'<ArticleTitle>(.*?)</ArticleTitle>', a, re.S); d['title'] = re.sub(r'<[^>]+>', '', m.group(1)).strip() if m else ''
        m = re.search(r'<Volume>(.*?)</Volume>', a, re.S); d['vol'] = m.group(1).strip() if m else ''
        m = re.search(r'<MedlinePgn>(.*?)</MedlinePgn>', a, re.S); d['page'] = m.group(1).strip() if m else ''
        m = re.search(r'<ISOAbbreviation>(.*?)</ISOAbbreviation>', a, re.S); d['journal'] = m.group(1).strip() if m else ''
        m = re.search(r'<Year>(\d{4})</Year>', a); d['year'] = m.group(1) if m else ''
        m = re.search(r'<ArticleId IdType="doi">(.*?)</ArticleId>', a, re.S); d['doi'] = m.group(1).strip().lower() if m else ''
        # first author last name
        m = re.search(r'<Author[^>]*>\s*<LastName>(.*?)</LastName>', a, re.S); d['first_author'] = m.group(1).strip() if m else ''
        out.append(d)
    return out

def crossref(doi):
    r = get('https://api.crossref.org/works/' + urllib.parse.quote(doi))
    if r.startswith('__ERROR__'): return None
    try:
        m = json.loads(r)['message']
        return {'title': (m.get('title') or [''])[0], 'journal': (m.get('container-title') or [''])[0],
                'vol': m.get('volume', ''), 'page': m.get('page', ''),
                'year': str((m.get('issued', {}).get('date-parts', [['']])[0][0]))}
    except Exception:
        return None

# ---- The 23 pending references: identifying metadata ----
# type: 'journal' -> PubMed; 'arxiv' -> assign arXiv DOI; 'flag' -> needs human decision
REFS = [
    dict(n=1, typ='journal', au='Cramer', yr='2008', ta='Nat Rev Neurosci', vol='9', pg='720', kw='stroke recovery mechanisms'),
    dict(n=3, typ='journal', au='Li', yr='2021', ta='Acta Neuropathol Commun', vol='9', pg='152', kw='single-cell ischemic stroke transcriptomic'),
    dict(n=4, typ='journal', au='Anrather', yr='2024', ta='Nat Immunol', vol='25', pg='294', kw='single-cell stroke brain blood immune'),
    dict(n=15, typ='arxiv', arxiv='1609.02907', desc='Kipf & Welling GCN ICLR 2017'),
    dict(n=16, typ='arxiv', arxiv='1710.10903', desc='Velickovic GAT ICLR 2018'),
    dict(n=17, typ='arxiv', arxiv='1706.02216', desc='Hamilton GraphSAGE NeurIPS 2017'),
    dict(n=18, typ='flag', desc='Dwivedi 2023 "GNN: a review of methods and applications" — title matches Zhou 2020 AI Open (10.1016/j.aiopen.2021.01.001); Dwivedi 2023 real work is "Benchmarking GNNs" JMLR. AMBIGUOUS.'),
    dict(n=19, typ='journal', au='Lotfollahi', yr='2019', ta='Nat Methods', vol='16', pg='1253', kw='scGen predicting cellular responses perturbations'),
    dict(n=25, typ='journal', au='Fancy', yr='2011', ta='Ann Neurol', vol='69', pg='579', kw='myelin regeneration cellular molecular'),
    dict(n=28, typ='journal', au='Doyle', yr='2008', ta='Neuropharmacology', vol='55', pg='310', kw='mechanisms ischemic brain damage repair'),
    dict(n=29, typ='journal', au='Jin', yr='2010', ta='J Leukoc Biol', vol='87', pg='779', kw='inflammatory mechanisms ischemic stroke'),
    dict(n=30, typ='journal', au='Barr', yr='2010', ta='Neurology', vol='75', pg='1009', kw='genomic biomarkers acute stroke'),
    dict(n=37, typ='flag', desc='Velickovic 2023 Nat Rev Phys "Message passing all you need" 5,343-356 — title likely fabricated; verify.'),
    dict(n=38, typ='arxiv', arxiv='1803.01422', desc='Zheng NO TEARS NeurIPS 2018'),
    dict(n=39, typ='journal', au='Lotfollahi', yr='2023', ta='Nat Biotechnol', vol='41', pg='1759', kw='biologically informed deep learning perturbation'),
    dict(n=40, typ='biorxiv', desc='Bereket & Karaletsos PerturbNet 2022'),
    dict(n=42, typ='flag', desc='Singh 2024 PerturbExpress Nat Microbiol — no vol/page given; may not exist. Verify.'),
    dict(n=44, typ='journal', au='Schraivogel', yr='2020', ta='Nat Methods', vol='17', pg='629', kw='targeted Perturb-seq genome-scale regulators'),
    dict(n=47, typ='flag', desc='Bhatt 2017 "Statins in stroke prevention" Circulation 135,1707-1720 — verify; Circulation 135 exists.'),
    dict(n=48, typ='flag', desc='Li 2020 "HDAC inhibitors in neurological diseases" Nat Rev Drug Discov 19,341-359 — verify title/journal.'),
    dict(n=49, typ='journal', au='Ahlmann-Eltze', yr='2025', ta='Nat Methods', vol='22', pg='322', kw='perturbation single-cell benchmark linear baselines'),
    dict(n=50, typ='journal', au='Kartha', yr='2023', ta='Nat Genet', vol='55', pg='1339', kw='functional inference gene regulation Perturb-seq'),
    dict(n=54, typ='flag', desc='Weider 2021 Nat Commun 12,4240 "Nfat/calcineurin oligodendrocyte" — high-cited version is 2018 Nat Commun 9 (10.1038/s41467-018-03336-3). Verify.'),
]

def norm(s):
    return re.sub(r'[^a-z0-9]', '', (s or '').lower())

results = []
for ref in REFS:
    n = ref['n']; typ = ref['typ']
    print(f"\n===== [{n}] type={typ} =====", flush=True)
    if typ == 'arxiv':
        doi = '10.48550/arXiv.' + ref['arxiv']
        cr = crossref(doi); time.sleep(0.5)
        ok = cr is not None
        print(f"  arXiv DOI {doi}  crossref_title={cr['title'] if cr else 'N/A'}")
        results.append(dict(n=n, status='arxiv_ok' if ok else 'arxiv_unchecked', doi=doi, note=ref['desc'], cr=cr))
        continue
    if typ in ('flag', 'biorxiv'):
        # still try a pubmed search to help
        au = ref.get('au'); found = []
        if 'Bhatt' in ref.get('desc',''):
            found = parse_articles(efetch(esearch('Bhatt[Author] AND 2017[PDAT] AND Circulation[TA] AND statin*[Title]')))
        elif 'HDAC' in ref.get('desc',''):
            found = parse_articles(efetch(esearch('HDAC inhibitor* neurological disease* review'))) 
        elif 'Weider' in ref.get('desc',''):
            found = parse_articles(efetch(esearch('Weider[Author] AND oligodendrocyte[Title]')))
        elif 'PerturbExpress' in ref.get('desc',''):
            found = parse_articles(efetch(esearch('PerturbExpress perturbation')))
        elif 'Velickovic 2023' in ref.get('desc','') or 'Nat Rev Phys' in ref.get('desc',''):
            found = parse_articles(efetch(esearch('Velickovic[Author] AND 2023[PDAT]')))
        elif 'PerturbNet' in ref.get('desc',''):
            found = parse_articles(efetch(esearch('PerturbNet perturbation Bereket')))
        time.sleep(0.5)
        print(f"  FLAG: {ref['desc']}")
        for f in found[:5]:
            print(f"    cand pmid={f['pmid']} {f['first_author']} {f['year']} {f['journal']} {f['vol']}:{f['page']} doi={f['doi']}  | {f['title'][:70]}")
        results.append(dict(n=n, status='flag', doi='', note=ref['desc'], candidates=found[:5]))
        continue
    # journal type: try vol+page precise, then author+year+journal
    au = ref['au']; yr = ref['yr']; ta = ref['ta']; vol = ref['vol']; pg = ref['pg']; kw = ref['kw']
    queries = [
        f'{au}[Author] AND {ta}[TA] AND {vol}[VI] AND {pg}[PG]',
        f'{au}[Author] AND {yr}[PDAT] AND {ta}[TA]',
        f'{au}[Author] AND {ta}[TA] AND {kw}',
    ]
    picked = None; cands = []
    for q in queries:
        ids = esearch(q); time.sleep(0.4)
        arts = parse_articles(efetch(ids)); time.sleep(0.4)
        cands = arts
        for a in arts:
            # match if page start matches or vol matches + author matches
            pg_start = re.split(r'[-–]', a['page'])[0] if a['page'] else ''
            if a['vol'] == vol and (pg_start == pg or pg in a['page']):
                picked = a; break
        if picked: break
    if picked:
        cr = crossref(picked['doi']) if picked['doi'] else None; time.sleep(0.4)
        print(f"  MATCH pmid={picked['pmid']} {picked['first_author']} {picked['year']} {picked['journal']} {picked['vol']}:{picked['page']}")
        print(f"        title={picked['title'][:80]}")
        print(f"        DOI={picked['doi']}  crossref={'OK' if cr else 'n/a'}")
        results.append(dict(n=n, status='matched', doi=picked['doi'], pmid=picked['pmid'],
                            meta=picked, cr=cr))
    else:
        print(f"  NO CONFIDENT MATCH for [{n}] {au} {yr} {ta} {vol}:{pg}")
        for a in cands[:5]:
            print(f"    cand pmid={a['pmid']} {a['first_author']} {a['year']} {a['journal']} {a['vol']}:{a['page']} doi={a['doi']} | {a['title'][:60]}")
        results.append(dict(n=n, status='nomatch', doi='', candidates=cands[:5],
                            expected=dict(au=au, yr=yr, ta=ta, vol=vol, pg=pg)))

with open('doi_pubmed_results.json', 'w', encoding='utf-8') as f:
    json.dump(results, f, ensure_ascii=False, indent=2)

print("\n\n========== SUMMARY ==========")
for r in results:
    print(f"[{r['n']:>2}] {r['status']:<12} {r.get('doi','')}")
print(f"\nSaved -> doi_pubmed_results.json")
