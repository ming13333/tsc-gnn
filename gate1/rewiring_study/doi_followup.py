#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Targeted follow-up searches for unresolved / flagged references. Robust per-query."""
import urllib.request, urllib.parse, ssl, json, time, re

ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
HDR = {'User-Agent': 'tsc-gnn-ref-verify/1.0 (mailto:hwenwu321@gmail.com)'}

def get(url):
    for _ in range(3):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers=HDR), timeout=30, context=ctx).read().decode('utf-8','replace')
        except Exception as e:
            last=e; time.sleep(2)
    return f"__ERROR__:{last}"

def search(term, n=6):
    u='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=%d&term=%s'%(n,urllib.parse.quote(term))
    r=get(u)
    if r.startswith('__ERROR__'): return []
    try: ids=json.loads(r)['esearchresult'].get('idlist',[])
    except Exception: return []
    time.sleep(0.4)
    if not ids: return []
    xml=get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&id='+','.join(ids))
    time.sleep(0.4)
    out=[]
    for a in re.split(r'<PubmedArticle>', xml)[1:]:
        d={}
        m=re.search(r'<PMID[^>]*>(\d+)</PMID>',a); d['pmid']=m.group(1) if m else ''
        m=re.search(r'<ArticleTitle>(.*?)</ArticleTitle>',a,re.S); d['title']=re.sub(r'<[^>]+>','',m.group(1)).strip() if m else ''
        m=re.search(r'<Volume>(.*?)</Volume>',a,re.S); d['vol']=m.group(1).strip() if m else ''
        m=re.search(r'<MedlinePgn>(.*?)</MedlinePgn>',a,re.S); d['page']=m.group(1).strip() if m else ''
        m=re.search(r'<ISOAbbreviation>(.*?)</ISOAbbreviation>',a,re.S); d['journal']=m.group(1).strip() if m else ''
        m=re.search(r'<Year>(\d{4})</Year>',a); d['year']=m.group(1) if m else ''
        m=re.search(r'<ArticleId IdType="doi">(.*?)</ArticleId>',a,re.S); d['doi']=m.group(1).strip().lower() if m else ''
        m=re.search(r'<Author[^>]*>\s*<LastName>(.*?)</LastName>',a,re.S); d['au']=m.group(1).strip() if m else ''
        out.append(d)
    return out

QUERIES = {
 '1  Cramer stroke recovery mechanisms': ['Cramer[Author] AND stroke[Title] AND recovery[Title] AND 2008[PDAT]',
                                          'Cramer SC repairing brain stroke mechanisms recovery'],
 '3  Li 2021 scRNA ischemic stroke': ['Li[Author] AND 2021[PDAT] AND ischemi*[Title] AND single-cell[Title]',
                                       'single-cell ischemic stroke temporal glial neuronal 2021',
                                       'Acta Neuropathol Commun[TA] AND 2021[PDAT] AND stroke[Title]'],
 '25 Fancy myelin regeneration 2011': ['Fancy[Author] AND myelin[Title]',
                                        'Fancy SP 2011 myelin regeneration oligodendrocyte'],
 '39 Lotfollahi 2023 biol-informed DL': ['Lotfollahi[Author] AND 2023[PDAT] AND (program*[Title] OR informed[Title])',
                                          'Lotfollahi biologically informed deep learning gene programs 2023'],
 '50 Kartha 2023 Perturb-seq gene reg': ['Kartha[Author] AND 2023[PDAT] AND Perturb*[Title]',
                                          'functional inference gene regulation single-cell Perturb-seq 2023'],
 '18 GNN review methods applications': ['Zhou[Author] AND graph neural network[Title] AND review',
                                         'Wu[Author] comprehensive survey graph neural networks'],
 '37 Velickovic 2023 GNN review': ['Velickovic[Author] AND 2023[PDAT] AND graph[Title]',
                                    'Velickovic everything is connected graph neural networks 2023'],
 '42 PerturbExpress / Singh': ['Singh[Author] AND perturbation[Title] AND single-cell',
                                'perturbation prediction benchmark platform single cell'],
 '47 Bhatt statin stroke': ['Bhatt[Author] AND statin*[Title] AND 2017[PDAT]',
                             'statins stroke prevention Circulation 2017'],
 '48 Li HDAC neuro 2020': ['histone deacetylase inhibitor* neurological disease* 2020 review',
                            'HDAC inhibitor neurodegenerative Nature Reviews Drug Discovery'],
 '54 Weider oligodendrocyte Nat Commun': ['Weider[Author] AND oligodendrocyte',
                                           'Weider Nfat calcineurin oligodendrocyte myelination'],
}

for label, qs in QUERIES.items():
    print(f"\n===== {label} =====", flush=True)
    seen=set()
    for q in qs:
        for a in search(q):
            key=a['pmid']
            if key in seen: continue
            seen.add(key)
            print(f"  {a['au']:<14} {a['year']} {a['journal']:<22} {a['vol']}:{a['page']:<10} doi={a['doi'] or '-'}")
            print(f"     {a['title'][:95]}")
