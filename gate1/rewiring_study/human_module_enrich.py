"""
Step C(2) — 人 TF 模块富集 / 跨物种模块一致性（本地、无网络依赖 · 逐 TF 版）
================================================================
修订：初版把 12 个广泛主调控因子的靶模块合并 → |H|=12106(占全集65%) 过宽，超几何被稀释、
置换经验 p≈0.5（伪阴性）。改为**逐 TF** 富集 + **靶基因数匹配**置换对照（每个真实 TF 与
靶数相近的随机人 TF 比较），才能识别特异保守链接（如 Sox10→髓鞘, Cebpb→炎症）。

设计：
1. M = 鼠三套分析 Top20(|ΔW|max) 中 >=2 套出现的 TF（Step C(1)）。
2. 正交 → 人同源 TF（symbol 大写）。
3. 逐 TF：其人 DoRothEA 靶集 Ht vs 文献策展参考集 R（myelin/oligo, neuroinflammation）。
4. 置换 null：从人 DoRothEA source 中抽「靶数匹配」(0.5–2x) 的随机 TF，2000 次，算经验 p。
5. 诚实边界：human DoRothEA 编码已知调控 → 富集证明「鼠识别的调控程序在人 GRN 保守」；
   置换排除「任意广靶 TF 都富集」的平庸解释。鼠 M 来自鼠数据独立识别，人验证是分离知识库。
"""
import pandas as pd, numpy as np, json, random
from scipy.stats import hypergeom
from collections import defaultdict, Counter

random.seed(42); np.random.seed(42)
WD = "."

full = pd.read_csv(f"{WD}/rewiring_full.csv")
c1   = pd.read_csv(f"{WD}/cc_gse174_rewiring_full.csv")
c2   = pd.read_csv(f"{WD}/cc_gse225_rewiring_full.csv")
FULL_TRANS=[("sham","24h"),("24h","2d"),("2d","14d"),("sham","14d")]
C1_TRANS=[("sham","24h")]; C2_TRANS=[("sham","2d"),("2d","14d")]

def tf_maxdw(df, trans):
    best={}
    for (a,b) in trans:
        col=f"dW_{a}_{b}"
        if col not in df.columns: continue
        g=df.groupby("tf")[col].apply(lambda s:s.abs().max())
        for tf,v in g.items(): best[tf]=max(best.get(tf,0.0),v)
    return pd.Series(best)

s_full,s_c1,s_c2=tf_maxdw(full,FULL_TRANS),tf_maxdw(c1,C1_TRANS),tf_maxdw(c2,C2_TRANS)
top_full,top_c1,top_c2=(set(s.sort_values(ascending=False).head(20).index) for s in (s_full,s_c1,s_c2))
cnt=Counter()
for t in top_full: cnt[t]+=1
for t in top_c1:   cnt[t]+=1
for t in top_c2:   cnt[t]+=1
M=sorted([tf for tf,c in cnt.items() if c>=2])
print(f"[M] 鼠可重现核心 TF 模块 (>=2/3 套 Top20), n={len(M)}: {M}")
print(f"   Sox10 in all3 Top20? {'Sox10' in top_full and 'Sox10' in top_c1 and 'Sox10' in top_c2}")

manual={'Sox10':'SOX10','Sox2':'SOX2','Sox9':'SOX9','Cebpb':'CEBPB','Erg':'ERG','Gata2':'GATA2',
 'Gata3':'GATA3','Nr2f2':'NR2F2','Pax5':'PAX5','Runx3':'RUNX3','E2f1':'E2F1','Ar':'AR'}
HUM_M=[manual.get(tf,tf.upper()) for tf in M]

HUMAN=f"{WD}/../data/dorothea/human_dorothea_regulon.tsv"
hd=pd.read_csv(HUMAN,sep="\t")
src2tgt=defaultdict(set)
for s,t in zip(hd['source'],hd['target']): src2tgt[s].add(t)
all_targets=set(hd['target'].unique()); N=len(all_targets)
all_sources=list(src2tgt.keys())
src_nt={s:len(v) for s,v in src2tgt.items()}
print(f"[human DoRothEA] sources={len(all_sources)} targets(universe)={N}")

HUM_M_present=[t for t in HUM_M if t in src2tgt]
missing=[t for t in HUM_M if t not in src2tgt]
print(f"[orthology] 人同源(在 DoRothEA source): {HUM_M_present}")
if missing: print(f"   缺失: {missing}")

REF={
 'myelin_oligodendrocyte': set(['PLP1','MBP','MOG','CNP','MAG','MOBP','OLIG1','OLIG2','CLDN11','MAL',
   'MPZ','PMP22','MYRF','NKX6-2','OPALIN','TSPAN2','UGT8','GJC2','KLK6','ASPA','NKX2-2','PDGFRA','ST18']),
 'neuroinflammation': set(['CEBPB','CEBPA','IL1B','IL6','TNF','NFKB1','NFKBIA','STAT1','STAT3','CCL2','CCL3',
   'CXCL10','CXCL1','TLR2','TLR4','IRF1','IRF7','NLRP3','PTGS2','NOS2','HMGB1','RELA','CXCL2','CCL5']),
}

def enrich(Ht, Rset):
    n=len(Rset & all_targets); k=len(Ht & Rset); Kh=len(Ht)
    if n==0 or k==0: return None
    p=hypergeom.sf(k-1,N,n,Kh)
    OR=(k/(Kh-k+1e-9))/((n-k)/(N-n-(Kh-k)+1e-9)+1e-9)
    return dict(k=k,n=n,Kh=Kh,p=p,OR=OR,overlap=sorted(Ht & Rset))

N_PERM=2000
print(f"\n=== 逐 TF 富集（人靶集 vs 参考集，{N_PERM}x 靶数匹配置换）===")
per_tf={}
for t in HUM_M_present:
    Ht=src2tgt[t]; nt=len(Ht)
    # 靶数匹配 null 候选
    cand=[s for s in all_sources if s!=t and 0.5*nt<=src_nt[s]<=2.0*nt]
    if len(cand)<50:  # 放宽
        cand=[s for s in all_sources if s!=t and 0.25*nt<=src_nt[s]<=4.0*nt]
    row={'nt':nt,'refs':{}}
    print(f"\n  ▸ {t} (人靶数={nt}, null候选={len(cand)})")
    for name,R in REF.items():
        Rr=R-set(HUM_M_present)  # 去除 M 同源自身
        r=enrich(Ht,Rr)
        if r is None:
            row['refs'][name]=dict(k=0,n=len(Rr & all_targets),OR=None,p=None,emp_p=None,
                                   rand_mean=None,rand_std=None,overlap=[])
            print(f"     {name:24s} 无重叠 (k=0)")
            continue
        # 置换
        emp=0; ks=[]
        for _ in range(N_PERM):
            rs=random.choice(cand); Hr=src2tgt[rs]
            ks.append(len(Hr & Rr))
            rr=enrich(Hr,Rr)
            if rr is not None and rr['p']<=r['p']: emp+=1
        ep=(1+emp)/(1+N_PERM)
        ks=np.array(ks)
        row['refs'][name]=dict(k=r['k'],n=r['n'],OR=r['OR'],p=r['p'],emp_p=ep,
                               rand_mean=float(ks.mean()),rand_std=float(ks.std()),
                               overlap=r['overlap'])
        flag="*** 显著" if ep<0.05 else ""
        print(f"     {name:24s} k={r['k']:2d}/{r['n']:2d} OR={r['OR']:.2f} p={r['p']:.2e} 经验p={ep:.4f} "
              f"(随机k̄={ks.mean():.1f}±{ks.std():.1f}) {flag}")
    per_tf[t]=row

out=dict(M=M, HUM_M_present=HUM_M_present, missing=missing, N_universe=N, per_tf=per_tf)
with open("human_module_enrich.json","w") as f: json.dump(out,f,indent=2,default=str)
print("\n[DONE] -> human_module_enrich.json")
