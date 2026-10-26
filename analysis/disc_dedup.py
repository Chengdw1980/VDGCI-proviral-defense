import pandas as pd, numpy as np, json
from scipy.stats import binomtest
from mechclass import mech
D='/tmp/claude-0/-home-claude/1d89e018-43ce-5a67-b16c-37382686b1fc/scratchpad/in/DC Viral Genomic Fossils for Bacterial Evolution Work Project 09-11-2026/'
s0=pd.read_csv(D+'defense_systems_detail.csv'); s=s0.drop_duplicates(['accession','sys_id']).copy()
g=pd.read_csv(D+'defense_per_genome.csv')[['accession','genome_bp','prophage_bp']]
c=pd.read_csv(D+'cohort48_VDGCI_table.csv')
miss=c[~c.accession.isin(g.accession)][['accession','bp','prov_bp']].rename(columns={'bp':'genome_bp','prov_bp':'prophage_bp'})
g=pd.concat([g,miss]); assert len(g)==48
s['mech']=s.type.map(mech); s['inp']=s.in_prophage.astype(str)=='True'
pv=g.prophage_bp.sum()/g.genome_bp.sum()
res=dict(n_rows=len(s0),n_unique=len(s),dups=len(s0)-len(s),dups_in_prophage=int((s0.duplicated(['accession','sys_id'])&(s0.in_prophage.astype(str)=='True')).sum()),pv=pv)
rng=np.random.default_rng(20261024); accs=g.accession.values; G=g.set_index('accession')
def est(x):
    K=x.inp.sum(); N=len(x); return K,N,(K/N)/pv
def boot(x,R=4000):
    per=x.groupby('accession').inp.agg(['sum','count']); out=[]
    for _ in range(R):
        ix=rng.choice(accs,len(accs)); sub=G.loc[ix]; pk=per.reindex(ix).fillna(0)
        if pk['count'].sum()>0: out.append((pk['sum'].sum()/pk['count'].sum())/(sub.prophage_bp.sum()/sub.genome_bp.sum()))
    return np.array(out)
K,N,f=est(s); b=boot(s)
res['overall']=dict(K=int(K),N=int(N),pct=100*K/N,fold=f,ci=list(np.percentile(b,[2.5,97.5])),binom_p=binomtest(int(K),int(N),pv,alternative='greater').pvalue)
# mechanism + ratio with paired bootstrap
per={m:x.groupby('accession').inp.agg(['sum','count']) for m,x in s.groupby('mech')}
M={}
for m,x in s.groupby('mech'):
    k,n,ff=est(x); M[m]=dict(K=int(k),N=int(n),fold=ff)
ratios=[]; mb={m:[] for m in per}
for _ in range(4000):
    ix=rng.choice(accs,len(accs)); sub=G.loc[ix]; p=sub.prophage_bp.sum()/sub.genome_bp.sum(); fv={}
    for m,pp in per.items():
        pk=pp.reindex(ix).fillna(0)
        fv[m]=(pk['sum'].sum()/pk['count'].sum())/p if pk['count'].sum()>0 else np.nan; mb[m].append(fv[m])
    a=fv['Abortive infection / toxin-antitoxin']; d=fv['Direct nucleic-acid targeting']
    if d and d>0 and not np.isnan(a): ratios.append(a/d)
for m in M: M[m]['ci']=list(np.nanpercentile(mb[m],[2.5,97.5]))
res['mech']=M
res['ratio']=dict(est=M['Abortive infection / toxin-antitoxin']['fold']/M['Direct nucleic-acid targeting']['fold'],ci=list(np.percentile(ratios,[2.5,97.5])))
T=s.groupby('type').agg(n=('inp','size'),k=('inp','sum')).sort_values('n',ascending=False)
res['types']={t:[int(r.n),int(r.k)] for t,r in T.iterrows()}
json.dump(res,open('disc_dedup.json','w'),indent=1,default=float)
print(json.dumps({k:res[k] for k in ['n_rows','n_unique','dups','dups_in_prophage','pv','overall','mech','ratio']},indent=1,default=float))
print({t:v for t,v in res['types'].items() if t in ['RM_Type_I','RM_Type_II','RM_Type_IIG','RM_Type_III','Cas','McrBC','Other_Type_IV','SNIPE','MazEF','RosmerTA','PfiAT','AbiD','AbiAlpha','AbiU','Retron','Gabija','Dag','Lamassu-Fam','MspJI','CoCoNut','GmrSD_RM_Type_IV']})
