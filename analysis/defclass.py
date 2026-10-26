import json, os, glob, numpy as np, pandas as pd
from scipy.stats import binomtest, chi2
B='/anvil/projects/x-bio260593/vdgci/HQ_5054/'; O=B+'results/defmap/'; OUT=O+'class/'; os.makedirs(OUT,exist_ok=True)
GN=os.environ.get('SCRATCH','/anvil/scratch/x-dcheng3')+'/vdgci/HQ_5054/genomad_out/'
rng=np.random.default_rng(20261025); MIN=30
S=pd.read_csv(O+'audit_system_ledger.tsv.gz',sep='\t'); D=pd.read_csv(O+'audit_per_genome_bp.tsv',sep='\t')
coh=pd.read_csv(B+'results/cohort_VDGCI.csv'); cls=dict(zip(coh.accession,coh['class'])); phy=dict(zip(coh.accession,coh.phylum)); fam=dict(zip(coh.accession,coh.family))
S=S[~S.cat.astype(str).str.startswith('unresolved')].copy()
S['inv']=S.cat.isin(['inside_provirus','inside_viral_record','straddle'])
S['mc']=np.where(S.mech.astype(str).str.startswith('Abortive'),'Abi',np.where(S.mech.astype(str).str.startswith('Direct'),'Direct','Other'))
def bh(p):
    p=np.asarray(p,float); m=len(p); o=np.argsort(p); q=np.empty(m); prev=1.0
    for r,i in enumerate(o[::-1]): prev=min(prev,p[i]*m/(m-r)); q[i]=prev
    return q
# ---------- per-genome table ----------
g=S.groupby('accession').agg(n=('inv','size'),k=('inv','sum'))
for mc in ['Abi','Direct']:
    t=S[S.mc==mc].groupby('accession').agg(**{'n_'+mc:('inv','size'),'k_'+mc:('inv','sum')}); g=g.join(t)
g=D.set_index('accession').join(g).fillna(0)
g['cls']=[cls.get(a,'NA') for a in g.index]; g['phy']=[phy.get(a,'NA') for a in g.index]; g['fam']=[fam.get(a,'NA') for a in g.index]
VB='viral_bp_recomputed'
def fold(df,k='k',n='n'): return (df[k].sum()/max(df[n].sum(),1))/(df[VB].sum()/df.genome_bp.sum())
def boot(df,k='k',n='n',R=2000):
    a=df[k].values;N=df[n].values;v=df[VB].values;G=df.genome_bp.values;out=[]
    for _ in range(R):
        i=rng.integers(0,len(df),len(df)); nn=N[i].sum()
        if nn==0 or v[i].sum()==0: continue
        out.append((a[i].sum()/nn)/(v[i].sum()/G[i].sum()))
    return [float(np.percentile(out,2.5)),float(np.percentile(out,97.5))] if out else [np.nan,np.nan]
def ratio_boot(df,R=2000):
    out=[]; idx=np.arange(len(df)); A=df.k_Abi.values;NA=df.n_Abi.values;Dd=df.k_Direct.values;ND=df.n_Direct.values
    for _ in range(R):
        i=rng.choice(idx,len(idx))
        if NA[i].sum()==0 or ND[i].sum()==0 or Dd[i].sum()==0: continue
        out.append((A[i].sum()/NA[i].sum())/(Dd[i].sum()/ND[i].sum()))
    return [float(np.percentile(out,2.5)),float(np.percentile(out,97.5))] if len(out)>100 else [np.nan,np.nan]
def rank_rows(level,col):
    rows=[]
    for c,df in g.groupby(col):
        if len(df)<MIN or c in('NA',''): continue
        n=int(df.n.sum()); k=int(df.k.sum()); f=df[VB].sum()/df.genome_bp.sum()
        P=binomtest(k,n,f).pvalue if n>0 else np.nan
        r=dict(level=level,taxon=c,phylum=df.phy.mode()[0],n_genomes=len(df),n_genomes_with_provirus=int((df[VB]>0).sum()),
               systems=n,systems_in_provirus=k,viral_bp_frac=f,expected=n*f,fold=fold(df),fold_ci=boot(df),P_binom=P)
        if df.k_Abi.sum()>0 and df.k_Direct.sum()>0:
            r['fold_Abi']=fold(df,'k_Abi','n_Abi'); r['fold_Direct']=fold(df,'k_Direct','n_Direct')
            r['Abi_Direct_ratio']=r['fold_Abi']/r['fold_Direct']; r['ratio_ci']=ratio_boot(df)
        rows.append(r)
    T=pd.DataFrame(rows); T['q_BH']=bh(T.P_binom.fillna(1)); return T.sort_values('fold',ascending=False)
C=rank_rows('class','cls'); F=rank_rows('family','fam')
# ---------- heterogeneity across classes (genome-label permutation keeps within-genome clustering) ----------
def Gstat(df,col):
    a=df.groupby(col).agg(k=('k','sum'),n=('n','sum'),v=(VB,'sum'),G=('genome_bp','sum'))
    f=a.v/a.G; th=a.k.sum()/(a.n*f).sum(); E=th*a.n*f; O_=a.k; N=a.n
    with np.errstate(divide='ignore',invalid='ignore'):
        t1=np.where(O_>0,O_*np.log(O_/E),0); t2=np.where(N-O_>0,(N-O_)*np.log((N-O_)/(N-E)),0)
    return float(2*np.nansum(t1+t2)), len(a)
keep=g[g.cls.isin(C.taxon)].copy(); G0,kc=Gstat(keep,'cls'); perm=[]
lab=keep.cls.values.copy()
for _ in range(2000):
    keep['pc']=rng.permutation(lab); perm.append(Gstat(keep,'pc')[0])
het=dict(n_classes=kc,G=G0,df=kc-1,P_LRT_chi2=float(chi2.sf(G0,kc-1)),P_genome_permutation=float((np.sum(np.array(perm)>=G0)+1)/(len(perm)+1)),
         n_classes_CI_above_1=int(sum(c[0]>1 for c in C.fold_ci)),n_classes_CI_below_1=int(sum(c[1]<1 for c in C.fold_ci)))
# ---------- class x defense type ----------
S['cls']=S.accession.map(cls); rows=[]
for c in C.taxon:
    sub=S[S.cls==c]; f=float(C.loc[C.taxon==c,'viral_bp_frac'].iloc[0])
    for t,d in sub.groupby('type'):
        n=len(d); k=int(d.inv.sum())
        if n<10 or k<3: continue
        rows.append(dict(taxon=c,type=t,mech=d.mc.mode()[0],systems=n,in_provirus=k,expected=n*f,fold=(k/n)/f,genomes_with_type=d.accession.nunique(),P_binom=binomtest(k,n,f).pvalue))
TT=pd.DataFrame(rows)
if len(TT): TT['q_BH']=bh(TT.P_binom); TT=TT.sort_values(['q_BH','fold'],ascending=[True,False])
# ---------- positions by class (null B, span-matched) ----------
Pp=pd.read_csv(O+'audit_positions.tsv.gz',sep='\t'); Pp=Pp[Pp.nullB.notna()].reset_index(drop=True); Pp['cls']=Pp.accession.map(cls)
arr=lambda s:np.array([float(x) for x in str(s).split(';') if x not in('','nan')])
nB=[arr(x) for x in Pp.nullB]; side=Pp.anchor_side.values; hasI=pd.notna(Pp.anchor_side).values; rel=Pp.rel.values
orient=np.where(side=='left',rel,1-rel)
pe=np.array([(np.minimum(n,1-n)<0.1).mean() if len(n) else np.nan for n in nB])
pf=np.array([((n if sd=='left' else 1-n)>=0.8).mean() if len(n) and isinstance(sd,str) else np.nan for n,sd in zip(nB,side)])
oe=(np.minimum(rel,1-rel)<0.1).astype(float); of=(orient>=0.8).astype(float)
def ptest(mask,obs,pr):
    m=mask&~np.isnan(pr)
    if m.sum()<15: return None
    o=obs[m];p=pr[m];prov=Pp.provirus.values[m]; up=np.unique(prov); idx={u:np.where(prov==u)[0] for u in up}; bd=[]
    sims=(rng.random((10000,m.sum()))<p).mean(1)
    for _ in range(2000):
        ii=np.concatenate([idx[u] for u in rng.choice(up,len(up))]); bd.append(o[ii].mean()-p[ii].mean())
    return dict(n_sys=int(m.sum()),n_prov=len(up),obs=float(o.mean()),exp=float(p.mean()),diff=float(o.mean()-p.mean()),
                P_perm=float(((sims>=o.mean()).sum()+1)/10001),ci_cluster=[float(np.percentile(bd,2.5)),float(np.percentile(bd,97.5))])
PR=[]
for c in ['All']+list(C.taxon):
    mk=np.ones(len(Pp),bool) if c=='All' else (Pp.cls==c).values
    for metric,obs,pr,base in [('edge10',oe,pe,np.ones(len(Pp),bool)),('distal20',of,pf,hasI)]:
        t=ptest(mk&base,obs,pr)
        if t: PR.append(dict(taxon=c,metric=metric,**t))
PR=pd.DataFrame(PR)
if len(PR):
    PR['q_BH']=np.nan
    for mtr in PR.metric.unique():
        m=(PR.metric==mtr)&(PR.taxon!='All'); PR.loc[m,'q_BH']=bh(PR.loc[m,'P_perm'])
# ---------- tRNA-adjacent integration by class ----------
W=int(os.environ.get('TRNA_WIN','2000'))
defprov=set(S.loc[S.cat.isin(['inside_provirus','straddle']),'provirus'].dropna())
trows=[]
for acc in coh.accession:
    d=GN+acc+'/'+acc+'_find_proviruses/'
    pv=d+acc+'_provirus.tsv'; ar=d+acc+'_provirus_aragorn.tsv'
    if not os.path.exists(pv): continue
    try: P=pd.read_csv(pv,sep='\t')
    except Exception: continue
    if not len(P): continue
    tr=pd.read_csv(ar,sep='\t',header=None,names=['name','s','e']) if os.path.exists(ar) and os.path.getsize(ar)>0 else pd.DataFrame(columns=['name','s','e'])
    tr['contig']=tr.name.astype(str).str.rsplit('_tRNA',n=1).str[0]; tr['aa']=tr.name.astype(str).str.rsplit('_',n=1).str[-1]
    for _,r in P.iterrows():
        t=tr[tr.contig==r.source_seq]; hit=None
        for _,x in t.iterrows():
            for b in (r.start,r.end):
                if min(abs(x.s-b),abs(x.e-b))<=W or (x.s<=b<=x.e): hit=x.aa
        trows.append(dict(accession=acc,cls=cls.get(acc),provirus=r.seq_name,length=r.length,trna_adjacent=hit is not None,trna_aa=hit,
                          carries_defense=r.seq_name in defprov,integrase=isinstance(r.integrases,str) and len(r.integrases)>0))
TR=pd.DataFrame(trows); TR.to_csv(OUT+'provirus_trna_adjacency.tsv.gz',sep='\t',index=False)
trs=[]
for c,d in [('All',TR)]+[(c,TR[TR.cls==c]) for c in C.taxon]:
    if len(d)<20: continue
    a=d[d.carries_defense]; b=d[~d.carries_defense]
    top=d[d.trna_adjacent].trna_aa.value_counts().head(3)
    r=dict(taxon=c,proviruses=len(d),frac_trna_adjacent=float(d.trna_adjacent.mean()),top_trna=';'.join(f'{k}:{v}' for k,v in top.items()),
           defense_proviruses=len(a),frac_trna_def=float(a.trna_adjacent.mean()) if len(a) else np.nan,frac_trna_nodef=float(b.trna_adjacent.mean()) if len(b) else np.nan)
    if len(a)>=5 and len(b)>=5:
        from scipy.stats import fisher_exact
        r['P_fisher_def_vs_nodef']=float(fisher_exact([[a.trna_adjacent.sum(),(~a.trna_adjacent).sum()],[b.trna_adjacent.sum(),(~b.trna_adjacent).sum()]])[1])
    trs.append(r)
TS=pd.DataFrame(trs)
TS['P_chi2_across_classes']=np.nan
from scipy.stats import chi2_contingency
sub=TR[TR.cls.isin(C.taxon)]; ct=pd.crosstab(sub.cls,sub.trna_adjacent)
trna_het=dict(chi2_P=float(chi2_contingency(ct)[1]),window_bp=W,n_proviruses=len(TR))
# ---------- save ----------
C.to_csv(OUT+'class_enrichment.tsv',sep='\t',index=False); F.to_csv(OUT+'family_enrichment.tsv',sep='\t',index=False)
TT.to_csv(OUT+'class_type_enrichment.tsv',sep='\t',index=False); PR.to_csv(OUT+'class_positions.tsv',sep='\t',index=False); TS.to_csv(OUT+'class_trna.tsv',sep='\t',index=False)
json.dump(dict(heterogeneity=het,trna_heterogeneity=trna_het,overall_fold=fold(g),n_classes=len(C),n_families=len(F)),open(OUT+'class_summary.json','w'),indent=1,default=float)
pd.set_option('display.width',250); pd.set_option('display.max_columns',30); pd.set_option('display.max_rows',200)
print(json.dumps(dict(heterogeneity=het,trna=trna_het),default=float))
print(C.round(4).to_string()); print(F.round(4).to_string()); print(TT.head(40).round(4).to_string()); print(PR.round(4).to_string()); print(TS.round(3).to_string())
