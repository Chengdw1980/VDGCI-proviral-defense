import json, numpy as np, pandas as pd
O='/anvil/projects/x-bio260593/vdgci/HQ_5054/results/defmap/'
led=json.load(open(O+'audit_ledger.json'))
S=pd.read_csv(O+'audit_system_ledger.tsv.gz',sep='\t'); D=pd.read_csv(O+'audit_per_genome_bp.tsv',sep='\t')
coh=pd.read_csv('/anvil/projects/x-bio260593/vdgci/HQ_5054/results/cohort_VDGCI.csv')
rng=np.random.default_rng(20261024); res={'ledger':led}
# enrichment with exact bp (resolved systems only)
res_sys=S[~S.cat.str.startswith('unresolved')]
g=res_sys.groupby('accession').cat.value_counts().unstack(fill_value=0)
g['n']=g.sum(1); g=D.set_index('accession').join(g).fillna(0)
for c in ['inside_provirus','inside_viral_record','straddle','flank']:
    if c not in g: g[c]=0
def fold(df,num,bp): return (df[num].sum()/df.n.sum())/(df[bp].sum()/df.genome_bp.sum())
def boot(df,num,bp,R=1000):
    a=df[num].values;n=df.n.values;b=df[bp].values;G=df.genome_bp.values;v=[]
    for _ in range(R):
        i=rng.integers(0,len(df),len(df)); v.append((a[i].sum()/n[i].sum())/(b[i].sum()/G[i].sum()))
    return [float(np.percentile(v,2.5)),float(np.percentile(v,97.5))]
g['inside']=g.inside_provirus+g.inside_viral_record; g['loose']=g.inside+g.straddle
E={}
for name,num,bp in [('inside','inside','viral_bp_recomputed'),('loose','loose','viral_bp_recomputed'),('flank_exact','flank','flank_bp')]:
    E[name]=dict(n=int(g[num].sum()),fold=float(fold(g,num,bp)),ci=boot(g,num,bp),bp=int(g[bp].sum()))
# side-specific flank
fs=res_sys[res_sys.cat=='flank'].flank_side.value_counts().to_dict(); E['flank_side_counts']=fs
tot=g.n.sum(); gb=g.genome_bp.sum()
for s in ['anchor','far']:
    E['flank_'+s+'_fold']=float((fs.get(s,0)/tot)/(led['bp']['flank_bp_'+s+'_side']/gb))
res['enrichment_exact']=E
# positional tests
Pp=pd.read_csv(O+'audit_positions.tsv.gz',sep='\t')
Pp=Pp[Pp.nullA.notna()].reset_index(drop=True)
def arr(s): return np.array([float(x) for x in str(s).split(';') if x not in('','nan')])
nA=[arr(x) for x in Pp.nullA]; nB=[arr(x) for x in Pp.nullB]
side=Pp.anchor_side.values; hasI=pd.notna(Pp.anchor_side).values
rel=Pp.rel.values; orient=np.where(side=='left',rel,1-rel)
def probs(nulls,fn):
    return np.array([fn(n).mean() if len(n) else np.nan for n in nulls])
edge=lambda v:(np.minimum(v,1-v)<0.1)
pe_A=probs(nA,edge); pe_B=probs(nB,edge)
far=[];nearA=[]
pf_A=np.array([((n if sd=='left' else 1-n)>=0.8).mean() if len(n) and isinstance(sd,str) else np.nan for n,sd in zip(nA,side)])
pf_B=np.array([((n if sd=='left' else 1-n)>=0.8).mean() if len(n) and isinstance(sd,str) else np.nan for n,sd in zip(nB,side)])
oe=(np.minimum(rel,1-rel)<0.1).astype(float); of=(orient>=0.8).astype(float)
Pp['phy']=np.where(Pp.phylum.isin(['Bacillota','Pseudomonadota']),Pp.phylum,'Other phyla')
cv=pd.read_csv(O+'sys_positions.tsv',sep='\t')[['provirus','checkv']].drop_duplicates('provirus')
Pp=Pp.merge(cv,on='provirus',how='left')
def test(mask,obs,pr,label):
    m=mask&~np.isnan(pr)
    if m.sum()<5: return None
    o=obs[m]; p=pr[m]; prov=Pp.provirus.values[m]
    diff=o.mean()-p.mean()
    sims=(rng.random((10000,m.sum()))<p).mean(1); P_sys=float(((sims>=o.mean()).sum()+1)/10001)
    # cluster bootstrap over proviruses
    up=np.unique(prov); idx={u:np.where(prov==u)[0] for u in up}; bd=[]
    for _ in range(2000):
        pick=rng.choice(up,len(up)); ii=np.concatenate([idx[u] for u in pick]); bd.append(o[ii].mean()-p[ii].mean())
    bd=np.array(bd)
    # cluster sign test: one random system per provirus, repeated
    one=[]
    for _ in range(500):
        ii=np.array([rng.choice(idx[u]) for u in up]); one.append(o[ii].mean()-p[ii].mean())
    return dict(label=label,n_sys=int(m.sum()),n_prov=int(len(up)),obs=float(o.mean()),exp=float(p.mean()),diff=float(diff),P_system_perm=P_sys,
                diff_ci_cluster=[float(np.percentile(bd,2.5)),float(np.percentile(bd,97.5))],P_cluster_boot=float((bd<=0).mean()),one_per_prov_mean_diff=float(np.mean(one)))
groups=[('All',np.ones(len(Pp),bool))]+[('mech:'+m,(Pp.mech==m).values) for m in Pp.mech.dropna().unique()]+[('phy:'+p,(Pp.phy==p).values) for p in ['Bacillota','Pseudomonadota','Other phyla']]+\
 [('checkv:HQ',Pp.checkv.isin(['Complete','High-quality']).values),('checkv:MQ',(Pp.checkv=='Medium-quality').values),('checkv:LQ',(Pp.checkv=='Low-quality').values)]
R={}
for metric,obs,pA,pB,maskx in [('edge10',oe,pe_A,pe_B,np.ones(len(Pp),bool)),('far20',of,pf_A,pf_B,hasI)]:
    for nullname,pr in [('nullA_gene',pA),('nullB_spanmatched',pB)]:
        out=[]
        for lab,mk in groups:
            t=test(mk&maskx,obs,pr,lab)
            if t: out.append(t)
        ps=np.array([t['P_system_perm'] for t in out[1:]]); order=np.argsort(ps); q=np.empty(len(ps)); m=len(ps)
        prev=1
        for r,i in enumerate(order[::-1]):
            rank=m-r; prev=min(prev,ps[i]*m/rank); q[i]=prev
        for t,qq in zip(out[1:],q): t['q_BH_subgroups']=float(qq)
        R[f'{metric}|{nullname}']=out
res['positions']=R
# interactions (cluster bootstrap of difference in diffs), null B
def did(maskA,maskB,obs,pr,base):
    a=maskA&base&~np.isnan(pr); b=maskB&base&~np.isnan(pr)
    prov=Pp.provirus.values; outs=[]
    ua=np.unique(prov[a]); ub=np.unique(prov[b])
    ia={u:np.where(a&(prov==u))[0] for u in ua}; ib={u:np.where(b&(prov==u))[0] for u in ub}
    obsd=(obs[a].mean()-pr[a].mean())-(obs[b].mean()-pr[b].mean())
    for _ in range(2000):
        x=np.concatenate([ia[u] for u in rng.choice(ua,len(ua))]); y=np.concatenate([ib[u] for u in rng.choice(ub,len(ub))])
        outs.append((obs[x].mean()-pr[x].mean())-(obs[y].mean()-pr[y].mean()))
    outs=np.array(outs); return dict(diff_in_diff=float(obsd),ci=[float(np.percentile(outs,2.5)),float(np.percentile(outs,97.5))],P_two_sided=float(2*min((outs<=0).mean(),(outs>=0).mean())))
bp=Pp.phy.isin(['Bacillota','Pseudomonadota']).values
res['interactions_far20_nullB']={
 'Bacillota+Pseudomonadota vs other phyla':did(bp,~bp,of,pf_B,hasI),
 'HQ+MQ vs LQ':did(Pp.checkv.isin(['Complete','High-quality','Medium-quality']).values,(Pp.checkv=='Low-quality').values,of,pf_B,hasI),
 'AbiTA vs Direct':did((Pp.mech=='Abortive infection / toxin-antitoxin').values,(Pp.mech=='Direct nucleic-acid targeting').values,of,pf_B,hasI)}
res['position_definition']='system position = midpoint of span from first-gene start to last-gene end (min/max coordinates of constituent DefenseFinder genes), relative to geNomad provirus start/end; strand ignored; anchor = integration-annotated gene whose midpoint is closest to either end (first in coordinate order on ties); proviruses with anchors at both ends are oriented by the closer one'
json.dump(res,open(O+'robust_stats.json','w'),indent=1)
print(json.dumps(res['enrichment_exact'],indent=0)); print(json.dumps({k:[(t['label'],round(t['obs'],3),round(t['exp'],3),round(t['P_system_perm'],4),[round(x,3) for x in t['diff_ci_cluster']]) for t in v] for k,v in R.items()})[:3000])
print(json.dumps(res['interactions_far20_nullB']))
