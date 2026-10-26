import os, glob, json, itertools, numpy as np, pandas as pd
B='/anvil/projects/x-bio260593/vdgci/HQ_5054'; P=B+'/per_genome'; O=B+'/results/defmap/'
coh=pd.read_csv(B+'/results/cohort_VDGCI.csv')
loc=pd.read_csv(B+'/results/extended/defense_systems_located.tsv',sep='\t')
mech=dict(zip(loc.sys_id,loc.mechanism_class)); typ=dict(zip(loc.sys_id,loc.type)); inv=dict(zip(loc.sys_id,loc.in_viral_dna.astype(str)=='True'))
FL=10000
per=[]; sysrows=[]
for acc in coh.accession:
    d=f'{P}/{acc}'
    hf=glob.glob(d+'/*_protein_headers.txt'); gf=glob.glob(d+'/*_defense_finder_genes.tsv'); vf=glob.glob(d+'/*_virus_summary.tsv')
    H={}; clen={}
    if hf:
        for l in open(hf[0]):
            if l.startswith('>'):
                x=l[1:].split(' # '); a,b=int(x[1]),int(x[2]); H[x[0]]=(a,b); c=x[0].rsplit('_',1)[0]; clen[c]=max(clen.get(c,0),b)
    prov=[]; wholec=set()
    if vf:
        v=pd.read_csv(vf[0],sep='\t')
        for sn in v.seq_name:
            if '|provirus_' in sn:
                c,r=sn.split('|provirus_'); a,b=map(int,r.split('_')); prov.append((c,a,b))
            else: wholec.add(sn)
    flank_bp=0
    for c,a,b in prov:
        flank_bp+=min(FL,a-1)+min(FL,max(0,clen.get(c,b)-b))
    n=dict(acc=acc,n_sys=0,loose=0,strict=0,straddle=0,flank=0)
    if gf:
        g=pd.read_csv(gf[0],sep='\t')
        for sid,gg in g.groupby('sys_id'):
            co=[H[h] for h in gg.hit_id if h in H]
            if not co: continue
            c=gg.hit_id.iloc[0].rsplit('_',1)[0]; s=min(min(x) for x in co); e=max(max(x) for x in co); m=(s+e)/2
            cat='elsewhere'
            if c in wholec: cat='inside'
            else:
                for pc,a,b in prov:
                    if pc==c and s>=a and e<=b: cat='inside'; break
                    if pc==c and e>=a and s<=b: cat='straddle'
                if cat=='elsewhere':
                    for pc,a,b in prov:
                        if pc==c and ((a-FL<=m<a) or (b<m<=b+FL)): cat='flank'; break
            n['n_sys']+=1; n['loose']+=cat in('inside','straddle'); n['strict']+=cat=='inside'; n['straddle']+=cat=='straddle'; n['flank']+=cat=='flank'
            sysrows.append((acc,sid,typ.get(sid,''),mech.get(sid,''),cat))
    n['flank_bp']=flank_bp; per.append(n)
D=pd.DataFrame(per).merge(coh[['accession','phylum','genome_bp','viral_bp']],left_on='acc',right_on='accession')
D.to_csv(O+'per_genome_strict_flank.tsv',sep='\t',index=False)
S=pd.DataFrame(sysrows,columns=['accession','sys_id','type','mech','category']); S.to_csv(O+'system_categories.tsv.gz',sep='\t',index=False)
rng=np.random.default_rng(20261024)
def fold(df,num,bp): return (df[num].sum()/df.n_sys.sum())/(df[bp].sum()/df.genome_bp.sum())
def boot(df,num,bp,R=1000):
    idx=np.arange(len(df)); v=[]
    a=df[num].values; ns=df.n_sys.values; b=df[bp].values; g=df.genome_bp.values
    for _ in range(R):
        i=rng.choice(idx,len(idx)); v.append((a[i].sum()/ns[i].sum())/(b[i].sum()/g[i].sum()))
    return [float(np.percentile(v,2.5)),float(np.percentile(v,97.5))]
res={}
for name,num,bp in [('loose_inside_or_straddle','loose','viral_bp'),('strict_fully_inside','strict','viral_bp'),('flank_10kb_outside','flank','flank_bp')]:
    res[name]=dict(n=int(D[num].sum()),fold=float(fold(D,num,bp)),ci=boot(D,num,bp),pct_seq=float(D[bp].sum()/D.genome_bp.sum()*100))
res['n_sys_total']=int(D.n_sys.sum()); res['straddle']=int(D.straddle.sum())
res['strict_by_mech']={}
SS=S.merge(D[['acc']],left_on='accession',right_on='acc')
for mclass,grp in S.groupby('mech'):
    c=grp.groupby('accession').category.agg(lambda x:(x=='inside').sum()).rename('strict'); t=grp.groupby('accession').size().rename('n_sys')
    dd=D[['acc','genome_bp','viral_bp']].set_index('acc').join(pd.concat([c,t],axis=1)).fillna(0)
    res['strict_by_mech'][mclass]=dict(n=int(dd.strict.sum()),fold=float(fold(dd,'strict','viral_bp')),ci=boot(dd,'strict','viral_bp',500))
# positional, inside-only
R=pd.read_csv(O+'sys_positions.tsv',sep='\t')
R=R[(R.sys_start>=R.prov_start)&(R.sys_end<=R.prov_end)].reset_index(drop=True)
nulls=[np.array([float(x) for x in str(s).split(';') if x not in('','nan')]) for s in R.null_rel]
keep=[len(n)>0 for n in nulls]; R=R[keep].reset_index(drop=True); nulls=[n for n in nulls if len(n)>0]
NP=5000
idx=[rng.integers(0,len(n),NP) for n in nulls]; samp=np.array([n[i] for n,i in zip(nulls,idx)])
side=R.integrase_side.values; hasI=pd.notna(R.integrase_side).values
orn=np.where((side=='left')[:,None],samp,1-samp)
def postest(mask,label):
    m=mask; out={'n':int(m.sum())}
    if m.sum()<5: return out
    e_obs=(R.edge.values[m]<0.1).mean(); e_null=(np.minimum(samp,1-samp)[m]<0.1).mean(0)
    out.update(edge10_obs=float(e_obs),edge10_null=float(e_null.mean()),p_edge=float(((e_null>=e_obs).sum()+1)/(NP+1)))
    mi=m&hasI
    if mi.sum()>=5:
        o_obs=R.orient.values[mi]; on=orn[mi]
        far_obs=(o_obs>=0.8).mean(); far_null=(on>=0.8).mean(0)
        near_obs=(o_obs<0.2).mean(); near_null=(on<0.2).mean(0)
        out.update(n_int=int(mi.sum()),far20_obs=float(far_obs),far20_null=float(far_null.mean()),p_far=float(((far_null>=far_obs).sum()+1)/(NP+1)),
                   near20_obs=float(near_obs),near20_null=float(near_null.mean()),p_near=float(((near_null>=near_obs).sum()+1)/(NP+1)),
                   mean_orient=float(o_obs.mean()),mean_orient_null=float(on.mean()))
    return out
pos={'all':postest(np.ones(len(R),bool),'all')}
for mc in R.mech.dropna().unique(): pos['mech:'+mc]=postest((R.mech==mc).values,mc)
top=R.phylum.value_counts(); 
for ph in top.index[:5]: pos['phylum:'+ph]=postest((R.phylum==ph).values,ph)
pos['phylum:other']=postest((~R.phylum.isin(top.index[:5])).values,'other')
for tier,lab in [(['Complete','High-quality'],'HQ'),(['Medium-quality'],'MQ'),(['Low-quality'],'LQ')]: pos['checkv:'+lab]=postest(R.checkv.isin(tier).values,lab)
res['position_inside_only']=pos
# orient bins inside only
bins=np.linspace(0,1,11); ob,_=np.histogram(R.orient.values[hasI],bins)
nb=np.zeros(10)
for n,sd,h in zip(nulls,side,hasI):
    if h: hh,_=np.histogram(n if sd=='left' else 1-n,bins); nb+=hh/len(n)
pd.DataFrame(dict(bin_lo=bins[:-1],obs=ob,null_expected=nb)).to_csv(O+'orient_bins_inside.tsv',sep='\t',index=False)
# islands: type co-occurrence within same provirus (inside-only)
pairs={}
for pv,g in R.groupby('provirus'):
    ts=sorted(set(g.type))
    for a,b in itertools.combinations(ts,2): pairs[(a,b)]=pairs.get((a,b),0)+1
res['island_proviruses']=int((R.groupby('provirus').size()>=2).sum()); res['island_systems']=int(R.groupby('provirus').size()[lambda x:x>=2].sum())
res['top_pairs']={f'{a}+{b}':v for (a,b),v in sorted(pairs.items(),key=lambda x:-x[1])[:12]}
json.dump(res,open(O+'model_stats.json','w'),indent=1)
print(json.dumps({k:res[k] for k in ['loose_inside_or_straddle','strict_fully_inside','flank_10kb_outside','n_sys_total','straddle']}))
print(json.dumps(res['position_inside_only']['all']))
