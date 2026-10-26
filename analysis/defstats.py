import json, numpy as np, pandas as pd
O='/anvil/projects/x-bio260593/vdgci/HQ_5054/results/defmap/'
R=pd.read_csv(O+'sys_positions.tsv',sep='\t')
G=pd.read_csv(O+'defense_provirus_genes.tsv.gz',sep='\t')
rng=np.random.default_rng(20261024)
nulls=[np.array([float(x) for x in str(s).split(';') if x not in ('','nan')]) for s in R.null_rel]
ok=np.array([len(n)>0 for n in nulls]); R=R[ok].reset_index(drop=True); nulls=[n for n in nulls if len(n)>0]
side=R.integrase_side.values
def orient(v,sd): return v if sd=='left' else 1-v
obs_edge=R.edge.values; hasI=pd.notna(R.integrase_side).values
obs_or=R.orient.values
NP=10000
me=np.empty(NP); f10=np.empty(NP); mo=np.empty(NP); h1=np.empty(NP)
idx=[rng.integers(0,len(n),NP) for n in nulls]
samp=np.array([n[i] for n,i in zip(nulls,idx)])  # systems x NP
edge_null=np.minimum(samp,1-samp)
or_null=np.where((side=='left')[:,None],samp,1-samp)[hasI]
me=edge_null.mean(0); f10=(edge_null<0.1).mean(0); mo=or_null.mean(0); h1=(or_null<0.5).mean(0)
res=dict(n_systems=int(len(R)),n_proviruses=int(R.provirus.nunique()),n_with_integrase=int(hasI.sum()),
 obs_mean_edge=float(obs_edge.mean()),null_mean_edge=float(me.mean()),p_edge_lower=float(((me<=obs_edge.mean()).sum()+1)/(NP+1)),p_edge_higher=float(((me>=obs_edge.mean()).sum()+1)/(NP+1)),
 obs_frac_edge10=float((obs_edge<0.1).mean()),null_frac_edge10=float(f10.mean()),null_f10_ci=[float(np.percentile(f10,2.5)),float(np.percentile(f10,97.5))],
 obs_mean_orient=float(np.nanmean(obs_or[hasI])),null_mean_orient=float(mo.mean()),p_or_lower=float(((mo<=np.nanmean(obs_or[hasI])).sum()+1)/(NP+1)),p_or_higher=float(((mo>=np.nanmean(obs_or[hasI])).sum()+1)/(NP+1)),
 obs_frac_integrase_half=float((obs_or[hasI]<0.5).mean()),null_frac_integrase_half=float(h1.mean()),null_h_ci=[float(np.percentile(h1,2.5)),float(np.percentile(h1,97.5))])
# multi-system proviruses (defense islands)
c=R.groupby('provirus').size(); res['prov_with_ge2_sys']=int((c>=2).sum()); res['sys_in_multi']=int(c[c>=2].sum())
# binned densities (10 bins), null weighted so each system counts 1
bins=np.linspace(0,1,11)
ob,_=np.histogram(obs_or[hasI],bins)
nb=np.zeros(10)
for n,sd,h in zip(nulls,side,hasI):
    if not h: continue
    v=orient(n,sd); hh,_=np.histogram(v,bins); nb+=hh/len(n)
eb,_=np.histogram(obs_edge,np.linspace(0,.5,6)); enb=np.zeros(5)
for n in nulls:
    hh,_=np.histogram(np.minimum(n,1-n),np.linspace(0,.5,6)); enb+=hh/len(n)
pd.DataFrame(dict(bin_lo=bins[:-1],obs=ob,null_expected=nb)).to_csv(O+'orient_bins.tsv',sep='\t',index=False)
pd.DataFrame(dict(bin_lo=np.linspace(0,.5,6)[:-1],obs=eb,null_expected=enb)).to_csv(O+'edge_bins.tsv',sep='\t',index=False)
# neighbour classes: nearest annotated non-defense gene on each side
nb_cls=[]
for _,s in R.iterrows():
    g=G[(G.prov==s.provirus)&(~G.is_def)&(G.cls!='unannotated')]
    L=g[g.end<s.sys_start].sort_values('end'); Rr=g[g.start>s.sys_end].sort_values('start')
    nb_cls += ([L.cls.iloc[-1]] if len(L) else []) + ([Rr.cls.iloc[0]] if len(Rr) else [])
res['neighbour_classes']=pd.Series(nb_cls).value_counts().to_dict()
gc=G[~G.is_def].cls.value_counts().to_dict(); res['all_nondef_gene_classes']=gc
# per type
T=R.groupby('type').agg(n=('edge','size'),mean_edge=('edge','mean'),mean_orient=('orient','mean')).sort_values('n',ascending=False)
T.head(20).to_csv(O+'per_type_positions.tsv',sep='\t')
res['checkv_of_def_proviruses']=R.drop_duplicates('provirus').checkv.value_counts().to_dict()
json.dump(res,open(O+'position_stats.json','w'),indent=1)
# slim table
S=R[['accession','phylum','type','mech','checkv','prov_len','n_prov_genes','n_struct','n_lysis','n_integ','n_def_genes_in_prov','edge','orient','integrase_side','provirus','sys_start','sys_end']]
S.to_csv(O+'sys_positions_slim.tsv',sep='\t',index=False)
# example shortlist
ex=[]
for t in ['Lamassu-Fam','AbiD','RM_Type_II','CBASS','Thoeris','Aristaios','RM_Type_I','Gabija','AbiE','Septu']:
    c=R[(R.type==t)&R.checkv.isin(['Complete','High-quality'])&(R.n_integ>=1)&(R.n_struct>=5)&(R.n_lysis>=1)]
    c=c.sort_values(['n_struct','n_lysis'],ascending=False).head(3); ex.append(c)
E=pd.concat(ex); res['n_example_candidates']=len(E)
EG=G[G.prov.isin(E.provirus)].copy(); EG['annotation_description']=EG.annotation_description.astype(str).str.slice(0,45)
EG.drop(columns=['gene']).to_csv(O+'example_genes.tsv',sep='\t',index=False)
E[['accession','phylum','type','checkv','provirus','prov_len','sys_start','sys_end','n_struct','n_lysis']].to_csv(O+'example_list.tsv',sep='\t',index=False)
print(json.dumps(res,indent=0)[:1500]); print('example genes',len(EG))
