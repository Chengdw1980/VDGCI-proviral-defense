import pandas as pd, numpy as np, gzip, treeswift, re
from physig import vcv, blomberg_k, pagel_lambda
S='/tmp/claude-0/-home-claude/1d89e018-43ce-5a67-b16c-37382686b1fc/scratchpad'
D=S+'/in/DC Viral Genomic Fossils for Bacterial Evolution Work Project 09-11-2026/'
tax=pd.read_csv(S+'/gtdb/bac120_taxonomy_r232.tsv.gz',sep='\t',header=None,names=['gid','tax'])
nwk=gzip.open(S+'/gtdb/bac120_r232.tree.gz','rt').read()
tips=set(re.findall(r'((?:RS|GB)_GC[AF]_\d+\.\d+)',nwk))
tax['sp']=tax.tax.str.extract(r's__(.*)$')
rep=tax[tax.gid.isin(tips)].set_index('sp').gid.to_dict()
h=pd.read_csv('gtdb_tax_48.tsv',sep='\t'); h['core']=h['core'].astype(str).str.zfill(9)
v=pd.read_csv(D+'cohort48_VDGCI_table.csv'); d=pd.read_csv(D+'defense_per_genome.csv')
sd=pd.read_csv(D+'defense_systems_detail.csv').drop_duplicates(['accession','sys_id'])
d['n_defense_systems']=d.accession.map(sd.groupby('accession').size()).fillna(0).astype(int)
v['core']=v.accession.str[4:].str.split('.').str[0]
v=v.merge(h[['core','species']].rename(columns={'species':'gtdb_sp'}),on='core',how='left')
v['tip']=v.gtdb_sp.map(rep)
print('mapped tips',v.tip.notna().sum(),'unique',v.tip.nunique())
v=v.merge(d[['accession','n_defense_systems','n_defense_systems_in_prophage']],on='accession',how='left')
v=v[v.tip.notna()].drop_duplicates('tip')
t=treeswift.read_tree_newick(nwk)
tl=list(v.tip); C=vcv(t,tl)
C=C/np.max(np.diag(C))
rows=[]
for col,lab,tr in [('pct_viral','% viral-derived DNA',lambda x:np.log1p(x)),('n_prov','Provirus count',lambda x:np.log1p(x)),
                   ('VDGCI','VDGCI',lambda x:x),('n_defense_systems','Defense systems (count)',lambda x:np.log1p(x)),
                   ('bp','Genome size (positive control)',lambda x:np.log(x))]:
    sub=v[v[col].notna()]; idx=[tl.index(x) for x in sub.tip]; Cs=C[np.ix_(idx,idx)]
    x=tr(sub[col].values.astype(float))
    K,pK=blomberg_k(x,Cs,nperm=999); lam,lr,pl=pagel_lambda(x,Cs)
    rows.append([lab,len(sub),round(K,3),pK,round(lam,3),round(lr,2),pl])
out=pd.DataFrame(rows,columns=['trait','n','Blomberg_K','P_K','Pagel_lambda','LR','P_lambda'])
print(out.to_string()); out.to_csv('../dedup/physig48_gtdb_dedup.csv',index=False)

