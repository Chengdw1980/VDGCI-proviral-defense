import pandas as pd, numpy as np, glob, os, json
from scipy.stats import binomtest
D='/tmp/claude-0/-home-claude/1d89e018-43ce-5a67-b16c-37382686b1fc/scratchpad/in/DC Viral Genomic Fossils for Bacterial Evolution Work Project 09-11-2026/'
acc=[l.strip() for l in open('acc.txt') if l.strip()]
pv=pd.read_csv(D+'cohort48_provirus_detail.csv'); cs=pd.read_csv(D+'consensus_proviruses.tsv',sep='\t')
orig=pd.read_csv(D+'defense_systems_detail.csv').drop_duplicates(['accession','sys_id'])
L={}; 
for a in acc:
    cl={}; cur=None
    for l in open(f'g/{a}.fna'):
        if l[0]=='>': cur=l[1:].split()[0]; cl[cur]=0
        else: cl[cur]+=len(l.strip())
    L[a]=cl
def ivs(a):
    out=[]
    for _,r in pv[pv.accession==a].iterrows():
        if isinstance(r.coordinates,str):
            c=r.element.split('|')[0]; s,e=map(int,r.coordinates.split('-')); out.append((c,s,e))
        else: out.append((r.element,1,L[a].get(r.element,r.length_bp)))
    return out
def civs(a):
    out=[]
    for _,r in cs[(cs.accession==a)&(cs.phispy_overlap_frac>=0.5)].iterrows():
        out.append((r.contig.split('__',1)[1],int(r.start),int(r.end)))
    return out
rows=[]; G=[]; missing=[]
for a in acc:
    gbp=sum(L[a].values()); P=ivs(a); C=civs(a)
    pbp=sum(e-s+1 for _,s,e in P); cbp=sum(e-s+1 for _,s,e in C)
    sf=glob.glob(f'df/{a}/*_systems.tsv'); gf=glob.glob(f'df/{a}/*_genes.tsv')
    if not sf: missing.append(a); continue
    H={}
    for l in open(f'faa/{a}.faa'):
        if l[0]=='>':
            x=l[1:].split(' # '); H[x[0]]=(x[0].rsplit('_',1)[0],int(x[1]),int(x[2]))
    s=pd.read_csv(sf[0],sep='\t').drop_duplicates('sys_id'); g=pd.read_csv(gf[0],sep='\t')
    g=g[g.sys_id.isin(s.sys_id)]
    def inside(hits,iv):
        for h in hits:
            c,x,y=H[h]
            for cc,s_,e_ in iv:
                if cc==c and y>=s_ and x<=e_: return True
        return False
    for _,r in s.iterrows():
        hits=g[g.sys_id==r.sys_id].hit_id.tolist()
        rows.append(dict(accession=a,sys_id=r.sys_id,type=r.type,inp=inside(hits,P),incons=inside(hits,C)))
    ug=g.drop_duplicates('hit_id').hit_id
    gin=sum(1 for h in ug if any(cc==H[h][0] and H[h][2]>=s_ and H[h][1]<=e_ for cc,s_,e_ in P))
    G.append(dict(accession=a,genome_bp=gbp,prov_bp=pbp,cons_bp=cbp,n_genes=len(ug),genes_in=gin,n_sys=len(s)))
R=pd.DataFrame(rows); G=pd.DataFrame(G); R.to_csv('recount_systems.csv',index=False); G.to_csv('recount_genomes.csv',index=False)
res={'missing':missing,'n_genomes':len(G)}
N=len(R); K=int(R.inp.sum()); p=G.prov_bp.sum()/G.genome_bp.sum()
res['overall']=dict(N=N,K=K,pct=100*K/N,seqfrac=100*p,fold=(K/N)/p,P=binomtest(K,N,p,alternative='greater').pvalue)
pc=G.cons_bp.sum()/G.genome_bp.sum(); Kc=int(R.incons.sum())
res['consensus']=dict(K=Kc,N=N,seqfrac=100*pc,fold=(Kc/N)/pc,P=binomtest(Kc,N,pc,alternative='greater').pvalue)
gin=G.genes_in.sum(); gout=G.n_genes.sum()-gin
res['density']=dict(genes_in=int(gin),genes_out=int(gout),per_Mb_in=gin/(G.prov_bp.sum()/1e6),per_Mb_out=gout/((G.genome_bp.sum()-G.prov_bp.sum())/1e6))
res['density']['ratio']=res['density']['per_Mb_in']/res['density']['per_Mb_out']
# reproduction check vs original (deduplicated) records
o=orig.groupby(['accession','type']).size().rename('orig'); n=R.groupby(['accession','type']).size().rename('new')
cmp=pd.concat([o,n],axis=1).fillna(0); res['repro']=dict(orig_total=int(o.sum()),new_total=int(n.sum()),type_genome_pairs_identical=int((cmp.orig==cmp.new).sum()),pairs_total=len(cmp),abs_diff=int((cmp.orig-cmp.new).abs().sum()),orig_in=int((orig.in_prophage.astype(str)=='True').sum()))
json.dump(res,open('recount.json','w'),indent=1,default=float); print(json.dumps(res,indent=1,default=float))
