# Audit + robustness analysis for exploratory Figure 8 (answers review notes A-G)
import os, re, glob, json, numpy as np, pandas as pd
B='/anvil/projects/x-bio260593/vdgci/HQ_5054'; P=B+'/per_genome'; O=B+'/results/defmap/'
FNA=os.environ.get('SCRATCH','')+'/vdgci/HQ_5054/genomes/'
FL=10000
coh=pd.read_csv(B+'/results/cohort_VDGCI.csv'); phy=dict(zip(coh.accession,coh.phylum))
loc=pd.read_csv(B+'/results/extended/defense_systems_located.tsv',sep='\t')
loc['prim']=loc.in_viral_dna.astype(str)=='True'
KW=[('I',r'integrase|recombinase|excisionase|resolvase'),('S',r'capsid|portal|terminase|tail|baseplate|head|prohead|scaffold|fiber|fibre|tape.?measure|neck|sheath|tube|virion|structural'),
    ('L',r'holin|lysin|lysozyme|spanin|peptidoglycan|lysis'),('R',r'repressor|regulator|transcription|cro\b|\bci\b|polymerase|primase|helicase|replicat|antiterminat|anti-terminat|dna.binding|recombination|exonuclease|ssb|single.strand'),('M',r'transposase|insertion')]
def gcl(x):
    if not isinstance(x,str) or x in ('NA',''): return 'U'
    x=x.lower()
    for c,p in KW:
        if re.search(p,x): return c
    return 'O'
def clens(acc):
    f=FNA+acc+'.fna'; L={}; cur=None
    if not os.path.exists(f): return None
    with open(f) as fh:
        for l in fh:
            if l[0]=='>': cur=l[1:].split()[0]; L[cur]=0
            else: L[cur]+=len(l.strip())
    return L
def merge(iv):
    iv=sorted(iv); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]+1: out[-1][1]=max(out[-1][1],b)
        else: out.append([a,b])
    return out
def subtract(iv,cut):
    res=[]
    for a,b in iv:
        segs=[[a,b]]
        for c,d in cut:
            new=[]
            for x,y in segs:
                if d<x or c>y: new.append([x,y])
                else:
                    if c>x: new.append([x,c-1])
                    if d<y: new.append([d+1,y])
            segs=new
        res+=segs
    return res
led={'located_total':int(len(loc)),'primary_in_viral':int(loc.prim.sum())}
reasons={}; sysrows=[]; per=[]; prows=[]; flank_side_bp={'anchor':0,'far':0}; nolen=0
LG=dict(tuple(loc.groupby('accession')))
for acc in coh.accession:
    gl=LG.get(acc,loc.iloc[0:0]); i0=len(sysrows)
    d=f'{P}/{acc}'
    hf=glob.glob(d+'/*_protein_headers.txt'); gf=glob.glob(d+'/*_defense_finder_genes.tsv'); vf=glob.glob(d+'/*_virus_summary.tsv'); vg=glob.glob(d+'/*_virus_genes.tsv')
    H={}
    if hf:
        for l in open(hf[0]):
            if l[0]=='>':
                x=l[1:].split(' # '); H[x[0]]=(int(x[1]),int(x[2]))
    CL=clens(acc)
    if CL is None:
        nolen+=1; CL={}
        for k,(a,b) in H.items():
            c=k.rsplit('_',1)[0]; CL[c]=max(CL.get(c,0),b)
    prov=[]; whole=set()
    if vf:
        v=pd.read_csv(vf[0],sep='\t')
        for sn in v.seq_name:
            if '|provirus_' in sn:
                c,r=sn.split('|provirus_'); a,b=map(int,r.split('_')); prov.append((c,a,b,sn))
            else: whole.add(sn)
    G=None
    if vg and prov:
        G=pd.read_csv(vg[0],sep='\t',usecols=['gene','start','end','annotation_description'])
        G=G[G.gene.str.contains(r'\|provirus_')].copy(); G['prov']=G.gene.str.rsplit('_',n=1).str[0]; G['cls']=G.annotation_description.map(gcl)
    # anchors per provirus
    anchor={}
    for c,a,b,sn in prov:
        if G is None: continue
        g=G[(G.prov==sn)&(G.cls=='I')]
        if len(g):
            r=(((g.start+g.end)/2-a)/(b-a)).values; k=np.argmin(np.minimum(r,1-r)); anchor[sn]='left' if r[k]<0.5 else 'right'
    # flank windows: union minus viral
    viral=[(c,a,b) for c,a,b,_ in prov]+[(w,1,CL.get(w,0)) for w in whole]
    flank_bp=0; fwin={}
    for c in set(x[0] for x in prov):
        L=CL.get(c,0); wins=[]; sidewins=[]
        for cc,a,b,sn in prov:
            if cc!=c: continue
            lw=[max(1,a-FL),a-1]; rw=[b+1,min(L,b+FL)] if L else [b+1,b+FL]
            for w,side in [(lw,'left'),(rw,'right')]:
                if w[1]>=w[0]:
                    wins.append(w)
                    s = None if sn not in anchor else ('anchor' if side==anchor[sn] else 'far')
                    sidewins.append((w,sn,s))
        cut=[[a,b] for cc,a,b in viral if cc==c]
        u=subtract(merge(wins),merge(cut)); flank_bp+=sum(y-x+1 for x,y in u); fwin[c]=(u,sidewins)
        for w,sn,s in sidewins:
            if s: flank_side_bp[s]+=sum(y-x+1 for x,y in subtract([w],merge(cut)))
    viral_bp=sum(b-a+1 for c,a,b,_ in prov)+sum(CL.get(w,0) for w in whole)
    gtab=pd.read_csv(gf[0],sep='\t',usecols=['hit_id','sys_id']) if gf else pd.DataFrame(columns=['hit_id','sys_id'])
    defco=[H[h] for h in gtab.hit_id if h in H]; defctg=[h.rsplit('_',1)[0] for h in gtab.hit_id if h in H]
    hits=gtab.groupby('sys_id').hit_id.apply(list).to_dict()
    for _,s in gl.iterrows():
        sid=s.sys_id; rec=dict(accession=acc,sys_id=sid,type=s.type,mech=s.mechanism_class,prim=bool(s.prim))
        if not hf: rec['cat']='unresolved:no_protein_headers'
        elif sid not in hits: rec['cat']='unresolved:system_not_in_gene_table'
        else:
            co=[H[h] for h in hits[sid] if h in H]
            if not co: rec['cat']='unresolved:hits_not_in_headers'
            else:
                c=hits[sid][0].rsplit('_',1)[0]; ss=min(min(x) for x in co); se=max(max(x) for x in co); m=(ss+se)/2
                rec.update(contig=c,sys_start=ss,sys_end=se,n_hits=len(co))
                cat='elsewhere'
                if c in whole: cat='inside_viral_record'
                else:
                    for cc,a,b,sn in prov:
                        if cc==c and ss>=a and se<=b: cat='inside_provirus'; rec.update(provirus=sn,prov_start=a,prov_end=b); break
                        if cc==c and se>=a and ss<=b: cat='straddle'; rec.update(provirus=sn,prov_start=a,prov_end=b)
                    if cat=='elsewhere' and c in fwin:
                        u,sw=fwin[c]
                        if any(x<=m<=y for x,y in u):
                            cat='flank'
                            for w,sn,sd in sw:
                                if w[0]<=m<=w[1]: rec['flank_side']=sd or 'unoriented'; break
                rec['cat']=cat
        sysrows.append(rec)
    per.append(dict(accession=acc,genome_bp=int(coh.set_index('accession').genome_bp.get(acc,0)),viral_bp_recomputed=int(viral_bp),flank_bp=int(flank_bp)))
    # positional records for inside_provirus systems
    if G is not None:
        for rec in sysrows[i0:]:
            if rec.get('cat')!='inside_provirus': continue
            sn=rec['provirus']; a,b=rec['prov_start'],rec['prov_end']; L=b-a
            g=G[G.prov==sn].sort_values('start')
            if not len(g): continue
            isdef=np.array([any(cg==sn.split('|')[0] and x<=e and y>=st for (x,y),cg in zip(defco,defctg)) for st,e in zip(g.start,g.end)])
            mids=((g.start.values+g.end.values)/2-a)/L
            k=int(((g.start<=rec['sys_end'])&(g.end>=rec['sys_start'])).sum()) or 1
            n=len(g); st=g.start.values; en=g.end.values
            win=np.array([((st[i]+en[i+k-1])/2-a)/L for i in range(0,max(1,n-k+1))]) if n>=k else np.array([0.5])
            rel=((rec['sys_start']+rec['sys_end'])/2-a)/L
            side=anchor.get(sn)
            prows.append(dict(accession=acc,phylum=phy.get(acc,''),provirus=sn,sys_id=rec['sys_id'],type=rec['type'],mech=rec['mech'],k_genes=k,n_prov_genes=n,
                rel=rel,anchor_side=side,nullA=';'.join(f'{x:.4f}' for x in mids[~isdef]),nullB=';'.join(f'{x:.4f}' for x in win)))
led['n_genomes_without_fasta_lengths']=nolen
S=pd.DataFrame(sysrows); S.to_csv(O+'audit_system_ledger.tsv.gz',sep='\t',index=False)
D=pd.DataFrame(per); D.to_csv(O+'audit_per_genome_bp.tsv',sep='\t',index=False)
led['category_counts_all']=S.cat.value_counts().to_dict()
led['category_counts_primary622']=S[S.prim].cat.value_counts().to_dict()
led['category_counts_nonprimary']=S[~S.prim].cat.value_counts().to_dict()
led['bp']=dict(genome_bp=int(D.genome_bp.sum()),viral_bp_recomputed=int(D.viral_bp_recomputed.sum()),viral_bp_cohort_file=int(coh.viral_bp.sum()),flank_bp_union=int(D.flank_bp.sum()),flank_bp_anchor_side=int(flank_side_bp['anchor']),flank_bp_far_side=int(flank_side_bp['far']))
json.dump(led,open(O+'audit_ledger.json','w'),indent=1)
pd.DataFrame(prows).to_csv(O+'audit_positions.tsv.gz',sep='\t',index=False)
print(json.dumps(led,indent=1))
