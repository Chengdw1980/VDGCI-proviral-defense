import os, re, glob, json, sys
import numpy as np, pandas as pd
B='/anvil/projects/x-bio260593/vdgci/HQ_5054'
P=B+'/per_genome'; OUT=B+'/results/defmap'; os.makedirs(OUT,exist_ok=True)
loc=pd.read_csv(B+'/results/extended/defense_systems_located.tsv',sep='\t')
loc=loc[loc.in_viral_dna.astype(str)=='True']
coh=pd.read_csv(B+'/results/cohort_VDGCI.csv')
phy=dict(zip(coh.accession,coh.phylum)) if 'phylum' in coh else {}
KW=[('integration',r'integrase|recombinase|excisionase|resolvase'),
    ('structural',r'capsid|portal|terminase|tail|baseplate|head|prohead|scaffold|fiber|fibre|tape.?measure|neck|sheath|tube|virion|structural'),
    ('lysis',r'holin|lysin|lysozyme|spanin|peptidoglycan|lysis'),
    ('replication_regulation',r'repressor|regulator|transcription|cro\b|\bci\b|polymerase|primase|helicase|replicat|antiterminat|anti-terminat|dna.binding|recombination|exonuclease|ssb|single.strand'),
    ('mobility',r'transposase|insertion')]
def gclass(desc):
    if not isinstance(desc,str) or desc in ('NA',''): return 'unannotated'
    d=desc.lower()
    for c,pat in KW:
        if re.search(pat,d): return c
    return 'other_annotated'
def hdr(acc):
    m={}
    f=glob.glob(f'{P}/{acc}/*_protein_headers.txt')
    if not f: return m
    for l in open(f[0]):
        if l.startswith('>'):
            x=l[1:].split(' # '); m[x[0]]=(int(x[1]),int(x[2]),int(x[3]))
    return m
rows=[]; genes_out=[]; nullpos=[]
for acc,g in loc.groupby('accession'):
    d=f'{P}/{acc}'
    vg=glob.glob(d+'/*_virus_genes.tsv'); dg=glob.glob(d+'/*_defense_finder_genes.tsv')
    if not vg or not dg: continue
    vg=pd.read_csv(vg[0],sep='\t'); vg=vg[vg.gene.str.contains(r'\|provirus_')]
    vg['prov']=vg.gene.str.rsplit('_',n=1).str[0]
    vg['cls']=vg.annotation_description.map(gclass)
    dfg=pd.read_csv(dg[0],sep='\t'); H=hdr(acc)
    cq={}
    cf=d+'/checkv_quality_summary.tsv'
    if os.path.exists(cf):
        c=pd.read_csv(cf,sep='\t'); cq=dict(zip(c.contig_id,c.checkv_quality))
    for _,s in g.iterrows():
        hits=dfg[dfg.sys_id==s.sys_id].hit_id.tolist()
        co=[H[h] for h in hits if h in H]
        if not co: continue
        contig=hits[0].rsplit('_',1)[0]
        ss=min(min(a,b) for a,b,_ in co); se=max(max(a,b) for a,b,_ in co)
        prov=None
        for pv in vg.prov.unique():
            ctg,rest=pv.split('|provirus_'); a,b=map(int,rest.split('_'))
            if ctg==contig and se>=a and ss<=b: prov=(pv,a,b); break
        if prov is None: continue
        pv,a,b=prov; L=b-a
        G=vg[vg.prov==pv].copy(); G['mid']=(G.start+G.end)/2; G['rel']=(G.mid-a)/L
        # flag defense genes (any system) in this provirus
        dco=[]
        for h in dfg.hit_id:
            if h in H and h.rsplit('_',1)[0]==contig:
                x,y,_=H[h]
                if y>=a and x<=b: dco.append((x,y))
        G['is_def']=[any(gs<=y and ge>=x for x,y in dco) for gs,ge in zip(G.start,G.end)]
        integ=G[G.cls=='integration']
        side=None
        if len(integ):
            r=integ.rel.values; k=np.argmin(np.minimum(r,1-r)); side='left' if r[k]<0.5 else 'right'
        rel=((ss+se)/2-a)/L
        orient=None if side is None else (rel if side=='left' else 1-rel)
        nd=G[~G.is_def]
        rows.append(dict(accession=acc,phylum=phy.get(acc,''),sys_id=s.sys_id,type=s.type,subtype=s.subtype,
            mech=s.mechanism_class,provirus=pv,prov_start=a,prov_end=b,prov_len=L,
            checkv=cq.get(pv,'NA'),sys_start=ss,sys_end=se,n_sys_genes=len(co),rel=round(rel,4),
            edge=round(min(rel,1-rel),4),integrase_side=side,orient=None if orient is None else round(orient,4),
            n_prov_genes=len(G),n_struct=int((G.cls=='structural').sum()),n_lysis=int((G.cls=='lysis').sum()),
            n_integ=len(integ),n_def_genes_in_prov=int(G.is_def.sum()),
            null_rel=';'.join(f'{v:.3f}' for v in nd.rel.values)))
        G['accession']=acc
        genes_out.append(G[['accession','prov','gene','start','end','strand','cls','is_def','annotation_description']])
R=pd.DataFrame(rows); R.to_csv(OUT+'/sys_positions.tsv',sep='\t',index=False)
GA=pd.concat(genes_out).drop_duplicates(subset=['gene'])
GA.to_csv(OUT+'/defense_provirus_genes.tsv.gz',sep='\t',index=False)
print('systems mapped',len(R),'of',len(loc),'proviruses',R.provirus.nunique(),'gene rows',len(GA))
