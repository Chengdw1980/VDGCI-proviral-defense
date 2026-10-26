import pandas as pd, gzip, numpy as np
S='/tmp/claude-0/-home-claude/1d89e018-43ce-5a67-b16c-37382686b1fc/scratchpad'
P=S+'/pipe/VDGCI_cluster_pipeline/config/'
m=pd.read_csv(P+'cohort_manifest.tsv',sep='\t'); a=pd.read_csv(P+'assembly_stats_HQ.tsv',sep='\t')
full=pd.read_csv(P+'cohort_manifest_full.tsv',sep='\t')
c=m.merge(a[['accession','n_contigs','contig_n50','ncbi_annotated']],on='accession')
print(len(m),len(c), c.phylum.nunique(), c['class'].nunique(), c['order'].nunique(), c.family.nunique())
print(c.ncbi_assembly_level.value_counts())
print('GCF',c.accession.str.startswith('GCF').sum())
print('annotated',(c.ncbi_annotated=='yes').sum())
print('N50<20kb',(c.contig_n50<20000).sum(), 'single contig',(c.n_contigs==1).sum())
print('size Mb', (c.genome_size_bp/1e6).describe().round(2).to_dict())
print('compl', c.checkm2_completeness.describe().round(2).to_dict(), 'contam', c.checkm2_contamination.describe().round(2).to_dict())
pc=c.phylum.value_counts(); print(pc.head(15).to_dict()); 
for k in [1,5,10,20,30,50,100]: print('phyla >=',k, (pc>=k).sum(), 'genomes', pc[pc>=k].sum())
# 48 overlap via taxonomy
tax=pd.read_csv(S+'/gtdb/bac120_taxonomy_r232.tsv.gz',sep='\t',header=None,names=['gid','tax'])
tax['acc']=tax.gid.str[3:]; tax['core']=tax.acc.str[4:].str.split('.').str[0]
t48=pd.read_csv(P+'cohort_manifest_48.tsv',sep='\t'); print(t48.columns.tolist())
acc48=t48['accession'] if 'accession' in t48 else t48.iloc[:,0]
core48=acc48.str[4:].str.split('.').str[0]
hit=tax[tax.core.isin(core48)].drop_duplicates('core'); print('48 in GTDB taxonomy:',len(hit))
hit=hit.assign(family=hit.tax.str.extract(r'f__([^;]*)'), phylum=hit.tax.str.extract(r'p__([^;]*)'), species=hit.tax.str.extract(r's__(.*)$'))
fam48=set(hit.family); print('families of 48:',len(fam48), 'in HQ cohort families:', len(fam48 & set(c.family)))
ccore=set(c.accession.str[4:].str.split('.').str[0]); print('identical genomes shared:', len(set(core48)&ccore))
hit.to_csv('gtdb_tax_48.tsv',sep='\t',index=False)
c.to_csv('cohort5054_annot.tsv',sep='\t',index=False)
print('full 6048 phyla',full.phylum.nunique())
