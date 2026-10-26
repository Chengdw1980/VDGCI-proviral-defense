import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, gzip, treeswift
from physig import vcv
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,
 'axes.edgecolor':'#52514e','axes.labelcolor':'#0b0b0b','xtick.color':'#52514e','ytick.color':'#52514e'})
B,O,A,Y,M,G,V,R='#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4','#008300','#6250d6','#e34948'
S='/tmp/claude-0/-home-claude/1d89e018-43ce-5a67-b16c-37382686b1fc/scratchpad'
c=pd.read_csv('cohort5054_annot.tsv',sep='\t'); v=pd.read_csv('cohort48_with_gtdb_tips.csv'); h=pd.read_csv('gtdb_tax_48.tsv',sep='\t')
pw=pd.read_csv('power_group_size.csv')
fig,ax=plt.subplots(1,3,figsize=(13,4.6),gridspec_kw={'width_ratios':[1.25,1,1]})
import pandas as _pd
_D='/tmp/claude-0/-home-claude/1d89e018-43ce-5a67-b16c-37382686b1fc/scratchpad/in/DC Viral Genomic Fossils for Bacterial Evolution Work Project 09-11-2026/'
_sd=_pd.read_csv(_D+'defense_systems_detail.csv').drop_duplicates(['accession','sys_id'])
_m=_sd.groupby('accession').size()
v['n_defense_systems']=[(_m.get(a,0) if not _pd.isna(n) else n) for a,n in zip(v.accession,v.n_defense_systems)]
# ---- phylo figure
t=treeswift.read_tree_newick(gzip.open(S+'/gtdb/bac120_r232.tree.gz','rt').read())
t=t.extract_tree_with(set(v.tip)); t.suppress_unifurcations()
order=[l.label for l in t.traverse_leaves()]; ypos={lab:i for i,lab in enumerate(order)}
xpos={}; 
for n in t.traverse_preorder(): xpos[n]=0 if n.is_root() else xpos[n.parent]+(n.edge_length or 0)
def gety(n):
    if n.is_leaf(): return ypos[n.label]
    ys=[gety(ch) for ch in n.children]; return (min(ys)+max(ys))/2
yy={n:gety(n) for n in t.traverse_postorder()}
fig=plt.figure(figsize=(11,9)); axt=fig.add_axes([0.03,0.06,0.33,0.80]); axh=fig.add_axes([0.62,0.06,0.25,0.80])
for n in t.traverse_preorder():
    if n.is_root(): continue
    p=n.parent; axt.plot([xpos[p],xpos[n]],[yy[n],yy[n]],color='#52514e',lw=0.8)
    axt.plot([xpos[p],xpos[p]],[yy[p],yy[n]],color='#52514e',lw=0.8)
phy=v.set_index('tip').phylum; ph_order=v.phylum.value_counts().index.tolist()
pal=[B,O,A,Y,M,G,V,R]+['#9a9890']*20; pcol={p:pal[i] if i<8 else '#9a9890' for i,p in enumerate(ph_order)}
mx=max(xpos.values())
for lab in order:
    axt.scatter(xpos[[l for l in t.traverse_leaves() if l.label==lab][0]],ypos[lab],s=18,color=pcol[phy[lab]],zorder=3)
axt.axis('off'); axt.set_ylim(-1,len(order))
sp=v.set_index('tip').species
for lab in order: fig.text(0.365,0.06+0.80*(ypos[lab]+1)/(len(order)+1),sp[lab][:38],fontsize=6.3,va='center',style='italic')
cols=[('pct_viral','% viral'),('n_prov','Proviruses'),('n_defense_systems','Defense'),('bp','Genome Mb')]
M_=np.array([[v.set_index('tip').loc[lab,c] for c,_ in cols] for lab in order],float); M_[:,3]/=1e6
Mn=np.nan_to_num((M_-np.nanmin(M_,0))/(np.nanmax(M_,0)-np.nanmin(M_,0)))
from matplotlib.colors import LinearSegmentedColormap
cm=LinearSegmentedColormap.from_list('b',['#f2f6fc','#9ec3ef','#2a78d6','#0d3f80'])
axh.imshow(Mn,aspect='auto',cmap=cm,origin='lower',extent=[-0.5,3.5,-0.5,len(order)-0.5])
for i in range(len(order)):
    for j in range(4):
        val=M_[i,j]; s='' if np.isnan(val) else (f'{val:.1f}' if j in (0,3) else f'{int(val)}')
        axh.text(j,i,s,ha='center',va='center',fontsize=5.8,color='white' if Mn[i,j]>0.55 else '#0b0b0b')
axh.set_xticks(range(4)); axh.set_xticklabels([l for _,l in cols],fontsize=8); axh.xaxis.tick_top(); axh.set_yticks([])
axh.set_ylim(-1,len(order))
for s in axh.spines.values(): s.set_visible(False)
# legend phyla
from matplotlib.lines import Line2D
hs=[Line2D([],[],marker='o',ls='',color=pcol[p],label=p) for p in ph_order[:8]]+[Line2D([],[],marker='o',ls='',color='#9a9890',label='Other phyla')]
fig.legend(handles=hs,loc='lower left',bbox_to_anchor=(0.88,0.06),frameon=False,fontsize=7,title='Phylum',title_fontsize=8)
ps=pd.read_csv('../dedup/physig48_gtdb_dedup.csv')
txt='Phylogenetic signal on the GTDB bac120 tree (log-transformed counts; 999 permutations for K; LRT for λ)\n'+'\n'.join(f"{r.trait.replace(' (count)','').replace(' (positive control)',' (control)')[:24]:24s} K={r.Blomberg_K:.2f} (P={r.P_K:.3f})  λ={r.Pagel_lambda:.2f} (P={r.P_lambda:.1g})" for r in ps.itertuples())
fig.text(0.40,0.905,txt,fontsize=6.6,family='DejaVu Sans Mono',va='bottom')
plt.savefig('../dedup/Figure_6_v10_Phylogeny_GTDB.png',dpi=300,facecolor='white')
print('ok')
