import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#52514e','xtick.color':'#52514e','ytick.color':'#52514e'})
B,O,A,GR='#2a78d6','#eb6834','#1baf7a','#b5b3ab'
D='/tmp/claude-0/-home-claude/1d89e018-43ce-5a67-b16c-37382686b1fc/scratchpad/in/DC Viral Genomic Fossils for Bacterial Evolution Work Project 09-11-2026/'
s=pd.read_csv(D+'defense_systems_detail.csv').drop_duplicates(['accession','sys_id']); mc=pd.read_csv('../dedup/mechanism_class_48.csv')
fig=plt.figure(figsize=(12,4.8)); g=fig.add_gridspec(1,3,width_ratios=[0.7,1.1,1.3])
a=fig.add_subplot(g[0]); a.bar([0,1],[1.91,7.19],color=[GR,O],width=0.6)
for x,v in [(0,1.91),(1,7.19)]: a.text(x,v+0.2,f'{v:.1f}%',ha='center',fontweight='bold')
a.set_xticks([0,1]); a.set_xticklabels(['Viral-derived\nsequence','Defense systems\nin viral DNA'],fontsize=8)
a.set_ylabel('% of cohort total'); a.set_ylim(0,11)
a.text(0.5,8.3,'3.8-fold\nbinomial P ≈ 6×10⁻¹⁰\ngenome bootstrap\n95% CI 2.3–5.3',ha='center',color='#0b0b0b',fontsize=8)
a.set_title('a  Overall enrichment (n = 48)',loc='left',fontweight='bold')
b=fig.add_subplot(g[1]); order=['Abortive infection / toxin-antitoxin','Other or mechanism unresolved','Direct nucleic-acid targeting']
m=mc.set_index('mechanism_class').loc[order]; y=np.arange(3)[::-1]; cols=[O,GR,B]
for yi,(nm,r),c in zip(y,m.iterrows(),cols):
    b.plot([r.boot_lo,r.boot_hi],[yi,yi],color=c,lw=2.2); b.plot(r.fold,yi,'o',color=c,ms=9,markeredgecolor='white',markeredgewidth=1.5)
    b.text(r.boot_hi+0.25,yi,f'{r.fold:.1f}× ({int(r.n_viral_derived)}/{int(r.n_systems)})',va='center',fontsize=8)
b.axvline(1,color='#52514e',lw=1,ls='--'); b.set_yticks(y); b.set_yticklabels(['Abortive infection /\ntoxin–antitoxin','Other / mechanism\nunresolved','Direct nucleic-acid\ntargeting (RM, CRISPR-Cas, …)'],fontsize=8)
b.set_xlim(0,11); b.set_xlabel('Fold-enrichment in viral-derived DNA\n(point; genome-level bootstrap 95% CI)')
b.set_title('b  By predicted mechanism class',loc='left',fontweight='bold')
c=fig.add_subplot(g[2]); t=s.groupby('type').in_prophage.agg(['count','sum']).sort_values('count',ascending=False).head(14)[::-1]
yy=np.arange(len(t)); c.barh(yy,t['count']-t['sum'],color=GR,label='Elsewhere in genome',edgecolor='white',linewidth=1)
c.barh(yy,t['sum'],left=t['count']-t['sum'],color=O,label='In viral-derived DNA',edgecolor='white',linewidth=1)
c.set_yticks(yy); c.set_yticklabels(t.index,fontsize=7.5); c.set_xlabel('Defense systems (cohort total)')
c.legend(frameon=False,fontsize=8,loc='lower right'); c.set_title('c  Most frequent defense types',loc='left',fontweight='bold')
plt.tight_layout(); plt.savefig('../dedup/Figure_4_v10_Defense.png',dpi=300,facecolor='white'); print('ok')
