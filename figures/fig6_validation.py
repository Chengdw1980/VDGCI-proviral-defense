import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#52514e','xtick.color':'#52514e','ytick.color':'#52514e'})
B,O,A,GR,V='#2a78d6','#eb6834','#1baf7a','#9a9890','#6250d6'
ph=pd.read_csv('enrichment_by_phylum.csv'); me=pd.read_csv('enrichment_by_mechanism.csv'); dc=pd.read_csv('decay_vs_defense.csv'); mg=pd.read_csv('mge_comparison.csv'); mg.loc[0,['boot95_lo','boot95_hi']]=[2.336,2.921]
fig=plt.figure(figsize=(12,8.2)); g=fig.add_gridspec(2,3,width_ratios=[1.25,1,0.8],height_ratios=[1,1],wspace=0.55,hspace=0.45)
a=fig.add_subplot(g[:,0]); ph=ph.iloc[::-1].reset_index(drop=True); XMAX=12
for i,r in ph.iterrows():
    if r.n_in_viral_dna==0:
        a.plot(0,i,'o',mfc='white',mec=GR,ms=5); continue
    sig=r.boot95_lo>1; c=B if sig else GR
    a.plot([r.boot95_lo,min(r.boot95_hi,XMAX)],[i,i],color=c,lw=1.8)
    if r.boot95_hi>XMAX: a.annotate('',xy=(XMAX+0.3,i),xytext=(XMAX-0.2,i),arrowprops=dict(arrowstyle='->',color=c,lw=1.5))
    a.plot(r.fold,i,'o',color=c,ms=6,mec='white',mew=1)
a.axvline(1,color='#52514e',ls='--',lw=0.9); a.axvline(2.624,color=O,lw=1.2,alpha=0.8)
a.text(2.7,len(ph)-0.3,'cohort 2.6×',color=O,fontsize=7.5)
a.set_yticks(range(len(ph))); a.set_yticklabels([f'{p} ({n})' for p,n in zip(ph.phylum,ph.n_genomes)],fontsize=7.5)
a.set_xlim(-0.3,XMAX+0.6); a.set_xlabel('Fold-enrichment of defense systems in viral DNA\n(genome bootstrap 95% CI)')
a.set_title('a  Replication within phyla (≥30 genomes)',loc='left',fontweight='bold')
a.text(0.98,0.01,'blue: CI excludes 1 (9 of 25)\ngrey: CI includes 1\nopen circle: no system in viral DNA',transform=a.transAxes,ha='right',va='bottom',fontsize=7,color='#52514e')
b=fig.add_subplot(g[0,1:]); classes=['Abortive infection / toxin-antitoxin','Other or mechanism unresolved','Direct nucleic-acid targeting','ALL']
labs=['Abortive infection /\ntoxin–antitoxin','Other / unresolved','Direct nucleic-acid\ntargeting','All systems']
sets=[('all','All 5,054',B),('contigN50>=20000','Scaffold N50 ≥ 20 kb (4,327)',A),('finished_assemblies','Finished (1,205)',V)]
disc={'Abortive infection / toxin-antitoxin':(5.51,2.83,8.58),'Other or mechanism unresolved':(3.58,1.63,6.17),'Direct nucleic-acid targeting':(2.83,1.17,5.09),'ALL':(3.68,2.35,5.22)}
for j,cl in enumerate(classes):
    for k,(s,lab,c) in enumerate(sets):
        r=me[(me.genome_set==s)&(me.mechanism_class==cl)].iloc[0]; x=j+(k-1.5)*0.18
        b.plot([x,x],[r.boot95_lo,r.boot95_hi],color=c,lw=2); b.plot(x,r.fold,'o',color=c,ms=6,mec='white',mew=1,label=lab if j==0 else None)
    d=disc[cl]; x=j+1.5*0.18; b.plot([x,x],[d[1],d[2]],color=GR,lw=2); b.plot(x,d[0],'D',color=GR,ms=5,mec='white',label='Discovery (48)' if j==0 else None)
b.axhline(1,color='#52514e',ls='--',lw=0.9); b.set_xticks(range(4)); b.set_xticklabels(labs,fontsize=8); b.set_ylabel('Fold-enrichment (95% CI)')
b.legend(frameon=False,fontsize=7.5,ncol=4,loc='upper center',bbox_to_anchor=(0.5,1.02)); b.set_ylim(0,9)
b.set_title('b  By mechanism class and genome set',loc='left',fontweight='bold',pad=18)
c=fig.add_subplot(g[1,1]); xs=[0,1]; cols=[O,B]
for x,(i,r),col in zip(xs,mg.iterrows(),cols):
    c.bar(x,r.fold,color=col,width=0.55); c.plot([x,x],[r.boot95_lo,r.boot95_hi],color='#0b0b0b',lw=1.4)
    c.text(x,r.boot95_hi+0.2,f'{r.fold:.1f}×',ha='center',fontweight='bold')
c.set_xticks(xs); c.set_xticklabels([f'Proviruses\n({mg.pct_sequence[0]:.2f}% of seq.)',f'geNomad plasmids\n({mg.pct_sequence[1]:.2f}% of seq.)'],fontsize=8)
c.axhline(1,color='#52514e',ls='--',lw=0.9); c.set_ylabel('Fold-enrichment (95% CI)'); c.set_ylim(0,6.5)
c.set_title('c  Viral DNA vs plasmids',loc='left',fontweight='bold')
d=fig.add_subplot(g[1,2]); dd=dc.iloc[:3]
d.bar(range(3),dd.systems_per_100kb,color=[A,'#7fcfae','#c5e9da'],width=0.6)
for i,r in dd.iterrows(): d.text(i,r.systems_per_100kb+0.02,f'{r.systems_per_100kb:.2f}\n(n={r.n_defense_systems})',ha='center',fontsize=7.5)
d.set_xticks(range(3)); d.set_xticklabels(['High','Medium','Low'],fontsize=8); d.set_xlabel('CheckV provirus quality')
d.set_ylabel('Defense systems per 100 kb'); d.set_ylim(0,0.85)
d.set_title('d  Across proviral completeness tiers',loc='left',fontweight='bold')
plt.savefig('Figure_7_v6_Validation.png',dpi=300,bbox_inches='tight',facecolor='white'); print('ok')
