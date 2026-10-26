import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False})
G='#009E73'; Gr='#9a9a9a'; O='#D55E00'
fig,(a,b)=plt.subplots(2,1,figsize=(7.5,10.2))
a.bar([0,1],[72,86],color=G,width=0.55,label='Corroborated by the other tool'); a.bar([0,1],[29,81],bottom=[72,86],color=Gr,width=0.55,label='Tool-specific')
for x,t,c,p in [(0,101,72,'71%'),(1,167,86,'51%')]:
    a.text(x,t+2,str(t),ha='center',fontweight='bold',fontsize=12); a.text(x,c/2,f'{c}\n({p})',ha='center',va='center',color='white',fontweight='bold',fontsize=11)
a.set_xticks([0,1]); a.set_xticklabels(['geNomad\nproviruses','PhiSpy\nprophages'],fontsize=11); a.set_ylabel('Number of prophage calls',fontsize=11)
a.set_ylim(0,195); a.legend(loc='upper left',fontsize=10); a.grid(axis='y',alpha=0.3); a.set_title('Two independent methods agree on most prophages',fontweight='bold',fontsize=12)
b.bar([0,1],[3.76,3.50],color=[O,G],width=0.5)
for x,v,t in [(0,3.76,'3.8×\np=6e-10'),(1,3.50,'3.5×\np=2e-06')]: b.text(x,v+0.08,t,ha='center',fontweight='bold',fontsize=12)
b.axhline(1,ls='--',color=Gr,lw=1); b.set_ylim(0,4.8); b.grid(axis='y',alpha=0.3)
b.set_xticks([0,1]); b.set_xticklabels(['All geNomad\nproviruses','Consensus\n(geNomad∩PhiSpy)'],fontsize=11); b.set_ylabel('Fold enrichment of defense systems',fontsize=11)
b.set_title('Defense-arsenal enrichment is robust to the prophage-calling method',fontweight='bold',fontsize=12)
for ax,l in [(a,'a'),(b,'b')]: ax.text(-0.12,1.08,l,transform=ax.transAxes,fontsize=20,fontweight='bold')
plt.tight_layout(h_pad=3); plt.savefig('Figure_5_v10_Consensus.png',dpi=200,facecolor='white'); print('ok')
