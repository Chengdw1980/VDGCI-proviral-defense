import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import poisson
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.linewidth':0.6,
                     'axes.spines.top':False,'axes.spines.right':False,
                     'xtick.major.width':0.6,'ytick.major.width':0.6,
                     'xtick.major.size':2.5,'ytick.major.size':2.5})
INK='#222222'; MUTED='#6b6b6b'; SIG='#D55E00'; NS='#9a9a9a'; POOL='#3B7DBF'
POOLED=2.603; POOL_CI=(2.323,2.867)
C=pd.read_csv('data/class_enrichment.tsv',sep='\t')
O=pd.read_csv('data/order_enrichment.tsv',sep='\t')

fig=plt.figure(figsize=(7.2,8.6))
gs=fig.add_gridspec(2,2,height_ratios=[1.95,1.0],width_ratios=[1.32,1.0],
                    hspace=0.30,wspace=0.42,left=0.175,right=0.975,top=0.955,bottom=0.065)
def lab(ax,t,x,y=1.012):
    ax.text(x,y,t,transform=ax.transAxes,fontsize=9,fontweight='bold',va='bottom')

def forest(ax,D,title,xmax,xlab='Fold-enrichment of defense systems in proviral DNA (95% CI)'):
    D=D.sort_values('fold').reset_index(drop=True)
    ax.axvspan(POOL_CI[0],POOL_CI[1],color=POOL,alpha=0.12,lw=0,zorder=0)
    ax.axvline(POOLED,color=POOL,lw=0.9,zorder=1)
    ax.axvline(1,ls='--',lw=0.6,color=MUTED,zorder=1)
    for i,r in D.iterrows():
        sig=r.ci_lo>1
        col=SIG if sig else NS
        hi=min(r.ci_hi,xmax)
        ax.plot([r.ci_lo,hi],[i,i],color=col,lw=1.25,solid_capstyle='butt',zorder=3)
        if r.ci_hi>xmax:
            ax.annotate('',xy=(xmax+0.22,i),xytext=(xmax-0.05,i),zorder=3,
                        arrowprops=dict(arrowstyle='-|>',color=col,lw=0.9,mutation_scale=5))
        ms=3.0+2.4*np.sqrt(r.n_genomes/416)
        ax.plot(r.fold,i,'o',ms=ms,color=col,mec='white',mew=0.7,zorder=4)
    ax.set_yticks(range(len(D)))
    ax.set_yticklabels([f"{t}" for t in D.taxon],fontsize=6.0)
    for t,sig in zip(ax.get_yticklabels(),D.ci_lo>1):
        t.set_color(INK if sig else MUTED)
    ax.set_ylim(-0.8,len(D)-0.2); ax.set_xlim(0,xmax+0.45)
    ax.set_xlabel(xlab,fontsize=6.8)
    ax.set_title(title,fontsize=7.5,loc='left',color=INK,pad=4)
    return D

ax=fig.add_subplot(gs[:,0]); lab(ax,'a',-0.40)
forest(ax,C,'Classes with ≥30 genomes (n = 30)',10.0)
ax.text(0.985,0.012,'orange, 95% CI excludes 1\nblue line, pooled estimate (2.60)\ncircle area ∝ genomes',
        transform=ax.transAxes,ha='right',va='bottom',fontsize=5.9,color=MUTED,linespacing=1.45)

ax=fig.add_subplot(gs[0,1]); lab(ax,'b',-0.30)
forest(ax,O,'Orders with ≥30 genomes (n = 19)',8.0,'Fold-enrichment (95% CI)')

# c  funnel: homogeneity check
ax=fig.add_subplot(gs[1,1]); lab(ax,'c',-0.30)
E=np.logspace(np.log10(0.18),np.log10(30),300)
for q,ls in [((0.025,0.975),'-')]:
    lo=poisson.ppf(q[0],POOLED*E)/E; hi=poisson.ppf(q[1],POOLED*E)/E
    ax.plot(E,lo,color=MUTED,lw=0.6,ls='--',alpha=0.8); ax.plot(E,hi,color=MUTED,lw=0.6,ls='--',alpha=0.8)
ax.axhline(POOLED,color=POOL,lw=0.9)
inb=0
for _,r in C.iterrows():
    lo=poisson.ppf(0.025,POOLED*r.expected)/r.expected
    hi=poisson.ppf(0.975,POOLED*r.expected)/r.expected
    ok=lo<=r.fold<=hi; inb+=ok
    ax.plot(r.expected,min(r.fold,13.5),'o',ms=3.4,mfc=(NS if ok else SIG),mec='white',mew=0.6)
ax.set_xscale('log'); ax.set_xlim(0.18,34); ax.set_ylim(-0.4,14)
ax.set_xlabel('Defense systems expected in proviral DNA',fontsize=6.8)
ax.set_ylabel('Fold-enrichment',fontsize=6.8)
ax.set_title('Homogeneity of class estimates',fontsize=7.5,loc='left',color=INK,pad=4)
ax.text(0.97,0.95,f'{inb}/30 classes within the 95%\ninterval expected under a\nsingle common enrichment\n\nPermutation test of\nheterogeneity: P = 0.58',
        transform=ax.transAxes,ha='right',va='top',fontsize=6.0,color=INK,linespacing=1.5)
fig.savefig('Figure_7_v11_Lineage.png',dpi=400)
fig.savefig('Figure_7_v11_Lineage.pdf')
print('in bounds',inb,'of',len(C))
print('class CI>1:',int((C.ci_lo>1).sum()),'order CI>1:',int((O.ci_lo>1).sum()))
print('class CI<1:',int((C.ci_hi<1).sum()),'order CI<1:',int((O.ci_hi<1).sum()))
