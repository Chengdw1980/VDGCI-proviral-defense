import pandas as pd, numpy as np
from scipy import stats
S='/tmp/claude-0/-home-claude/1d89e018-43ce-5a67-b16c-37382686b1fc/scratchpad'
D=S+'/in/DC Viral Genomic Fossils for Bacterial Evolution Work Project 09-11-2026/'
d=pd.read_csv(D+'defense_per_genome.csv'); sysd=pd.read_csv(D+'defense_systems_detail.csv')
d=d[d.n_defense_systems.notna()]
Sx=d.n_defense_systems.values; K=d.n_defense_systems_in_prophage.values; bp=d.genome_bp.values; vb=d.prophage_bp.values
def E(ix): return (K[ix].sum()/Sx[ix].sum())/(vb[ix].sum()/bp[ix].sum())
n0=len(d); print('genomes with defense data',n0,'E=%.2f'%E(np.arange(n0)))
rng=np.random.default_rng(11)
# genome-clustered bootstrap CI on the 48
bs=np.array([E(rng.integers(0,n0,n0)) for _ in range(5000)]); bs=bs[np.isfinite(bs)]
print('48-genome cluster-bootstrap 95%% CI for fold: %.2f-%.2f'%tuple(np.percentile(bs,[2.5,97.5])))
# overdispersion: variance of K share vs binomial
p=K.sum()/Sx.sum(); pis=K/np.maximum(Sx,1); 
phi=np.sum(Sx*(pis-p)**2)/(p*(1-p))/(n0-1); print('overdispersion phi=%.2f'%phi)
# power for cohort of n genomes, effect diluted by factor f (draft genomes under-call / weaker effect)
def power(n,f=1.0,reps=400,inner=300,alpha=0.05):
    hits=0
    for _ in range(reps):
        ix=rng.integers(0,n0,n)
        Ki=rng.binomial(Sx[ix], np.clip((K[ix]/np.maximum(Sx[ix],1))*f,0,1)) if f!=1 else K[ix]
        tot=lambda j:(Ki[j].sum()/Sx[ix][j].sum())/(vb[ix][j].sum()/bp[ix][j].sum())
        b=np.array([tot(rng.integers(0,n,n)) for _ in range(inner)]); b=b[np.isfinite(b)]
        hits+= np.percentile(b,100*alpha/2)>1
    return hits/reps
rows=[]
for n in [10,20,30,50,100]:
    for f in [1.0,0.6]:
        rows.append([n,f,power(n,f)])
pw=pd.DataFrame(rows,columns=['genomes_in_group','effect_scale','power_CI_excludes_1']); print(pw.to_string())
pw.to_csv('power_group_size.csv',index=False)
# per-defense-type expected counts at 5,054 (naive scaling) and binomial power w/ overdispersion + Bonferroni
pv=vb.sum()/bp.sum(); bt=sysd.groupby('type').in_prophage.agg(['count','sum']).reset_index()
scale=5054/n0; ntypes=bt.shape[0]; alpha=0.05/ntypes
def binpow(m,fold):
    m_eff=max(int(m/phi),1); q=min(fold*pv,0.99)
    crit=stats.binom.isf(alpha,m_eff,pv)  # reject if X>crit
    return stats.binom.sf(crit,m_eff,q)
bt['exp_5054']=(bt['count']*scale).round(0)
for fold in [2,3,5]: bt[f'power_fold{fold}']=[round(binpow(m,fold),2) for m in bt.exp_5054]
bt=bt.sort_values('count',ascending=False); print('types',ntypes,'viral seq frac %.4f'%pv); print(bt.head(25).to_string())
bt.to_csv('power_by_type.csv',index=False)
# minimum systems needed for 80% power at fold 3
for fold in [2,3,5]:
    m=1
    while binpow(m,fold)<0.8: m+=5
    print(f'fold {fold}: ~{m} systems of a type needed (after overdispersion phi) for 80% power at Bonferroni alpha')
