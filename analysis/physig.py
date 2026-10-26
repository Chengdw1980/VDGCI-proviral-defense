"""Phylogenetic signal: Blomberg's K (permutation test) and Pagel's lambda (ML, LRT vs lambda=0).
Reusable for the 5,054 cohort. Implemented with numpy (no R)."""
import numpy as np, treeswift
from scipy import optimize, stats

def vcv(tree_nwk_or_obj, tips):
    t = treeswift.read_tree_newick(tree_nwk_or_obj) if isinstance(tree_nwk_or_obj,str) else tree_nwk_or_obj
    keep=set(tips); t=t.extract_tree_with(keep); t.suppress_unifurcations()
    root=t.root; dist={}
    for n in t.traverse_preorder():
        dist[n]= 0.0 if n is root else dist[n.parent]+(n.edge_length or 0.0)
    leaves={l.label:l for l in t.traverse_leaves()}
    # shared path length = depth of MRCA
    anc={}
    for lab,l in leaves.items():
        path=[]; n=l
        while n is not None: path.append(n); n=n.parent
        anc[lab]=path
    n=len(tips); C=np.zeros((n,n))
    for i,a in enumerate(tips):
        sa=set(anc[a])
        for j in range(i,n):
            b=tips[j]
            if i==j: C[i,i]=dist[leaves[a]]; continue
            for x in anc[b]:
                if x in sa: C[i,j]=C[j,i]=dist[x]; break
    return C

def _gls(x,C):
    Ci=np.linalg.inv(C); one=np.ones(len(x))
    a=(one@Ci@x)/(one@Ci@one); r=x-a
    return a,r,Ci

def blomberg_k(x,C,nperm=999,rng=np.random.default_rng(1)):
    x=np.asarray(x,float); n=len(x)
    def K(xx):
        a,r,Ci=_gls(xx,C)
        mse0=(r@r)/(n-1); mse=(r@Ci@r)/(n-1)
        exp=(np.trace(C)-n/np.sum(Ci))/(n-1)
        return (mse0/mse)/exp
    k=K(x); null=np.array([K(rng.permutation(x)) for _ in range(nperm)])
    return k,(1+np.sum(null>=k))/(nperm+1)

def pagel_lambda(x,C):
    x=np.asarray(x,float); n=len(x); D=np.diag(np.diag(C))
    def nll(lam):
        V=lam*C+(1-lam)*D
        sign,logdet=np.linalg.slogdet(V); Vi=np.linalg.inv(V); one=np.ones(n)
        a=(one@Vi@x)/(one@Vi@one); r=x-a; s2=(r@Vi@r)/n
        return 0.5*(n*np.log(2*np.pi*s2)+logdet+n)
    res=optimize.minimize_scalar(nll,bounds=(0,1),method='bounded')
    lam=res.x; l0=nll(0.0); lr=2*(l0-res.fun)
    p=0.5*stats.chi2.sf(max(lr,0),1)  # boundary-corrected
    return lam,lr,p
