import shutil
S='/anvil/projects/x-bio260593/VDGCI_cluster_pipeline/scripts/'
def patch(fn, after_sg, after_loop, ind):
    p=S+fn; shutil.copy(p,p+'.bak_dup_20261024')
    L=open(p).read().split('\n')
    assert 'sg[r.get("sys_id", "")].append' in L[after_sg-1] and 'DictReader(open(sf[0])' in L[after_loop-1], fn
    L.insert(after_loop, ' '*(ind+4)+'if r.get("sys_id", "") in seen: continue  # DefenseFinder can list a system twice')
    L.insert(after_loop+1, ' '*(ind+4)+'seen.add(r.get("sys_id", ""))')
    L.insert(after_sg, ' '*ind+'sg = {k: list(dict.fromkeys(v)) for k, v in sg.items()}; seen = set()  # de-duplicate (v3.2)')
    open(p,'w').write('\n'.join(L)); print('patched',fn)
patch('05_analyze.py',86,88,4)
patch('06_extended_analysis.py',118,119,8)
