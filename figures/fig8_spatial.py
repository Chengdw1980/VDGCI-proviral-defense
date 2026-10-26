import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle, FancyBboxPatch, Polygon
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.linewidth':0.6,'axes.spines.top':False,'axes.spines.right':False,
                     'xtick.major.width':0.6,'ytick.major.width':0.6,'xtick.major.size':2.5,'ytick.major.size':2.5})
S=json.load(open('data/model_stats_v10.json'))
INK='#222222'; MUTED='#6b6b6b'; GRID='#e6e6e6'
C={'D':'#D55E00','S':'#3B7DBF','L':'#2E9E5B','I':'#6A3D9A','R':'#C9A227','M':'#9a9a9a','O':'#cfcfcf','U':'#ffffff'}
fig=plt.figure(figsize=(7.2,9.2))
gs=fig.add_gridspec(4,2,height_ratios=[1,1.05,1.3,0.8],hspace=0.72,wspace=0.30,left=0.115,right=0.98,top=0.97,bottom=0.02)
def lab(ax,t,x=-0.12,y=1.06): ax.text(x,y,t,transform=ax.transAxes,fontsize=9,fontweight='bold',va='bottom')
# a
ax=fig.add_subplot(gs[0,0]); lab(ax,'a')
rows=[('Fully\ninside',S['strict'],C['D']),('Inside or\ncrossing edge',S['loose'],'#e8a27a'),('Host DNA\n≤10 kb away',S['flank'],MUTED)]
for i,(name,d,col) in enumerate(rows):
    y=len(rows)-1-i
    ax.plot(d['ci'],[y,y],color=col,lw=1.6,solid_capstyle='round'); ax.plot(d['fold'],y,'o',ms=6,color=col,mec='white',mew=1)
    ax.text(d['ci'][1]+0.08,y,f"{d['fold']:.2f} ({d['ci'][0]:.2f}–{d['ci'][1]:.2f})\nn = {d['n']}",va='center',fontsize=6.3,color=INK)
ax.axvline(1,ls='--',lw=0.6,color=MUTED); ax.set_yticks(range(3)); ax.set_yticklabels([r[0] for r in rows][::-1])
ax.set_xlim(0.8,4.6); ax.set_ylim(-0.6,2.6); ax.set_xlabel('Defense-system fold-enrichment (95% CI)')
ax.set_title('Proviruses and their integration neighbourhood',fontsize=7.5,loc='left',color=INK)
# b
ax=fig.add_subplot(gs[0,1]); lab(ax,'b')
ob=np.array(S['orient_bins']['obs']); nb=np.array(S['orient_bins']['null']); x=np.arange(10)*0.1+0.05
ax.bar(x,ob,width=0.085,color=C['D'],label='Observed')
ax.step(np.r_[0,np.arange(10)*0.1+0.1],np.r_[nb[0],nb],where='pre',color=INK,lw=1,label='Expected')
ax.set_xlim(0,1); ax.set_ylim(0,44); ax.set_xticks([0,0.5,1]); tl=ax.set_xticklabels(['Anchor end','Middle','Distal end']); tl[0].set_ha('left'); tl[2].set_ha('right')
ax.set_ylabel('Defense systems'); ax.legend(frameon=False,fontsize=6,loc='upper left',ncol=1)
p=S['pos']['All systems']
ax.text(0.98,0.97,f"distal 20%: {p['f_o']*100:.0f}% vs {p['f_n']*100:.0f}% expected\nP = {p['p_f']:.3f}; n = {p['ni']}",transform=ax.transAxes,ha='right',va='top',fontsize=6.2)
ax.set_title('Position along proviruses with an integration anchor',fontsize=7.5,loc='left',color=INK)
# c
ax=fig.add_subplot(gs[1,:]); lab(ax,'c',x=-0.045)
groups=[('All systems',None),('Direct targeting','Mechanism'),('Abi / TA','Mechanism'),('Other / unresolved','Mechanism'),
        ('Bacillota','Phylum'),('Pseudomonadota','Phylum'),('Other phyla','Phylum'),
        ('CheckV high/complete','Completeness'),('CheckV medium','Completeness'),('CheckV low','Completeness')]
xs=[];xl=[];k=0;prev=None
for g,cat in groups:
    if cat!=prev and prev is not None: k+=0.6
    prev=cat; d=S['pos'][g]
    ax.plot([k,k],[d['f_n']*100,d['f_o']*100],color=GRID,lw=2,zorder=1)
    ax.plot(k,d['f_n']*100,'o',ms=5,mfc='white',mec=INK,mew=0.9,zorder=2)
    ax.plot(k,d['f_o']*100,'o',ms=6,color=C['D'] if d['p_f']<0.05 else '#c7c7c7',mec='white',mew=0.8,zorder=3)
    star='*' if d['p_f']<0.05 else ''
    ax.text(k,max(d['f_o'],d['f_n'])*100+2.5,f"{star}\nn={d['ni']}",ha='center',va='bottom',fontsize=5.6,color=MUTED,linespacing=0.9)
    xs.append(k); xl.append(g); k+=1
ax.set_xticks(xs); ax.set_xticklabels(xl,rotation=25,ha='right',fontsize=6.3); ax.set_ylabel('% of systems in distal 20%\n(end opposite anchor)'); ax.set_ylim(0,52)
for name,lo,hi in [('Mechanism class',1,3),('Phylum',4,6),('CheckV completeness',7,9)]:
    ax.text((xs[lo]+xs[hi])/2,51,name,ha='center',va='top',fontsize=6.5,color=INK)
ax.plot([],[],'o',color=C['D'],ms=6,label='Observed (P < 0.05)'); ax.plot([],[],'o',color='#c7c7c7',ms=6,label='Observed (n.s.)'); ax.plot([],[],'o',mfc='white',mec=INK,ms=5,label='Expected')
ax.legend(frameon=False,fontsize=6,loc='lower left',bbox_to_anchor=(0.0,0.0))
ax.set_title('Distal placement by subgroup (unadjusted permutation P; systems are not independent replicates)',fontsize=7.5,loc='left',color=INK)
# d gene maps
ax=fig.add_subplot(gs[2,:]); lab(ax,'d',x=-0.045,y=1.2)
maps={}
for line in open('data/maps.txt'):
    k_,v=line.strip().split(':',1)
    maps[k_]=[(int(a),int(b),s,c,d=='D') for a,b,s,c,d in (t.split(',') for t in v.strip(';').split(';'))]
info={'A':('GCF_004103715.1 · Bacillota · complete genome · CheckV high-quality',3332019,3386816,[(3326749,3333365,'RM type I*'),(3334887,3340604,'RM type II + Lamassu')],True,[0.42,0.62]),
      'B':('GCA_002753565.1 · Pseudomonadota · draft · CheckV high-quality',43362,121636,[(38407,46520,'RM type I*'),(107414,114449,'RM type I + AbiD')],True,[0.42,0.42]),
      'C':('GCF_011009015.1 · Pseudomonadota · complete genome · CheckV high-quality',1279604,1400248,[(1360152,1363007,'Lamassu'),(1365582,1366986,'PhiR47 GevAT')],True,[0.42,0.62])}
W=122000; ys={'A':2.45,'B':1.25,'C':0.0}
for key in ['A','B','C']:
    title,a,b,systems,flip,lh=info[key]; y=ys[key]
    def X(p): return (b-p) if flip else (p-a)
    ax.plot([X(a)-6000 if not flip else -6000, (X(b)+6000) if not flip else X(a)+6000],[y,y],color='#bdbdbd',lw=0.8,zorder=0)
    ax.add_patch(Rectangle((0,y-0.2),b-a,0.4,fc='#f4f1ec',ec='none',zorder=0))
    for s,e,st,c,dd in maps[key]:
        x0,x1=sorted([X(s),X(e)]); L=x1-x0; fw=(st=='+')!=flip
        col=C['D'] if dd else C.get(c,'#ffffff')
        hl=min(L*0.35,900); h=0.26
        if fw: pts=[(x0,y-h/2),(x1-hl,y-h/2),(x1-hl,y-h/2-0.06),(x1,y),(x1-hl,y+h/2+0.06),(x1-hl,y+h/2),(x0,y+h/2)]
        else: pts=[(x1,y-h/2),(x0+hl,y-h/2),(x0+hl,y-h/2-0.06),(x0,y),(x0+hl,y+h/2+0.06),(x0+hl,y+h/2),(x1,y+h/2)]
        ax.add_patch(Polygon(pts,closed=True,fc=col,ec='#555555',lw=0.3,zorder=2))
    for (s,e,name),hh in zip(systems,lh):
        x0,x1=sorted([X(s),X(e)])
        if x0<0 or x1>b-a:  # crosses boundary: draw host-side part hatched
            hx0,hx1=(x0,0) if x0<0 else (b-a,x1)
            ax.add_patch(Rectangle((hx0,y-0.13),hx1-hx0,0.26,fc='white',ec=C['D'],hatch='/////',lw=0.6,zorder=2))
        ax.annotate(name,xy=((x0+x1)/2,y+0.2),xytext=((x0+x1)/2,y+hh),ha='center',fontsize=6,color=C['D'],arrowprops=dict(arrowstyle='-',lw=0.5,color=C['D']))
    ax.text(0,y-0.33,title,fontsize=6,color=MUTED,va='top')
ax.plot([0,10000],[-0.62,-0.62],color=INK,lw=1); ax.text(5000,-0.67,'10 kb',ha='center',va='top',fontsize=6)
ax.set_xlim(-9000,W+6000); ax.set_ylim(-0.8,3.25); ax.axis('off')
handles=[Rectangle((0,0),1,1,fc=C[k_],ec='#555555',lw=0.3) for k_ in ['D','I','S','L','R','O','U']]+[Rectangle((0,0),1,1,fc='white',ec=C['D'],hatch='/////')]
ax.legend(handles,['Defense system','Integrase/recombinase (anchor)','Capsid, tail, portal, terminase','Lysis','Replication / regulation','Other annotated','Unannotated','Beyond boundary'],
          frameon=False,fontsize=6,ncol=4,loc='upper center',bbox_to_anchor=(0.5,1.17))
ax.text(0.5,1.23,'Example proviruses carrying defense systems, drawn with the anchor end at left; shaded band = predicted provirus',transform=ax.transAxes,ha='center',fontsize=7.5,color=INK)
# e model
ax=fig.add_subplot(gs[3,:]); lab(ax,'e',x=-0.045); ax.axis('off'); ax.set_xlim(0,100); ax.set_ylim(6,42)
ax.plot([2,98],[20,20],color='#9e9e9e',lw=3,solid_capstyle='round',zorder=0)
ax.text(2,23.5,'Host chromosome',fontsize=6.5,color=MUTED)
ax.add_patch(FancyBboxPatch((22,15),56,10,boxstyle='round,pad=0.3,rounding_size=1.5',fc='#f4f1ec',ec='#8d8d8d',lw=0.8,zorder=1))
ax.add_patch(Rectangle((22.6,16.5),4.5,7,fc=C['I'],ec='none',zorder=2)); ax.text(24.85,20,'int',ha='center',va='center',color='white',fontsize=6,zorder=3)
ax.add_patch(Rectangle((28,16.5),12,7,fc=C['R'],ec='none',zorder=2)); ax.text(34,20,'regulation /\nreplication',ha='center',va='center',fontsize=5.8,zorder=3)
ax.add_patch(Rectangle((40.6,16.5),22,7,fc=C['S'],ec='none',zorder=2)); ax.text(51.6,20,'head · tail · lysis modules',ha='center',va='center',color='white',fontsize=6,zorder=3)
ax.add_patch(Rectangle((63.5,16.5),14,7,fc=C['D'],ec='none',zorder=2)); ax.text(70.5,20,'defense\nsystems',ha='center',va='center',color='white',fontsize=6,fontweight='bold',zorder=3)
for x0 in [12,84]:
    ax.add_patch(Rectangle((x0 if x0<50 else 88,17.5),5,5,fc='white',ec=C['D'],hatch='/////',lw=0.7,zorder=2))
ax.text(14.5,13,'flanking DNA ≤10 kb\n(pooled 2.7-fold)',ha='center',va='top',fontsize=6,color=INK)
ax.text(91,13,'flanking DNA ≤10 kb\n(pooled 2.7-fold)',ha='center',va='top',fontsize=6,color=INK)
ax.text(2,39,f"Fully contained systems: 2.4-fold  ·  {S['islands']['systems']} of 407 host-flanked systems co-localize with another system ({S['islands']['proviruses']} proviruses)",ha='left',fontsize=6.3,color=INK)
ax.annotate('',xy=(70.5,25.6),xytext=(70.5,28.6),arrowprops=dict(arrowstyle='->',lw=0.7,color=C['D']))
ax.text(70.5,29,'distal placement among 191 orientable systems: 26% vs 19% expected;\nclearest in Bacillota, Pseudomonadota and high/medium-completeness elements',ha='center',va='bottom',fontsize=6,color=INK)
ax.set_title('Schematic summary of pooled spatial associations (not a demonstrated mechanism or module order)',fontsize=7.5,loc='left',color=INK)
fig.savefig('Figure_8_v10_Spatial_Defense.png',dpi=300); fig.savefig('Figure_8_v10_Spatial_Defense.pdf')
print('ok')
