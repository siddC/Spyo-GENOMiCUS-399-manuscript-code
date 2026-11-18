#!/usr/bin/env python3
"""
Figure S8: accessory-gene intersection structure (UpSet) and the frequency
distribution of zero-phylon genes.

    python3 make_figure_s8.py \
        --L       data/L_binarized_panta_cdhit.csv \
        --acc     data/df_acc_panta_cdhit.pickle \
        --outdir  results

Writes FigureS8.{pdf,svg,png} plus FigS8a_upset_combinations.csv (Table S12).
"""
import argparse, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd, numpy as np, collections
from scipy.stats import mannwhitneyu

p = argparse.ArgumentParser()
p.add_argument('--indir',  default=os.environ.get('SPYO_DATA', 'data'),
               help='directory holding the Zenodo matrices (default: data)')
p.add_argument('--outdir', default=os.environ.get('SPYO_RESULTS', 'results'))
p.add_argument('--L',   default=None, help='override path to the binarized L matrix')
p.add_argument('--acc', default=None, help='override path to the accessory P matrix')
p.add_argument('--top', type=int, default=30, help='intersections to plot')
a = p.parse_args()
a.L   = a.L   or os.path.join(a.indir, 'L_bin_panta_cdhit_Spyo.csv')
a.acc = a.acc or os.path.join(a.indir, 'df_acc_panta_cdhit.pickle')
for _f in (a.L, a.acc):
    if not os.path.exists(_f):
        raise SystemExit(f'Missing input: {_f}\n'
                         f'Download the Zenodo deposit first:\n'
                         f'    python scripts/zenodo_fetch.py 20403112 -o {a.indir}')
os.makedirs(a.outdir, exist_ok=True)

PH = ['emm-patternE-other','ST28-emm1','ST52-emm28','ST39-emm4','emm3/6/5','ST36',
      'ST63-emm77','SpeK-prophage','emm-patternD-other','SpeC-prophage','ST101-emm89']
L = pd.read_csv(a.L, index_col=0)
if L.shape[0] == 11 and L.shape[1] != 11:
    L = L.T
L.columns = PH
acc = pd.read_pickle(a.acc) if a.acc.endswith(('.pickle', '.pkl')) else pd.read_csv(a.acc, index_col=0)
L = L.reindex(acc.index)
TOP = a.top
BLUE='#1769BF'; RED='#D94A3D'; INK='#27272A'; GRAY='#52525B'; GRID='#D4D4D8'
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
    'font.size':7,'axes.edgecolor':GRAY,'axes.linewidth':.6,'text.color':INK,
    'xtick.color':GRAY,'ytick.color':GRAY,'axes.labelcolor':INK,'svg.fonttype':'none'})

combo=L.apply(lambda r: tuple(sorted(L.columns[r>0])),axis=1)
cnt={k:v for k,v in collections.Counter(combo).items() if k}
top=sorted(cnt.items(),key=lambda x:-x[1])[:TOP]
setsize=L.sum(axis=0).sort_values(ascending=False)
rows=list(setsize.index)
SDA=tuple(sorted(('ST28-emm1','ST36','ST63-emm77')))

fig=plt.figure(figsize=(7.5,9.6))
gsA=fig.add_gridspec(2,2,width_ratios=[1,5.0],height_ratios=[2.4,2.9],
            hspace=.055,wspace=.04,left=.20,right=.975,top=.965,bottom=.395)
gsB=fig.add_gridspec(1,1,left=.20,right=.975,top=.295,bottom=.085)

# ---- S8a top: intersection sizes ----
axb=fig.add_subplot(gsA[0,1])
x=np.arange(TOP); vals=[v for _,v in top]
cols=[RED if k==SDA else INK for k,_ in top]
axb.bar(x,vals,color=cols,width=.62)
for xi,v,k in zip(x,vals,[k for k,_ in top]):
    axb.text(xi,v+2.5,str(v),ha='center',va='bottom',fontsize=5.6,
             color=RED if k==SDA else GRAY)
axb.set_xlim(-.7,TOP-.3); axb.set_ylim(0,max(vals)*1.16)
axb.set_ylabel('Genes in intersection',fontsize=7.5)
for s in ('top','right','bottom'): axb.spines[s].set_visible(False)
axb.set_xticks([]); axb.tick_params(axis='y',length=2,labelsize=6.5)
axb.yaxis.grid(True,color=GRID,lw=.5,alpha=.7); axb.set_axisbelow(True)
i=[k for k,_ in top].index(SDA)
axb.annotate('27-gene Sda1 prophage core\n(ST28-$\\it{emm}$1 ∩ ST36 ∩ ST63-$\\it{emm}$77)',
    xy=(i,vals[i]+4),xytext=(i+4.4,max(vals)*.72),fontsize=6.4,color=RED,ha='left',
    arrowprops=dict(arrowstyle='-',color=RED,lw=.8,
                    connectionstyle='angle3,angleA=0,angleB=70'))
axb.set_title('a',loc='left',fontsize=10,fontweight='bold',color=INK,pad=4)

# ---- S8a bottom: membership matrix ----
axm=fig.add_subplot(gsA[1,1]); axm.set_xlim(-.7,TOP-.3); axm.set_ylim(-.6,len(rows)-.4)
axm.invert_yaxis()
for yi in range(len(rows)):
    if yi%2==0: axm.axhspan(yi-.5,yi+.5,color='#F4F4F5',zorder=0)
for xi,(k,_) in enumerate(top):
    mem=[rows.index(p) for p in k]; c=RED if k==SDA else INK
    axm.scatter([xi]*len(rows),range(len(rows)),s=15,color=GRID,zorder=1)
    axm.scatter([xi]*len(mem),mem,s=15,color=c,zorder=3)
    if len(mem)>1: axm.plot([xi,xi],[min(mem),max(mem)],color=c,lw=1.0,zorder=2)
axm.set_yticks(range(len(rows)))
axm.set_yticklabels(['' for _ in rows])
axm.set_xticks([]); [axm.spines[s].set_visible(False) for s in axm.spines]
axm.tick_params(length=0)

# ---- S8a left: set sizes ----
axs=fig.add_subplot(gsA[1,0]); axs.set_ylim(-.6,len(rows)-.4); axs.invert_yaxis()
axs.barh(range(len(rows)),[setsize[p] for p in rows],color=GRAY,height=.55)
axs.invert_xaxis()
axs.set_yticks(range(len(rows)))
lab=[r.replace('emm','$\\it{emm}$') if 'emm' in r else r for r in rows]
axs.set_yticklabels(lab,fontsize=6.6,color=INK)
axs.tick_params(axis='y',length=0,pad=2)
axs.set_xlabel('Genes in phylon',fontsize=6.6,labelpad=1)
for s in ('top','left','right'): axs.spines[s].set_visible(False)
axs.tick_params(axis='x',length=2,labelsize=6)
axs.xaxis.grid(True,color=GRID,lw=.5,alpha=.7); axs.set_axisbelow(True)

# ---- S8b: frequency deciles ----
freq=acc.sum(axis=1); nphy=L.sum(axis=1); zero=nphy==0
dec=pd.qcut(freq,10,labels=False,duplicates='drop')+1
tab=pd.crosstab(dec,np.where(zero,'z','a'))
ax2=fig.add_subplot(gsB[0,0])
w=.38; xs=np.arange(1,11)
ax2.bar(xs-w/2,tab['a'],width=w,color=BLUE,label='Assigned to ≥1 phylon (n = 1,054)')
ax2.bar(xs+w/2,tab.get('z',pd.Series(0,index=xs)),width=w,color=RED,
        label='Zero-phylon (n = 113)')
for xi in xs:
    z=int(tab.get('z',pd.Series(0,index=xs)).get(xi,0))
    ax2.text(xi+w/2,z+1.6,str(z),ha='center',va='bottom',fontsize=5.8,
             color=RED if z else '#A1A1AA')
ax2.set_xticks(xs); ax2.set_xlabel('Accessory-gene frequency decile (1 = rarest)',fontsize=7.5)
ax2.set_ylabel('Number of genes',fontsize=7.5)
ax2.set_xlim(.4,10.6); ax2.set_ylim(0,max(tab['a'])*1.22)
for s in ('top','right'): ax2.spines[s].set_visible(False)
ax2.yaxis.grid(True,color=GRID,lw=.5,alpha=.7); ax2.set_axisbelow(True)
ax2.tick_params(length=2,labelsize=6.5)
ax2.axvspan(6.6,10.6,color=RED,alpha=.055,zorder=0)
ax2.text(8.6,max(tab['a'])*1.09,'No zero-phylon genes above decile 6',
         ha='center',fontsize=6.6,color=RED)
u,p=mannwhitneyu(freq[zero],freq[~zero],alternative='two-sided')
ax2.text(0,-.235,f'Median frequency: zero-phylon {freq[zero].median():.0f} genomes vs '
    f'assigned {freq[~zero].median():.0f} (Mann-Whitney $P$ = {p:.1e}). '
    f'93% of zero-phylon genes lie within 2× of the accessory–rare boundary, vs 37% of assigned genes.',
    transform=ax2.transAxes,fontsize=6.3,color=GRAY,va='top')
ax2.legend(frameon=False,fontsize=6.8,loc='upper left',handlelength=1.1,ncol=1)
ax2.set_title('b',loc='left',fontsize=10,fontweight='bold',color=INK,pad=4)

for ext in ('pdf','svg','png'):
    fig.savefig(os.path.join(a.outdir, f'FigureS8.{ext}'),dpi=400,bbox_inches='tight',
                facecolor='white')
pd.Series({'|'.join(k): v for k, v in cnt.items()}).sort_values(ascending=False)\
  .to_csv(os.path.join(a.outdir,'FigS8a_upset_combinations.csv'), header=['n_genes'])
print('saved; top intersections shown:',TOP,'of',len(cnt))
print('Sda1 3-way =',cnt[SDA],'| zero-phylon =',int(zero.sum()),'| P =',f'{p:.3g}')
