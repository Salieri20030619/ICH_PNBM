from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.stats import wilcoxon,mannwhitneyu
b=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930');r=b/'validation_v2'
loo=pd.read_csv(r/'results/gene_leave_one_patient_out.csv')
print('ROBUST',loo[loo.loo_fraction_fdr_significant==1][['gene','contrast','fdr_limma_joint','loo_fdr_max']].to_json(orient='records'))
print('EARLY_LOO',loo[loo.contrast=='72h_vs_base'][['gene','fdr_limma_joint','loo_fdr_max','loo_fraction_fdr_significant']].to_json(orient='records'))
def bh(a):
 a=np.asarray(a);bad=~np.isfinite(a);a=np.where(bad,1,a);o=np.argsort(a);q=np.empty(len(a));q[o]=np.minimum.accumulate((a[o]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];q=np.minimum(q,1);q[bad]=np.nan;return q
d=pd.read_csv(r/'raw/external.xls.gz',sep='\t');d=d[~d.symbol.duplicated(False)&d.symbol.notna()].set_index('symbol')
panels=json.loads((b/'metadata/ich_prespecified_panels.json').read_text())
log=np.log2(d[[f'ICH{i}' for i in range(1,23)]]+1)
scores={name:log.loc[genes].mean(axis=0).to_numpy() for name,genes in panels.items()}
scores['presentation_minus_neutrophil']=scores['antigen_presentation']-scores['neutrophil_associated']
rows=[]
for name,a in scores.items():
 v=a[1::2]-a[::2]
 rows.append(dict(panel=name,n_pairs=11,median_change=np.median(v),fraction_up=np.mean(v>0),p_raw=wilcoxon(v).pvalue))
o=pd.DataFrame(rows);o['fdr_four_prespecified_panels']=bh(o.p_raw);o.to_csv(r/'results/original_panels_external.csv',index=False)
print('ORIGINAL_PANELS_EXTERNAL',o.to_json(orient='records'))
d=pd.read_excel(b/'raw/ijms-3877263-supplementary.xlsx',header=1);d.columns=d.columns.str.strip();features=d.columns[3:18]
c=d[d.Database==1].sort_values('Time from admission, days').drop_duplicates('Patients ID').set_index('Patients ID')
js=[]
for pid in c.index:
 cc=c.drop(index=pid);ps=[]
 for m in features:
  x=cc.loc[cc.group==1,m].dropna().to_numpy();y=cc.loc[cc.group==0,m].dropna().to_numpy()
  ps.append(mannwhitneyu(x,y,method='asymptotic').pvalue if len(x)>1 and len(y)>1 and np.ptp(np.r_[x,y])>0 else np.nan)
 qs=bh(ps)
 for m in ['p-HPhLA','PhLA']:
  j=list(features).index(m);js.append(dict(excluded_patient=int(pid),feature=m,p_raw=ps[j],fdr=qs[j]))
j=pd.DataFrame(js);j.to_csv(r/'results/chronic_infection_loo.csv',index=False)
print('CHRONIC_LOO',j.groupby('feature').fdr.agg(['min','max']).to_json())
