from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.stats import wilcoxon
r=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930/validation_v3')
data=pd.read_csv(r/'results/sorted_patient_gene_values.csv');data['day']=data.bio_key.str.split('|').str[2]
targets=['ITGA2B','ITGB3','GRAMD1C','DEFA3','ALOX15','PRSS33','PLXNB3']
rows=[];pairs=[]
def bh(a):
 a=np.asarray(a);bad=~np.isfinite(a);a=np.where(bad,1.,a);o=np.argsort(a);q=np.empty(len(a));q[o]=np.minimum.accumulate((a[o]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];q=np.minimum(q,1);q[bad]=np.nan;return q
for cell in ['monos','neuts']:
 for g in targets:
  for stage in ['early','late']:
   d=data[(data.cell==cell)&(data.gene==g)&(data.stage==stage)]
   k=d.groupby(['patientid','day','source'])[['logcpm','cpm','hrs after symptom onset']].mean().unstack('source')
   if len(k) and all(('logcpm',x) in k for x in ['blood','hematoma']):
    k=k.dropna();k=k[abs(k[('hrs after symptom onset','blood')]-k[('hrs after symptom onset','hematoma')])<=6]
    vv=(k[('logcpm','hematoma')]-k[('logcpm','blood')]).groupby('patientid').mean()
    detect=np.mean(np.r_[k[('cpm','hematoma')],k[('cpm','blood')]]>=1) if len(k) else 0
    for pid,val in vv.items():pairs.append(dict(cell=cell,gene=g,stage=stage,patientid=pid,delta=val))
   else:vv=pd.Series(dtype=float);detect=0
   v=vv.to_numpy();p=wilcoxon(v).pvalue if len(v)>=3 and detect>=.2 and np.any(v) else np.nan
   rows.append(dict(cell=cell,gene=g,stage=stage,n_pairs=len(v),detectable_fraction=detect,median_change=np.median(v) if len(v) else np.nan,fraction_up=np.mean(v>0) if len(v) else np.nan,p_raw=p))
out=pd.DataFrame(rows);out['fdr_28']=bh(out.p_raw);out.to_csv(r/'results/strict_day_matched_tissue.csv',index=False)
pd.DataFrame(pairs).to_csv(r/'results/strict_day_matched_patient_deltas.csv',index=False)
print('STRICT',out[out.gene.isin(['PLXNB3','ITGB3'])].to_json(orient='records'))
a=pd.read_csv(r/'results/discovery_proxy_adjusted_all.csv.gz');e=pd.read_csv(r/'results/external_proxy_adjusted_all.csv.gz');cc=pd.read_csv(r/'results/target_proxy_correlations.csv')
print('PROXY',a[a.gene.isin(['ITGA2B','ITGB3','GRAMD1C','PLXNB3'])].to_json(orient='records'))
print('EXT_PROXY',e[e.gene.isin(['ITGA2B','ITGB3','GRAMD1C','PLXNB3'])].to_json(orient='records'))
print('CORR',cc[cc.gene.isin(['ITGA2B','ITGB3','PLXNB3'])].to_json(orient='records'))
