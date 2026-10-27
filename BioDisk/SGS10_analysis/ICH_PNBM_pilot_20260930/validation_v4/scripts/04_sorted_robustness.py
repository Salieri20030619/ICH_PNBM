from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.stats import wilcoxon,spearmanr,binomtest
p=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930/validation_v4');v=p.parent/'validation_v3'
m=pd.read_csv(v/'metadata/monos_collapsed_samples.csv',index_col=0)
c=pd.read_csv(p/'results/sorted_marker_cpm.csv',index_col=0);log=pd.read_csv(p/'results/sorted_marker_logcpm.csv',index_col=0)
def bh(a):
 a=np.asarray(a,float);bad=~np.isfinite(a);a=np.where(bad,1.,a);o=np.argsort(a);q=np.empty(len(a));q[o]=np.minimum.accumulate((a[o]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];q=np.minimum(q,1);q[bad]=np.nan;return q
panels={'astro':['GFAP','AQP4','SLC1A3','ALDH1L1','HEPACAM'],'oligo':['PLP1','MBP','MOG','SOX10'],'opc':['OLIG1','OLIG2','PDGFRA','CSPG4','GPR17'],'platelet':['PF4','PPBP','GP1BA','GP9','TUBB1']}
z=pd.DataFrame(index=c.index)
for name,genes in panels.items():
 genes=[g for g in genes if g in c]
 z[name]=np.log2(c[genes]+1).mean(axis=1)
z['PLXNB3']=log.PLXNB3
d=m.join(z);d=d[d.stage=='early']
k=d.groupby(['patientid','day of trial','source'])[list(z.columns)+['hrs after symptom onset']].mean().unstack('source').dropna()
k=k[abs(k[('hrs after symptom onset','blood')]-k[('hrs after symptom onset','hematoma')])<=6]
k.to_csv(p/'results/sorted_strict_day_markers.csv')
paired={}
for source in ['blood','hematoma']:
 paired[source]=k.xs(source,axis=1,level=1).groupby('patientid').mean()
paired['delta']=paired['hematoma']-paired['blood']
rows=[]
for source,a in paired.items():
 a.to_csv(p/'results'/('sorted_patient_'+source+'.csv'))
 for panel in panels:
  r,pv=spearmanr(a.PLXNB3,a[panel]);rows.append(dict(source=source,panel=panel,n=len(a),rho=r,p_raw=pv))
o=pd.DataFrame(rows);o['fdr_12']=bh(o.p_raw);o.to_csv(p/'results/sorted_source_correlations.csv',index=False)
# Prespecified descriptive robustness, not new primary evidence.
rows=[]
for pseudo in [.1,1,5]:
 q=m[['patientid','day of trial','source','stage','hrs after symptom onset']].join(c.PLXNB3.rename('cpm'))
 q=q[q.stage=='early'];q['transformed']=np.log2(q.cpm+pseudo)
 a=q.groupby(['patientid','day of trial','source'])[['transformed','hrs after symptom onset']].mean().unstack('source').dropna()
 a=a[abs(a[('hrs after symptom onset','blood')]-a[('hrs after symptom onset','hematoma')])<=6]
 dif=(a[('transformed','hematoma')]-a[('transformed','blood')]).groupby('patientid').mean()
 rows.append(dict(pseudocount=pseudo,n=len(dif),median_delta=dif.median(),up=int((dif>0).sum()),p_raw=wilcoxon(dif).pvalue))
pd.DataFrame(rows).to_csv(p/'results/pseudocount_sensitivity.csv',index=False)
orig=pd.read_csv(v/'results/strict_day_matched_tissue.csv')
pv=pd.read_csv(v/'results/strict_day_matched_patient_deltas.csv')
target=pv[(pv.cell=='monos')&(pv.gene=='PLXNB3')&(pv.stage=='early')]
loos=[]
for exclude in target.patientid:
 rr=[]
 for a in orig.to_dict('records'):
  x=pv[(pv.cell==a['cell'])&(pv.gene==a['gene'])&(pv.stage==a['stage'])&(pv.patientid!=exclude)].delta
  a['p_raw']=wilcoxon(x).pvalue if len(x)>=3 and np.isfinite(a['p_raw']) and np.any(x) else np.nan
  rr.append(a)
 rr=pd.DataFrame(rr);rr['fdr_28']=bh(rr.p_raw)
 t=rr[(rr.cell=='monos')&(rr.gene=='PLXNB3')&(rr.stage=='early')].iloc[0]
 loos.append(dict(excluded=int(exclude),p_raw=t.p_raw,fdr_28=t.fdr_28))
pd.DataFrame(loos).to_csv(p/'results/sorted_leave_one_patient_out.csv',index=False)
print('CORR',o.to_json(orient='records'))
print('PSEUDO',json.dumps(rows));print('LOO',pd.DataFrame(loos).describe().to_json())
# patient averaged raw CPM and individual sample count evidence.
audit=m.join(c[['PLXNB3']].rename(columns={'PLXNB3':'PLXNB3_CPM'})).join(pd.read_csv(p/'results/sorted_marker_counts.csv',index_col=0)[['PLXNB3']].rename(columns={'PLXNB3':'PLXNB3_raw_count'}))
audit.to_csv(p/'results/sorted_plxnb3_sample_audit.csv')
print('CPM',audit.groupby(['source','stage']).PLXNB3_CPM.agg(['count','median','min','max']).to_json())
