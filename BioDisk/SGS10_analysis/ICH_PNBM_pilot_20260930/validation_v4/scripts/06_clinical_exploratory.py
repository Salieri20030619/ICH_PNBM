from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.stats import spearmanr,rankdata
import statsmodels.api as sm
p=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930/validation_v4')
m=pd.read_csv(p.parent/'validation_v3/metadata/monos_collapsed_samples.csv',index_col=0)
x=pd.read_csv(p/'results/sorted_marker_logcpm.csv',index_col=0)
d=m.join(x.PLXNB3);d=d[(d.source=='hematoma')&(d.stage=='early')].sort_values('hrs after symptom onset').groupby('patientid').first()
cols=['modified rankin score at 1 year post-ich','initial ich volume','hrs after symptom onset']
rows=[]
for col in cols:
 a=d[['PLXNB3',col]].apply(pd.to_numeric,errors='coerce').dropna();r,pv=spearmanr(a.PLXNB3,a[col]);rows.append(dict(target=col,n=len(a),rho=r,p_raw=pv))
def bh(a):
 a=np.asarray(a);o=np.argsort(a);q=np.empty(len(a));q[o]=np.minimum.accumulate((a[o]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];return np.minimum(q,1)
o=pd.DataFrame(rows);o['fdr_3']=bh(o.p_raw);o.to_csv(p/'results/clinical_correlations_exploratory.csv',index=False)
d[['PLXNB3']+cols+['age']].to_csv(p/'results/clinical_patient_inputs.csv')
a=d[['PLXNB3',cols[0],'initial ich volume','age']].apply(pd.to_numeric,errors='coerce').dropna().rank()
fit=sm.OLS(a[cols[0]],sm.add_constant(a[['PLXNB3','initial ich volume','age']])).fit()
r=dict(n=len(a),model='rank(mRS) ~ rank(PLXNB3) + rank(initial ICH volume) + rank(age); exploratory sensitivity, earliest early hematoma per patient',target='PLXNB3',coefficient=fit.params.PLXNB3,p_raw=fit.pvalues.PLXNB3,ci=fit.conf_int().loc['PLXNB3'].tolist())
(p/'results/clinical_adjusted_sensitivity.json').write_text(json.dumps(r,indent=2))
print(o.to_json(orient='records'));print(r)
