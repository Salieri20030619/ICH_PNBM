from pathlib import Path
import pandas as pd,numpy as np,itertools,json,re,subprocess
from scipy.stats import rankdata,mannwhitneyu
b=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930');r=b/'validation_v2'
d=pd.read_excel(b/'raw/ijms-3877263-supplementary.xlsx',header=1);d.columns=d.columns.str.strip()
d['patient_key']=d.Database.astype(str)+'_'+d['Patients ID'].astype(str)
features=d.columns[3:18].tolist()
def bh(v):
 a=np.asarray(v);missing=~np.isfinite(a);a=np.where(missing,1.,a);o=np.argsort(a);q=np.empty(len(a));q[o]=np.minimum.accumulate((a[o]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];q=np.minimum(q,1);q[missing]=np.nan;return q
acute=d[d.Database==2].groupby('patient_key')[features+['group']].median()
chronic=d[d.Database==1].sort_values('Time from admission, days').drop_duplicates('patient_key').set_index('patient_key')
rows=[]
for cohort,p in [('acute_postoperative',acute),('chronic_brain_injury',chronic)]:
 for m in features:
  x=p.loc[p.group==1,m].dropna().to_numpy(float);y=p.loc[p.group==0,m].dropna().to_numpy(float)
  if min(len(x),len(y))<2 or np.ptp(np.r_[x,y])==0:
   rows.append(dict(cohort=cohort,feature=m,n_infected=len(x),n_noninfected=len(y),p_raw=np.nan,method='not_testable_missing_or_constant'));continue
  ranks=rankdata(np.r_[x,y]);n=len(x);obs=abs(ranks[:n].sum()-n*ranks.mean())
  if cohort=='acute_postoperative':
   v=np.array([ranks[list(c)].sum() for c in itertools.combinations(range(len(ranks)),n)])
   pv=float(np.mean(abs(v-n*ranks.mean())>=obs-1e-9));method='exact_label_permutation_rank_sum'
  else:
   pv=float(mannwhitneyu(x,y,alternative='two-sided',method='asymptotic').pvalue);method='two_sided_MW_tie_corrected'
  auc=float(((x[:,None]>y).sum()+.5*(x[:,None]==y).sum())/(len(x)*len(y)))
  rows.append(dict(cohort=cohort,feature=m,n_infected=len(x),n_noninfected=len(y),median_infected=np.median(x),median_noninfected=np.median(y),
 median_ratio=np.median(x)/np.median(y) if np.median(y)>0 else np.nan,p_raw=pv,auc_descriptive=auc,method=method))
out=pd.DataFrame(rows);out['fdr_within_cohort_15']=out.groupby('cohort').p_raw.transform(lambda x:bh(x));out['fdr_all_30_sensitivity']=bh(out.p_raw)
out.to_csv(r/'results/infection_all_analytes.csv',index=False)
text=(r/'metadata/paired_infection_S1.txt').read_text()
rec={};ids=[]
for line in text.splitlines():
 toks=line.strip().split()
 if toks and toks[0]=='Sample' and len(toks)>1 and all(t.isdigit() for t in toks[1:]):
  ids=[int(t) for t in toks[1:]]
  for k in ids:rec.setdefault(k,{})
 elif toks and toks[0].endswith('_csf') and len(toks)==len(ids)+1:
  try:
   for k,v in zip(ids,toks[1:]):rec[k][toks[0][:-4]]=float(v)
  except ValueError:pass
new=pd.DataFrame.from_dict(rec,orient='index')
met=features[:11];matches=[]
for sample,row in new.iterrows():
 cols=[x for x in met if x in row.index and pd.notna(row[x])]
 if len(cols)<8:continue
 mask=np.isclose(d[cols].to_numpy(float),row[cols].to_numpy(float),rtol=1e-5,atol=.0001).all(axis=1)
 for _,z in d.loc[mask].iterrows():
  matches.append(dict(external_sample=sample,source_sample=str(z['Sample ID']),source_patient=z.patient_key,matched_analytes=len(cols),cohort=int(z.Database)))
pd.DataFrame(matches).to_csv(r/'results/infection_external_overlap.csv',index=False)
summary={'all_analytes':features,'acute_fdr_hits':out[(out.cohort=='acute_postoperative')&(out.fdr_within_cohort_15<.05)].to_dict('records'),
'chronic_fdr_hits':out[(out.cohort=='chronic_brain_injury')&(out.fdr_within_cohort_15<.05)].to_dict('records'),
'pdf_samples_extracted':len(new),'pdf_matched_samples':len(set(x['external_sample'] for x in matches)),
'pdf_matched_patients':len(set(x['source_patient'] for x in matches)),
'warning':'Chronic sensitivity cohort is not acute postoperative cohort. External paired-CSF/serum source overlapping exact multianalyte profiles cannot be called independent replication.'}
(r/'results/infection_v2_summary.json').write_text(json.dumps(summary,indent=2))
(r/'logs/infection_v2.exit').write_text('0')
print(json.dumps(summary,indent=2),flush=True)
