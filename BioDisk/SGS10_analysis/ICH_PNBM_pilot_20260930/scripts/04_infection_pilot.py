from pathlib import Path
import pandas as pd,numpy as np,json,itertools,hashlib,sys
from scipy.stats import rankdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930')
d=pd.read_excel(r/'raw/ijms-3877263-supplementary.xlsx',header=1)
d.columns=d.columns.str.strip()
d['patient_key']=d.Database.astype(str)+'_'+d['Patients ID'].astype(str)
assert len(d)==96
d['source_excel_row']=np.arange(len(d))+3
d[d['Sample ID'].duplicated(False)].to_csv(r/'results/infection_duplicate_sample_ids.csv',index=False)
d.to_csv(r/'metadata/infection_samples.csv',index=False)
d.isna().mean().rename('missing_fraction').to_csv(r/'results/infection_missingness.csv')
markers=['p-HPhLA','PhLA','3ILA','IL_6']
a=d[d.Database==2].copy()
assert a['Sample ID'].is_unique
assert a.groupby('patient_key').group.nunique().max()==1
p=a.groupby('patient_key')[markers+['group','Positive CSF culture']].median()
p.to_csv(r/'results/infection_patient_level.csv')
rng=np.random.default_rng(20260930)
def auc(x,y):return float(((x[:,None]>y).sum()+0.5*(x[:,None]==y).sum())/(len(x)*len(y)))
rows=[]
for m in markers:
 x=p.loc[p.group==1,m].to_numpy();y=p.loc[p.group==0,m].to_numpy()
 values=np.r_[x,y];ranks=rankdata(values);n=len(x);expected=n*ranks.mean()
 obs=abs(ranks[:n].sum()-expected)
 perm=np.array([ranks[list(c)].sum() for c in itertools.combinations(range(len(values)),n)])
 pv=float(np.mean(abs(perm-expected)>=obs-1e-9))
 boot=np.array([auc(rng.choice(x,len(x)),rng.choice(y,len(y))) for _ in range(2000)])
 drop=[auc(np.delete(x,i),y) for i in range(len(x))]+[auc(x,np.delete(y,i)) for i in range(len(y))]
 rows.append(dict(marker=m,n_infected=len(x),n_noninfected=len(y),median_infected=np.median(x),median_noninfected=np.median(y),
 median_ratio=np.median(x)/np.median(y),p_raw=pv,auc_descriptive=auc(x,y),auc_bootstrap_lo=np.quantile(boot,.025),auc_bootstrap_hi=np.quantile(boot,.975),auc_leave_one_out_min=min(drop),auc_leave_one_out_max=max(drop)))
out=pd.DataFrame(rows)
order=np.argsort(out.p_raw);q=np.empty(len(out));q[order]=np.minimum.accumulate((out.p_raw.to_numpy()[order]*len(out)/np.arange(1,len(out)+1))[::-1])[::-1];out['fdr_bh']=np.minimum(q,1)
out['test_family']='4_preselected_markers_in_17_postneurosurgical_patients'
out.to_csv(r/'results/infection_marker_pilot.csv',index=False)
first=a.drop_duplicates('patient_key').set_index('patient_key')
sens=[]
for m in markers:
 x=first.loc[first.group==1,m].to_numpy();y=first.loc[first.group==0,m].to_numpy()
 sens.append(dict(marker=m,auc_first_source_row=auc(x,y),note='Source-row sensitivity only; timestamps unavailable; not earliest clinical sample.'))
pd.DataFrame(sens).to_csv(r/'results/infection_sensitivity.csv',index=False)
fig,axs=plt.subplots(1,4,figsize=(12,3.8))
for ax,m in zip(axs,markers):
 for g,col in [(0,'#3987a3'),(1,'#c5664b')]:
  z=p.loc[p.group==g,m].to_numpy()
  ax.scatter(np.full(len(z),g)+rng.uniform(-.08,.08,len(z)),z,color=col,s=30,alpha=.85)
  ax.hlines(np.median(z),g-.18,g+.18,color=col,lw=2)
 ax.set_yscale('log');ax.set_xticks([0,1],['Noninfected\nn=7','Infected\nn=10'])
 ax.set_title(m);ax.spines[['top','right']].set_visible(False)
 ax.set_ylabel('IL-6 (pg/mL, log scale)' if m=='IL_6' else 'Concentration (nmol/L, log scale)')
fig.suptitle('Post-neurosurgical CSF pilot: patient medians, exploratory')
fig.tight_layout();fig.savefig(r/'results/infection_pilot.png',dpi=160);plt.close(fig)
hem=d[d['Type of Diagnosis']=='Haemorrage']
summary=dict(source='10.3390/ijms262110522; correction 10.3390/ijms262411856',
 n_rows=len(d),reported_patients=53,observed_patient_ids=d['Patients ID'].nunique(),observed_cohort_patient_keys=d.patient_key.nunique(),
 cohort_counts=d.groupby('Database').agg(samples=('Sample ID','size'),patient_keys=('patient_key','nunique')).reset_index().to_dict('records'),
 acute_patients=len(p),acute_infected=int((p.group==1).sum()),acute_noninfected=int((p.group==0).sum()),
 acute_culture_positive_patients=int((a.groupby('patient_key')['Positive CSF culture'].max()>0).sum()),
 acute_sample_times_missing=int(a['Time from admission, days'].isna().sum()),
 hemorrhage_patient_keys=hem.patient_key.nunique(),hemorrhage_infected_patient_keys=hem.loc[hem.group==1,'patient_key'].nunique(),
 hemorrhage_noninfected_patient_keys=hem.loc[hem.group==0,'patient_key'].nunique(),
 label_switch_patients=d.groupby('patient_key').group.nunique().loc[lambda x:x>1].index.tolist(),
 il6_5000_samples=int((a.IL_6==5000).sum()),marker_results=out.to_dict('records'),
 limits=['Exploratory reanalysis of published data; no independent validation.',
 'Acute cohort consists of tumor surgery patients; chronic cohort is a different population.',
 'One patient ID is reused across cohorts with different age and diagnosis; cohort-specific keys required.',
 '53 reported patients versus 52 cohort-specific keys: unresolved source discrepancy.',
 'Two chronic-cohort sample IDs are reused on different collection days; preserved with source row identifiers; not part of acute analysis.',
 'No acute sampling time; medians summarize collected samples, not a prospective diagnostic timepoint.',
 'Only one culture-positive acute patient; infection labels also rely on clinical/CSF criteria.',
 'IL6 value 5000 is repeated; upper-limit censoring cannot be ruled out.',
 'Correction concerns reference concentrations; no use of reference thresholds in this pilot.',
 'One chronic tumor patient changes infection label; not evidence of a hemorrhage infection-before/after cohort.'])
(r/'results/infection_summary.json').write_text(json.dumps(summary,indent=2))
(r/'logs/infection_pilot.exit').write_text('0')
print(json.dumps(summary,indent=2),flush=True)
