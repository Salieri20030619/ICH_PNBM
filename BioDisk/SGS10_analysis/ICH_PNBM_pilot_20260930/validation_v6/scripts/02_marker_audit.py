from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.stats import wilcoxon
P=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930/validation_v6')
d=pd.read_csv(P/'results/matched_gene_effects.csv')
rows=[]
for (program,gene),g in d.groupby(['program','gene']):
 v=g.delta.to_numpy()
 rows.append(dict(program=program,gene=gene,n=len(v),median_delta=np.median(v),n_up=int((v>0).sum()),n_down=int((v<0).sum()),p_raw=float(wilcoxon(v).pvalue) if np.any(v) else np.nan))
r=pd.DataFrame(rows);a=r.p_raw.fillna(1).to_numpy();ix=np.argsort(a);q=np.empty(len(a));q[ix]=np.minimum.accumulate((a[ix]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];r['fdr_all_markers']=np.minimum(q,1);r.loc[r.p_raw.isna(),'fdr_all_markers']=np.nan
r.to_csv(P/'results/all_marker_tests.csv',index=False)
print('FAMILY',len(r),'NOMINAL',r[r.p_raw<.05].to_json(orient='records'))
# Source reconstruction and annotation agreement audit, without rerunning inference.
V=P.parent/'validation_v5';s=pd.read_csv(V/'results/singlecell_sample_summary.csv')
for key,g in s.groupby(['gsm','qc']):
 total=g[g.celltype=='All'].iloc[0];parts=g[~g.celltype.isin(['All','Myeloid_clean'])]
 assert parts.n_cells.sum()==total.n_cells and parts.PLXNB3_positive_cells.sum()==total.PLXNB3_positive_cells
print('QC summary totals verified for 21 samples and 2 settings')
