from pathlib import Path
import pandas as pd,numpy as np,subprocess,json
from scipy.stats import t as tdist
b=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930');r=b/'validation_v3'
targets=['ITGA2B','ITGB3','GRAMD1C','DEFA3','ALOX15','PRSS33','PLXNB3']
plt=['PF4','PPBP','GP1BA','GP9','TUBB1']
rs="""suppressPackageStartupMessages(library(edgeR))
a<-commandArgs(trailingOnly=TRUE);x<-as.matrix(read.csv(gzfile(a[1]),row.names=1,check.names=FALSE));d<-calcNormFactors(DGEList(x))
gg<-c("ITGA2B","ITGB3","GRAMD1C","DEFA3","ALOX15","PRSS33","PLXNB3","PF4","PPBP","GP1BA","GP9","TUBB1")
write.csv(cpm(d,log=TRUE,prior.count=1)[intersect(gg,rownames(x)),,drop=FALSE],a[2])
"""
rp=r/'scripts/sorted_proxy_normalize.R';rp.write_text(rs)
rows=[]
base_test=pd.read_csv(r/'results/strict_day_matched_tissue.csv')
def bh(a):
 a=np.asarray(a);bad=~np.isfinite(a);a=np.where(bad,1.,a);o=np.argsort(a);q=np.empty(len(a));q[o]=np.minimum.accumulate((a[o]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];q=np.minimum(q,1);q[bad]=np.nan;return q
for cell in ['monos','neuts']:
 lp=r/'results'/(cell+'_target_and_platelet_logcpm.csv')
 z=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(rp),str(r/'raw'/(cell+'_collapsed_counts.csv.gz')),str(lp)],capture_output=True,text=True)
 assert z.returncode==0,z.stderr
 x=pd.read_csv(lp,index_col=0)
 m=pd.read_csv(r/'metadata'/(cell+'_collapsed_samples.csv'),dtype={'patientid':str,'day of trial':str}).set_index('bio_key')
 m['platelet_proxy']=x.loc[plt,m.index].mean(axis=0).to_numpy()
 for gene in targets:
  for stage in ['early','late']:
   row=dict(cell=cell,gene=gene,stage=stage,n_pairs=0,adjusted_effect=np.nan,p_raw=np.nan)
   eligibility=base_test[(base_test.cell==cell)&(base_test.gene==gene)&(base_test.stage==stage)].iloc[0]
   if gene not in x.index or eligibility.detectable_fraction<.2:rows.append(row);continue
   mm=m[m.stage==stage].copy();mm['target']=x.loc[gene,mm.index].to_numpy()
   d=mm.groupby(['patientid','day of trial','source'])[['target','platelet_proxy','hrs after symptom onset']].mean().unstack('source').dropna()
   d=d[abs(d[('hrs after symptom onset','blood')]-d[('hrs after symptom onset','hematoma')])<=6]
   diffs=pd.DataFrame({'target':d[('target','hematoma')]-d[('target','blood')],'proxy':d[('platelet_proxy','hematoma')]-d[('platelet_proxy','blood')]}).groupby('patientid').mean()
   v=diffs.target.to_numpy();z=diffs.proxy.to_numpy();n=len(v);X=np.column_stack([np.ones(n),z])
   if n>3 and np.linalg.matrix_rank(X)==2:
    coef=np.linalg.lstsq(X,v,rcond=None)[0];res=v-X@coef;cov=(res@res/(n-2))*np.linalg.inv(X.T@X);se=np.sqrt(cov[0,0]);pv=2*tdist.sf(abs(coef[0]/se),n-2) if se>0 else np.nan
    row.update(n_pairs=n,adjusted_effect=coef[0],p_raw=pv,proxy_rho=diffs.corr(method='spearman').loc['target','proxy'],proxy_range_min=z.min(),proxy_range_max=z.max())
   rows.append(row)
o=pd.DataFrame(rows);o['fdr_28']=bh(o.p_raw);o.to_csv(r/'results/sorted_platelet_proxy_sensitivity.csv',index=False)
print(o[o.gene.isin(['PLXNB3','ITGB3'])].to_json(orient='records'))
(r/'logs/sorted_proxy.exit').write_text('0')
