from pathlib import Path
import pandas as pd,numpy as np,json,subprocess
from scipy.stats import wilcoxon
b=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930');r=b/'validation_v3'
targets=['ITGA2B','ITGB3','GRAMD1C','DEFA3','ALOX15','PRSS33','PLXNB3']
plan={'targets':targets,'time':'early <96h vs late 96-171h after onset','QC':'fastq_frag_count>=4e6, gene_count>=9000 monos/8500 neuts; exclude missing hours and >171h; healthy controls excluded',
'technical_replicates':'sum raw counts for same patient, cell type, source, trial day and onset hours before edgeR TMM normalization',
'analysis':'mean logCPM across days within patient/stage/compartment; paired Wilcoxon across patients. Temporal family 7x2x2=28. Tissue family 7x2x2=28. Source/tissue support, not matched day7 whole-blood replication.',
'detectability':'CPM>=1 in at least20% of samples in each tested stratum; otherwise no inferential P'}
(r/'metadata/sorted_validation_plan.json').write_text(json.dumps(plan,indent=2))
meta=pd.read_csv(r/'metadata/sorted_samples.csv',dtype=str)
for c in ['hrs after symptom onset','fastq_frag_count','gene_count']:meta[c]=pd.to_numeric(meta[c],errors='coerce')
rs="""suppressPackageStartupMessages(library(edgeR))
args<-commandArgs(trailingOnly=TRUE)
x<-as.matrix(read.csv(gzfile(args[1]),row.names=1,check.names=FALSE))
d<-calcNormFactors(DGEList(x))
genes<-c("ITGA2B","ITGB3","GRAMD1C","DEFA3","ALOX15","PRSS33","PLXNB3")
genes<-intersect(genes,rownames(x))
write.csv(cpm(d,log=TRUE,prior.count=1)[genes,,drop=FALSE],args[2])
write.csv(cpm(d,log=FALSE)[genes,,drop=FALSE],args[3])
"""
rp=r/'scripts/sorted_normalize.R';rp.write_text(rs)
allrows=[];counts=[]
for cell in ['monos','neuts']:
 x=pd.read_csv(r/'raw'/(cell+'_counts.csv.gz'),index_col=0)
 m=meta[meta['cell type']==cell].set_index('matrix_id').loc[x.columns].rename_axis('matrix_id').reset_index()
 keep=m.source.isin(['blood','hematoma'])&m['hrs after symptom onset'].between(0,171)&(m.fastq_frag_count>=4e6)&(m.gene_count>=(9000 if cell=='monos' else 8500))
 m=m[keep].copy();x=x[m.matrix_id]
 m['bio_key']=m.patientid+'|'+m.source+'|'+m['day of trial']+'|'+m['hrs after symptom onset'].astype(str)
 xx=x.T.groupby(m.set_index('matrix_id').bio_key).sum().T
 mm=m.drop_duplicates('bio_key').set_index('bio_key').loc[xx.columns].copy()
 assert len(xx.columns)==len(mm)
 mm['stage']=np.where(mm['hrs after symptom onset']<96,'early','late');mm.to_csv(r/'metadata'/(cell+'_collapsed_samples.csv'))
 fp=r/'raw'/(cell+'_collapsed_counts.csv.gz');xx.to_csv(fp)
 lp=r/'results'/(cell+'_target_logcpm.csv');cp=r/'results'/(cell+'_target_cpm.csv')
 z=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(rp),str(fp),str(lp),str(cp)],capture_output=True,text=True)
 assert z.returncode==0,z.stderr
 logs=pd.read_csv(lp,index_col=0);cpm=pd.read_csv(cp,index_col=0)
 for gene in targets:
  if gene not in logs.index:continue
  temp=mm[['patientid','source','stage','hrs after symptom onset']].copy()
  temp['gene']=gene;temp['cell']=cell;temp['logcpm']=logs.loc[gene,temp.index].to_numpy();temp['cpm']=cpm.loc[gene,temp.index].to_numpy()
  allrows.append(temp.reset_index())
 counts.append(dict(cell=cell,raw_libraries=len(x.columns),collapsed_samples=len(xx.columns),patients=mm.patientid.nunique()))
data=pd.concat(allrows,ignore_index=True);data.to_csv(r/'results/sorted_patient_gene_values.csv',index=False)
def bh(a):
 a=np.asarray(a);bad=~np.isfinite(a);a=np.where(bad,1.,a);o=np.argsort(a);q=np.empty(len(a));q[o]=np.minimum.accumulate((a[o]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];q=np.minimum(q,1);q[bad]=np.nan;return q
rows=[]
def compare(sub,gene,cell,kind,st,a,b,field):
 detected=float(np.mean(sub.cpm>=1)) if len(sub) else 0.
 z=sub.groupby(['patientid',field]).logcpm.mean().unstack(field)
 if a not in z or b not in z:z=pd.DataFrame(columns=[a,b])
 z=z[[a,b]].dropna();v=(z[b]-z[a]).to_numpy()
 pv=float(wilcoxon(v,method='auto').pvalue) if len(v)>=3 and detected>=.2 and np.any(v) else np.nan
 rows.append(dict(gene=gene,cell=cell,test_family=kind,stratum=st,contrast=b+'_vs_'+a,n_pairs=len(v),detectable_fraction=detected,median_change=float(np.median(v)) if len(v) else np.nan,fraction_up=float(np.mean(v>0)) if len(v) else np.nan,p_raw=pv,reason='tested' if np.isfinite(pv) else 'insufficient_pairs_or_detection_or_constant'))
for cell in ['monos','neuts']:
 for gene in targets:
  dd=data[(data.cell==cell)&(data.gene==gene)]
  for source in ['blood','hematoma']:compare(dd[dd.source==source],gene,cell,'temporal_28',source,'early','late','stage')
  for stage in ['early','late']:compare(dd[dd.stage==stage],gene,cell,'tissue_28',stage,'blood','hematoma','source')
out=pd.DataFrame(rows);out['fdr']=out.groupby('test_family').p_raw.transform(lambda a:bh(a))
out.to_csv(r/'results/sorted_validation_tests.csv',index=False)
summary={'cohort_counts':counts,'significant':out[out.fdr<.05].to_dict('records'),'GRAMD1C':out[out.gene=='GRAMD1C'].to_dict('records')}
(r/'results/sorted_summary.json').write_text(json.dumps(summary,indent=2))
(r/'logs/sorted.exit').write_text('0');print(json.dumps(summary,indent=2),flush=True)
