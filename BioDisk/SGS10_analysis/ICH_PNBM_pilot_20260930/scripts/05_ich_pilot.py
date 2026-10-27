from pathlib import Path
import pandas as pd,numpy as np,json,re,gzip,hashlib
from scipy.stats import wilcoxon
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930')
modules={
'antigen_presentation':['HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','CD74','CIITA'],
'interferon_response':['ISG15','IFIT1','IFIT3','MX1','OAS1','STAT1'],
'neutrophil_associated':['FCGR3B','CSF3R','S100A8','S100A9','CEACAM8','MMP8']}
(r/'metadata/ich_prespecified_panels.json').write_text(json.dumps(modules,indent=2))
def bh(p):
 p=np.asarray(p);o=np.argsort(p);q=np.empty(len(p));q[o]=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];return np.minimum(q,1)
def canon(c):
 pid,t=c.split('_',1);return pid+'_'+{'bas':'base','base':'base','72h':'72h','7da':'7day','7day':'7day'}[t]
mr=pd.read_csv(r/'raw/mrna.csv.gz',index_col=0);mi=pd.read_csv(r/'raw/mirna.csv.gz',index_col=0)
original=mr.columns.tolist();mr.columns=[canon(c) for c in mr.columns];mi.columns=[canon(c) for c in mi.columns]
assert mr.columns.is_unique and mi.columns.is_unique and mr.index.is_unique
assert np.isfinite(mr.to_numpy()).all() and (mr.to_numpy()>=0).all()
assert np.isfinite(mi.to_numpy()).all() and (mi.to_numpy()>=0).all()
meta=pd.DataFrame([dict(sample=c,source_column=src,patient=c.split('_')[0],time=c.split('_')[1]) for c,src in zip(mr.columns,original)])
soft=gzip.open(r/'raw/geo_family.soft.gz','rt').read()
sm=[]
for block in soft.split('^SAMPLE = ')[1:]:
 title=re.search(r'!Sample_title = (.+)',block).group(1).strip()
 typ=title.split('-')[0];pid=re.search(r'(?:mRNA|miRNA)-(\d+)',title).group(1)
 tp=re.search(r'[-\u2013\u2014]\s*(base|72h|7day)\s*$',title).group(1)
 bt=re.search(r'!Sample_characteristics_ch1 = batch: (.+)',block)
 sm.append(dict(accession=block.splitlines()[0],assay=typ,patient=pid,time=tp,batch=bt.group(1).strip() if bt else 'unreported'))
smeta=pd.DataFrame(sm);smeta.to_csv(r/'metadata/geo_sample_metadata.csv',index=False)
meta=meta.merge(smeta[smeta.assay=='mRNA'].drop(columns='assay'),on=['patient','time'],how='left',validate='one_to_one')
assert meta.accession.notna().all()
meta.to_csv(r/'metadata/ich_sample_map.csv',index=False)
pd.crosstab(meta.batch,meta.time).to_csv(r/'results/ich_batch_by_time.csv')
cover=[];log=np.log2(mr+1);base_cols=meta.loc[meta.time=='base','sample']
centered=log.sub(log[base_cols].median(axis=1),axis=0)
scores=pd.DataFrame(index=mr.columns)
gene_rows=[]
contrasts=[('base','72h'),('72h','7day'),('base','7day')]
for k,genes in modules.items():
 present=[g for g in genes if g in mr.index]
 cover.append(dict(module=k,expected=len(genes),present=len(present),missing=';'.join(set(genes)-set(present))))
 assert len(present)>=4
 scores[k]=centered.loc[present].mean(axis=0)
 for g in present:
  for t0,t1 in contrasts:
   pids=sorted(set(meta.loc[meta.time==t0,'patient'])&set(meta.loc[meta.time==t1,'patient']))
   delta=log.loc[g,[p+'_'+t1 for p in pids]].to_numpy()-log.loc[g,[p+'_'+t0 for p in pids]].to_numpy()
   gene_rows.append(dict(module=k,gene=g,contrast=t1+'_vs_'+t0,n_pairs=len(pids),median_log2_change=np.median(delta),fraction_increased=np.mean(delta>0)))
scores['presentation_minus_neutrophil']=scores.antigen_presentation-scores.neutrophil_associated
scores=scores.join(meta.set_index('sample')[['patient','time','batch']]);scores.to_csv(r/'results/ich_patient_module_scores.csv')
pd.DataFrame(cover).to_csv(r/'results/ich_module_coverage.csv',index=False)
pd.DataFrame(gene_rows).to_csv(r/'results/ich_gene_effects.csv',index=False)
rng=np.random.default_rng(20260930)
rows=[]
def test_row(name,contrast,delta,family):
 p=float(wilcoxon(delta,alternative='two-sided',method='auto').pvalue) if np.any(delta!=0) else 1.
 med=np.median(delta)
 bs=np.median(rng.choice(delta,(2000,len(delta))),axis=1)
 return dict(feature=name,contrast=contrast,n_pairs=len(delta),median_change=med,median_change_ci_lo=np.quantile(bs,.025),median_change_ci_hi=np.quantile(bs,.975),
 fraction_increased=np.mean(delta>0),fraction_decreased=np.mean(delta<0),p_raw=p,test_family=family)
for k in list(modules)+['presentation_minus_neutrophil']:
 for t0,t1 in contrasts:
  a=scores[scores.time==t0].set_index('patient')[k];b=scores[scores.time==t1].set_index('patient')[k];pids=a.index.intersection(b.index)
  rows.append(test_row(k,t1+'_vs_'+t0,(b.loc[pids]-a.loc[pids]).to_numpy(),'12_panel_contrasts'))
out=pd.DataFrame(rows);out['fdr_bh']=bh(out.p_raw);out.to_csv(r/'results/ich_module_pilot.csv',index=False)
mi_rows=[]
for g in ['Hsa-Mir-3613_3p*','Hsa-Mir-3613_5p','Hsa-Mir-3690_3p*','Hsa-Mir-3690_5p']:
 assert g in mi.index
 for t0,t1 in contrasts:
  pids=sorted(set(c.split('_')[0] for c in mi if c.endswith('_'+t0))&set(c.split('_')[0] for c in mi if c.endswith('_'+t1)))
  delta=np.log2(mi.loc[g,[p+'_'+t1 for p in pids]].to_numpy()+1)-np.log2(mi.loc[g,[p+'_'+t0 for p in pids]].to_numpy()+1)
  mi_rows.append(test_row(g,t1+'_vs_'+t0,delta,'12_published_miRNA_arm_contrasts'))
mir=pd.DataFrame(mi_rows);mir['fdr_bh']=bh(mir.p_raw);mir.to_csv(r/'results/ich_mirna_reproduction.csv',index=False)
fig,axs=plt.subplots(1,4,figsize=(13,3.8),sharex=True)
for ax,k in zip(axs,list(modules)+['presentation_minus_neutrophil']):
 for pid,part in scores.groupby('patient'):
  v=part.set_index('time')[k].reindex(['base','72h','7day'])
  ax.plot([0,1,2],v,color='#7fa1b4',alpha=.3,lw=.7)
 med=scores.groupby('time')[k].median().reindex(['base','72h','7day'])
 ax.plot([0,1,2],med,color='#a74936',lw=2,marker='o')
 ax.set_title(k.replace('_','\n'),fontsize=10)
 ax.set_xticks([0,1,2],['<24h','72h','7d']);ax.spines[['top','right']].set_visible(False)
 ax.set_ylabel('Baseline-centered mean log2 expression')
fig.suptitle('ICH longitudinal pilot: paired patients, descriptive immune panels')
fig.tight_layout();fig.savefig(r/'results/ich_pilot.png',dpi=160);plt.close(fig)
checks=dict(mrna_rows=mr.shape[0],mrna_samples=mr.shape[1],mirna_rows=mi.shape[0],mirna_samples=mi.shape[1],n_patients=meta.patient.nunique(),
 time_counts=meta.time.value_counts().to_dict(),complete_three_time_patients=int((meta.groupby('patient').time.nunique()==3).sum()),
 assay_sample_sets_match=set(mr.columns)==set(mi.columns),patients_with_multiple_batches=int((meta.groupby('patient').batch.nunique()>1).sum()),module_coverage=cover,
 nominal_panel_tests=int((out.p_raw<.05).sum()),fdr_panel_tests=int((out.fdr_bh<.05).sum()),results=out.to_dict('records'),
 limits=['Processed normalized values; log2(value+1) paired pilot, not raw-count differential expression.',
 'Panels are predeclared descriptive proxies, not validated cell activities or direct cell proportions.',
 'Whole-blood changes can reflect composition; no cell-intrinsic mechanism established.',
 'No healthy controls or verified patient-level outcome/infection labels in the downloaded sample metadata.',
 'Time samples are biologically paired; missing timepoints not imputed.',
 'miRNA arms are reproduction checks of published candidates; not new discoveries.',
 'Novelty of antigen-presentation/neutrophil balance requires separate literature and independent-cohort confirmation.'])
(r/'results/ich_summary.json').write_text(json.dumps(checks,indent=2))
(r/'logs/ich_pilot.exit').write_text('0')
print(json.dumps(checks,indent=2),flush=True)
