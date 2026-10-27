from pathlib import Path
import json,itertools
import pandas as pd,numpy as np
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon,spearmanr
from threadpoolctl import threadpool_limits
threadpool_limits(2)
B=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930')
P=B/'validation_v5';V=B/'validation_v6_external'
for f in ['scripts','metadata','results','logs']:(V/f).mkdir(parents=True,exist_ok=True)
programs=json.loads((P/'metadata/gene_programs.json').read_text())
plan={'purpose':'Independent hematoma cohort state replication using identical matching rules; sparse detection limits eligible patients','matching':'Within patient and OL class; nearest 5 negative cells with replacement in standardized log1p(library minus PLXNB3) and log1p(n_genes minus PLXNB3 detection); maximum distance 0.5; >=5 matched positives and >=20 negative candidates','settings':['standard','strict_QC','lineage_specific'], 'lineage_specific':'SOX10 detectable OR >=2 of MOG MAG MOBP detectable, applied within original OL class','multiplicity':'Six program tests within each sensitivity setting; gene-deletion sensitivity does not select a new primary result; retained all six programs','gene_deletion':'On fixed standard matches, delete each marker in each program and recompute patient effect and the full six-test family','no_network':True}
(V/'metadata/audit_plan.json').write_text(json.dumps(plan,indent=2))
def bh(x):
 x=np.asarray(x,float);bad=~np.isfinite(x);x=np.where(bad,1,x);o=np.argsort(x);q=np.empty(len(x));q[o]=np.minimum.accumulate((x[o]*len(x)/np.arange(1,len(x)+1))[::-1])[::-1];q=np.minimum(q,1);q[bad]=np.nan;return q
effects=[];diag=[];gene_effects=[];unmatched=[]
for gsm in sorted(pd.read_csv(P/'results/extension/all_sample_celltype_metrics.csv').query("series==333172 and source=='hematoma'").gsm.unique()):
 audit=pd.read_csv(P/'results'/f'{gsm}_cell_audit.csv.gz').set_index('barcode')
 counts=pd.read_csv(P/'results'/f'{gsm}_marker_counts.csv.gz',index_col=0)
 audit=audit.loc[counts.index]
 lineage=(counts.SOX10>=1)|((counts[['MOG','MAG','MOBP']]>=1).sum(axis=1)>=2)
 for setting in ['standard','strict_QC','lineage_specific']:
  keep=audit.qc&(audit.label=='Oligodendrocyte')
  if setting=='strict_QC':keep=keep&audit.strict_qc
  if setting=='lineage_specific':keep=keep&lineage
  a=audit[keep];c=counts.loc[a.index];y=a.PLXNB3_UMI.to_numpy();pos=np.flatnonzero(y>=1);neg=np.flatnonzero(y==0)
  if len(pos)<5 or len(neg)<20:
   diag.append(dict(gsm=gsm,setting=setting,evaluable=False,n_positive=len(pos),n_negative=len(neg)));continue
  cov=np.column_stack([np.log1p(a.library.to_numpy()-y),np.log1p(a.n_genes.to_numpy()-(y>=1))])
  scale=cov.std(axis=0,ddof=1);scale=np.where(scale==0,1,scale);z=(cov-cov.mean(axis=0))/scale
  dist,ix=cKDTree(z[neg]).query(z[pos],k=5)
  ok=dist[:,-1]<=.5;pp=pos[ok];nn=neg[ix[ok]]
  if len(pp)<5:
   diag.append(dict(gsm=gsm,setting=setting,evaluable=False,n_positive=len(pos),n_negative=len(neg),n_matched=len(pp)));continue
  row=dict(gsm=gsm,setting=setting,evaluable=True,n_positive=len(pos),n_negative=len(neg),n_matched=len(pp),match_fraction=len(pp)/len(pos),n_unique_controls=len(np.unique(nn)),max_distance=float(dist[ok].max()))
  for j,name in enumerate(['library','genes']):
   row[name+'_SMD_before']=float((cov[pos,j].mean()-cov[neg,j].mean())/scale[j])
   row[name+'_SMD_after']=float((cov[pp,j].mean()-cov[nn,j].mean())/scale[j])
  diag.append(row)
  for name,genes in programs.items():
   scores=np.log1p(c[genes].to_numpy()*10000/a.library.to_numpy()[:,None])
   per_gene=(scores[pp]-scores[nn].mean(axis=1)).mean(axis=0)
   effects.append(dict(gsm=gsm,setting=setting,program=name,delta=float(per_gene.mean()),n_matched=len(pp)))
   if setting=='standard':
    for g,delta in zip(genes,per_gene):gene_effects.append(dict(gsm=gsm,program=name,gene=g,delta=float(delta)))
 print('MATCHED',gsm,flush=True)
ef=pd.DataFrame(effects);dg=pd.DataFrame(diag);ge=pd.DataFrame(gene_effects)
ef.to_csv(V/'results/matched_patient_effects.csv',index=False);dg.to_csv(V/'results/matching_diagnostics.csv',index=False);ge.to_csv(V/'results/matched_gene_effects.csv',index=False)
def stats(v):
 v=np.asarray(v,float);v=v[np.isfinite(v)]
 return dict(n=len(v),median_delta=float(np.median(v)) if len(v) else None,mean_delta=float(np.mean(v)) if len(v) else None,n_up=int((v>0).sum()),n_down=int((v<0).sum()),p_raw=float(wilcoxon(v).pvalue) if len(v)>=3 and np.any(v) else None)
rows=[]
for setting in ['standard','strict_QC','lineage_specific']:
 for name in programs:rows.append(dict(setting=setting,program=name,**stats(ef.loc[(ef.setting==setting)&(ef.program==name),'delta'])))
out=pd.DataFrame(rows);out['fdr_6']=out.groupby('setting').p_raw.transform(lambda q:bh(q));out.to_csv(V/'results/matched_state_tests.csv',index=False)
# Leave-one-patient-out, preserve six-program family.
loo=[]
for gsm in ef.gsm.unique():
 r=[]
 for name in programs:r.append(dict(excluded=gsm,program=name,**stats(ef.loc[(ef.setting=='standard')&(ef.program==name)&(ef.gsm!=gsm),'delta'])))
 r=pd.DataFrame(r);r['fdr_6']=bh(r.p_raw);loo.append(r)
pd.concat(loo).to_csv(V/'results/leave_one_patient_out.csv',index=False)
# Delete each marker from its program; recompute family with remaining five unchanged.
delete=[]
for name,genes in programs.items():
 for gene in genes:
  rows=[]
  for other in programs:
   vals=ge[(ge.program==other)&((ge.gene!=gene) if other==name else True)].groupby('gsm').delta.mean()
   rows.append(dict(program=other,**stats(vals)))
  rr=pd.DataFrame(rows);rr['fdr_6']=bh(rr.p_raw);t=rr[rr.program==name].iloc[0].to_dict();t.update(deleted_gene=gene);delete.append(t)
pd.DataFrame(delete).to_csv(V/'results/leave_one_marker_out.csv',index=False)

print(out.to_json(orient='records'),flush=True)
(V/'logs/analysis_complete.json').write_text(json.dumps({'ok':True,'eligible_patients':len(ef.gsm.unique())}))
