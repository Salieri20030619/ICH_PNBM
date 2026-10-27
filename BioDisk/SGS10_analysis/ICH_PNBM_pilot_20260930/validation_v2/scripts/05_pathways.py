from pathlib import Path
import pandas as pd,numpy as np,json,re,gzip
from scipy.stats import wilcoxon,ttest_1samp
b=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930');r=b/'validation_v2'
def bh(v):
 v=np.asarray(v);bad=~np.isfinite(v);v=np.where(bad,1.,v);o=np.argsort(v);q=np.empty(len(v));q[o]=np.minimum.accumulate((v[o]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1];q=np.minimum(q,1);q[bad]=np.nan;return q
x=pd.read_csv(b/'raw/mrna.csv.gz',index_col=0);x.columns=[re.sub('_bas$','_base',re.sub('_7da$','_7day',c)) for c in x]
d=pd.read_csv(r/'raw/external.xls.gz',sep='\t');dup=d.symbol.duplicated(False)|d.symbol.isna()
ex=d.loc[~dup].set_index('symbol')[[f'ICH{i}' for i in range(1,23)]]
# Verify metadata rather than relying solely on odd/even naming.
soft=gzip.open(r/'raw/external.soft.gz','rt').read();mapping=[]
for block in soft.split('^SAMPLE = ')[1:]:
 title=re.search(r'!Sample_title = (.+)',block).group(1)
 mt=re.fullmatch(r'Patient(\d+)-(first|second)',title.strip())
 col=re.search(r'!Sample_description = (ICH\d+)',block).group(1)
 pid=int(mt.group(1));time=mt.group(2)
 assert col=='ICH'+str(pid*2-(time=='first'))
 mapping.append(dict(patient=pid,time=time,column=col))
pd.DataFrame(mapping).to_csv(r/'metadata/external_pair_map.csv',index=False)
primary=set((r/'metadata/ich_gene_universe.txt').read_text().splitlines());external=set((r/'metadata/external_gene_universe.txt').read_text().splitlines())
common=primary&external
log=np.log2(x+1);elog=np.log2(ex+1);defs=[]
for line in (r/'raw/ReactomePathways.gmt').read_text().splitlines():
 parts=line.split('\t');genes=sorted(set(parts[2:])&common)
 if 15<=len(genes)<=500:defs.append((parts[0],parts[1],genes))
(r/'metadata/pathway_definitions.json').write_text(json.dumps([dict(name=n,id=i,genes=g) for n,i,g in defs],indent=2))
score=np.array([log.loc[g].mean(axis=0) for n,i,g in defs]);escore=np.array([elog.loc[g].mean(axis=0) for n,i,g in defs])
pd.DataFrame(score,index=[i for n,i,g in defs],columns=log.columns).to_csv(r/'results/pathway_scores_discovery.csv.gz')
pd.DataFrame(escore,index=[i for n,i,g in defs],columns=elog.columns).to_csv(r/'results/pathway_scores_external.csv.gz')
rows=[];delta_store={}
for t0,t1 in [('base','72h'),('72h','7day'),('base','7day')]:
 ids=sorted(set(c.split('_')[0] for c in x if c.endswith('_'+t0))&set(c.split('_')[0] for c in x if c.endswith('_'+t1)))
 i0=[x.columns.get_loc(p+'_'+t0) for p in ids];i1=[x.columns.get_loc(p+'_'+t1) for p in ids]
 delta=score[:,i1]-score[:,i0];contrast=t1+'_vs_'+t0;delta_store[contrast]=delta
 for j,(n,i,g) in enumerate(defs):
  v=delta[j];pw=float(wilcoxon(v,method='auto').pvalue) if np.any(v) else np.nan
  rows.append(dict(pathway=n,id=i,contrast=contrast,n_pairs=len(ids),genes=len(g),mean_change=np.mean(v),median_change=np.median(v),fraction_up=np.mean(v>0),p_rank=pw,p_t=float(ttest_1samp(v,0).pvalue)))
out=pd.DataFrame(rows);out['fdr_rank_joint']=bh(out.p_rank);out['fdr_t_joint_sensitivity']=bh(out.p_t)
out.to_csv(r/'results/pathway_discovery_all.csv',index=False)
deltae=escore[:,1::2]-escore[:,::2];er=[]
for j,(n,i,g) in enumerate(defs):
 v=deltae[j];er.append(dict(pathway=n,id=i,n_pairs=11,mean_change=np.mean(v),median_change=np.median(v),fraction_up=np.mean(v>0),p_rank=float(wilcoxon(v,method='auto').pvalue) if np.any(v) else np.nan,p_t=float(ttest_1samp(v,0).pvalue)))
ee=pd.DataFrame(er);ee['fdr_rank_all']=bh(ee.p_rank);ee['fdr_t_all_sensitivity']=bh(ee.p_t);ee.to_csv(r/'results/pathway_external_all.csv',index=False)
merged=out.merge(ee.drop(columns='pathway'),on='id',suffixes=('_discovery','_external'))
merged['same_direction']=np.sign(merged.mean_change_discovery)==np.sign(merged.mean_change_external)
merged.to_csv(r/'results/pathway_replication.csv',index=False)
g=pd.read_csv(r/'results/ich_all_gene_tests.csv.gz');e=pd.read_csv(r/'results/external_gene_tests.csv.gz')
ge=g.merge(e,on='gene',suffixes=('_discovery','_external'))
ge['same_direction']=np.sign(ge.mean_log2_change_discovery)==np.sign(ge.mean_log2_change_external)
ge.to_csv(r/'results/gene_replication_all.csv.gz',index=False)
sel=ge[(ge.fdr_limma_joint<.05)&ge.same_direction&(ge.p_limma_external<.05)].sort_values('fdr_limma_joint')
sel.to_csv(r/'results/gene_discovery_fdr_external_nominal.csv',index=False)
sm={'pathways':len(defs),'common_measured_genes':len(common),
'pathways_discovery_rank_fdr':int((out.fdr_rank_joint<.05).sum()),'pathways_external_rank_fdr':int((ee.fdr_rank_all<.05).sum()),
'pathways_external_t_fdr_sensitivity':int((ee.fdr_t_all_sensitivity<.05).sum()),
'top_pathways':merged[merged.fdr_rank_joint<.05].sort_values('fdr_rank_joint').head(20).to_dict('records'),
'discovery_fdr_external_nominal_genes':sel.to_dict('records'),
'cautions':['Different sampling window: external second draw ~72h after first, versus discovery 72h after onset.',
'External older FPKM vs newer normalized counts: compare direction and evidence, not pooled absolute expression.',
'Same research group; different BioProjects and earlier publication support a separate cohort but patient-level linkage is unavailable.',
'All 3 discovery contrasts corrected jointly within each modality; external genes/pathways corrected over full universes.',
'Rank test primary for pathway scores; t-test results are sensitivity only; no choice of test by smallest FDR.']}
(r/'results/ich_v2_summary.json').write_text(json.dumps(sm,indent=2))
print(json.dumps(sm,indent=2),flush=True)
(r/'logs/pathways.exit').write_text('0')
