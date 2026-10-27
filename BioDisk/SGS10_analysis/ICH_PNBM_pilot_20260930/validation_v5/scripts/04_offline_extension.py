from pathlib import Path
import json, itertools, math, gzip, io
import numpy as np
import pandas as pd
from scipy.stats import rankdata, wilcoxon
from threadpoolctl import threadpool_limits
threadpool_limits(2)
P = Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930/validation_v5')
O = P/'results/extension'; O.mkdir(exist_ok=True)
programs = json.loads((P/'metadata/gene_programs.json').read_text())
summary = pd.read_csv(P/'results/singlecell_sample_summary.csv')
assert summary.gsm.nunique() == 21
frames = [pd.read_csv(P/'results'/f'{g}_programs.csv') for g in summary.gsm.unique()]
pr = pd.concat(frames, ignore_index=True)
d = summary.merge(pr, on=['gsm','title','series','qc','celltype','n_cells','total_UMI'], how='left', validate='one_to_one')
d['hours'] = pd.to_numeric(d.title.str.extract(r'ICH_(\d+)h', expand=False), errors='coerce')
d['patient'] = pd.to_numeric(d.title.str.extract(r'Patient (\d+)', expand=False), errors='coerce')
d['source'] = np.where(d.title.str.contains('Hematoma'), 'hematoma', np.where(d.title.str.contains('Blood'), 'blood', 'perihematoma'))
d.to_csv(O/'all_sample_celltype_metrics.csv', index=False)

def bh(vals):
    vals=np.asarray(vals,float); bad=~np.isfinite(vals); a=np.where(bad,1,vals); ix=np.argsort(a)
    out=np.empty(len(a)); out[ix]=np.minimum.accumulate((a[ix]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1]
    out=np.minimum(out,1); out[bad]=np.nan; return out

perms={}
def spearman_exact(x,y):
    a=pd.DataFrame({'x':x,'y':y}).dropna(); n=len(a)
    if n<3 or a.x.nunique()<2 or a.y.nunique()<2: return n,np.nan,np.nan
    assert n<=9, 'Exact enumeration bound exceeded'
    x=rankdata(a.x); y=rankdata(a.y); x-=x.mean(); y-=y.mean()
    denom=np.linalg.norm(x)*np.linalg.norm(y); obs=x@y
    if n not in perms:
        perms[n]=np.fromiter(itertools.chain.from_iterable(itertools.permutations(range(n))),dtype=np.int16).reshape(-1,n)
    null=x[perms[n]]@y
    return n,float(obs/denom),float(np.mean(np.abs(null)>=abs(obs)-1e-10))

def paired_result(x):
    a=np.asarray(x,float); a=a[np.isfinite(a)]; n=len(a)
    return dict(n=n,median_delta=float(np.median(a)) if n else None,mean_delta=float(a.mean()) if n else None,
                n_up=int((a>0).sum()),n_down=int((a<0).sum()),n_zero=int((a==0).sum()),
                p_raw=float(wilcoxon(a).pvalue) if n>=3 and np.any(a) else None)

source_tests=[]; time_tests=[]; corr=[]; matched=[]
for qc in ['standard','strict']:
    z=d[(d.series==266873)&(d.qc==qc)].copy()
    ol=z[(z.celltype=='Oligodendrocyte')&(z.n_cells>=30)].set_index('gsm')
    for ct in ['Myeloid','Astrocyte']:
        ref=z[(z.celltype==ct)&(z.n_cells>=30)].set_index('gsm')
        ids=ol.index.intersection(ref.index)
        delta=np.log2(ol.loc[ids,'pseudobulk_CPM']+1)-np.log2(ref.loc[ids,'pseudobulk_CPM']+1)
        source_tests.append(dict(qc=qc,contrast='Oligodendrocyte_vs_'+ct,**paired_result(delta)))
    n,r,p=spearman_exact(ol.hours,ol.pseudobulk_CPM)
    time_tests.append(dict(qc=qc,endpoint='OL_PLXNB3_CPM_vs_onset_hours',n=n,rho=r,p_raw=p,fdr=None,family='single_prespecified_temporal_endpoint'))
    for name in programs:
        n,r,p=spearman_exact(ol.pseudobulk_CPM,ol[name])
        corr.append(dict(qc=qc,program=name,n=n,rho=r,p_raw=p))
    ol.to_csv(O/f'perihematoma_ol_patient_values_{qc}.csv')
src=pd.DataFrame(source_tests); src['fdr_2']=src.groupby('qc').p_raw.transform(lambda x:bh(x));src.to_csv(O/'source_localization_tests.csv',index=False)
pd.DataFrame(time_tests).to_csv(O/'temporal_tests.csv',index=False)
co=pd.DataFrame(corr);co['fdr_6']=co.groupby('qc').p_raw.transform(lambda x:bh(x));co.to_csv(O/'program_correlations.csv',index=False)

# Within-patient, UMI-stratified state contrasts. Biological inference uses patient effects only.
for gsm in d[d.series==266873].gsm.unique():
    cells=pd.read_csv(P/'results'/f'{gsm}_cell_audit.csv.gz')
    a=cells[cells.qc & (cells.label=='Oligodendrocyte')].copy()
    a['positive']=a.PLXNB3_UMI>=1
    if a.positive.sum()<5 or (~a.positive).sum()<20: continue
    a['depth_bin']=pd.qcut(a.library.rank(method='first'),5,labels=False)
    for name in programs:
        effects=[]; weights=[]
        for _,b in a.groupby('depth_bin'):
            pos=b[b.positive]; neg=b[~b.positive]
            if len(pos)==0 or len(neg)<10:continue
            effects.append(pos[name+'_score'].mean()-neg[name+'_score'].mean());weights.append(len(pos))
        if weights:matched.append(dict(gsm=gsm,program=name,n_positive=int(a.positive.sum()),n_negative=int((~a.positive).sum()),n_positive_used=sum(weights),delta=float(np.average(effects,weights=weights))))
mt=pd.DataFrame(matched);mt.to_csv(O/'within_ol_depth_stratified_patient_effects.csv',index=False)
rr=[]
for name in programs:
    x=mt.loc[mt.program==name,'delta'] if len(mt) else []
    rr.append(dict(program=name,**paired_result(x)))
rr=pd.DataFrame(rr);rr['fdr_6']=bh(rr.p_raw);rr.to_csv(O/'within_ol_state_tests.csv',index=False)

# Six label-consistent pairs: no inferential test when gene detectability is inadequate.
pairs=[]; tests=[]
z=d[(d.series==333172)&(d.qc=='standard')]
for ct in ['All','Myeloid','Myeloid_clean','Oligodendrocyte']:
    a=z[(z.celltype==ct)&(z.n_cells>=30)].pivot(index='patient',columns='source',values='pseudobulk_CPM').dropna()
    if not {'blood','hematoma'}.issubset(a.columns):continue
    a['delta']=np.log2(a.hematoma+1)-np.log2(a.blood+1)
    for pid,row in a.iterrows():pairs.append(dict(patient=int(pid),celltype=ct,**row.to_dict()))
    result=paired_result(a.delta); detected=int(((a.blood>=1)|(a.hematoma>=1)).sum())
    evaluable=len(a)>=3 and detected>=3
    if not evaluable:result['p_raw']=None
    tests.append(dict(celltype=ct,detectable_patients=detected,evaluable=evaluable,**result))
pd.DataFrame(pairs).to_csv(O/'independent_hematoma_blood_pairs.csv',index=False)
tests=pd.DataFrame(tests);tests['fdr_4']=bh(list(tests.p_raw)+[np.nan]*(4-len(tests)))[:len(tests)];tests.to_csv(O/'independent_pair_tests.csv',index=False)

# GSE24265: preserve contralateral grey and white matter separately, no pooled-control inference.
meta=pd.read_csv(P/'metadata/gse24265_samples.csv')
values=pd.read_csv(P/'results/gse24265_plxnb3_probe_values.csv',index_col=0)
assert list(values.index)==['205957_at']
a=meta.assign(log2_RMA=lambda q:q.gsm.map(values.loc['205957_at'])).pivot(index='patient',columns='source',values='log2_RMA')
a.to_csv(O/'brain_array_patient_values.csv')
rows=[]
for ref in ['CW','CG']:
    q=a[['PH',ref]].dropna();diff=q.PH-q[ref]
    rows.append(dict(contrast='PH_vs_'+ref,probe='205957_at',**paired_result(diff)))
rows=pd.DataFrame(rows);rows['fdr_2']=bh(rows.p_raw);rows.to_csv(O/'brain_array_tests.csv',index=False)

# Positive-cell evidence: retain rather than discard mixed lineage RNA.
pos=[]
for gsm in d.gsm.unique():
    c=pd.read_csv(P/'results'/f'{gsm}_cell_audit.csv.gz');c=c[c.qc&(c.PLXNB3_UMI>0)].copy();c.insert(0,'gsm',gsm);pos.append(c)
pd.concat(pos,ignore_index=True).to_csv(O/'all_positive_cell_audit.csv.gz',index=False)
print('Offline core analysis finished',flush=True)
