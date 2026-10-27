from pathlib import Path
import json,gzip,re,os,ctypes,select,time
import numpy as np,pandas as pd
from scipy.io import mmread
from threadpoolctl import threadpool_limits
threadpool_limits(4)
p=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930/validation_v4');raw=p/'raw/singlecell'
samples=[s for s in json.loads((p/'metadata/sc_samples.json').read_text()) if s['gsm'] in ['GSM9758327','GSM9758337','GSM9758339','GSM9758345','GSM8255347']]
panels={'Myeloid':['LST1','TYROBP','FCER1G','AIF1','CD68','CSF1R','LYZ'],'Astrocyte':['GFAP','AQP4','SLC1A3','ALDH1L1','GJA1','HEPACAM'],'Oligodendrocyte':['PLP1','MBP','MOG','MOBP','MAG','SOX10'],'OPC':['PDGFRA','CSPG4','OLIG1','OLIG2','GPR17'],'Platelet':['PF4','PPBP','GP9','GP1BA','TUBB1'],'Neuron':['SNAP25','SYT1','RBFOX3','STMN2'],'Endothelial':['CLDN5','PECAM1','VWF','KDR','EMCN'],'T_NK':['CD3D','CD3E','TRAC','NKG7','GNLY'],'B':['MS4A1','CD79A','CD79B','CD19'],'Neutrophil':['CSF3R','FCGR3B','CEACAM8','CXCR2']}
genes=sorted(set(sum(panels.values(),[])+['PLXNB3','SEMA5A','CD14','FCGR3A','P2RY12','TMEM119','SALL1','GPR34','C1QA','ARHGDIA','RAC1','MET']))
plan={'target':'PLXNB3','primary':'GSE333172: four hematoma specimens; source localization only; no inferential testing','qc':'200<=detected genes<=6000; mitochondrial fraction<=20%; UMI<=40000','strict_qc':'500<=genes<=6000; mito<=10%; UMI<=40000','annotation':'Provisional marker-based broad class, no PLXNB3 used. >=2 detected markers, mean log1p(10000*UMI/library)>=0.5. Top two scores within 0.15 flagged mixed. Myeloid has >=3 markers and no glial or platelet panel with >=2 detected markers. This does not distinguish monocytes from microglia.','secondary':'GSE266873: one perihematomal 28h specimen; source localization only, no inferential testing','selection':'Download-availability convenience pilot, chosen by file size before cell-expression results; not representative of full cohorts. 503 failures prevent planned 21-sample analysis in this round.','panels':panels}
(p/'metadata/sc_analysis_plan.json').write_text(json.dumps(plan,indent=2))
def run(s):
 gsm=s['gsm'];files={u.split('/')[-1].split('_')[-1]:raw/u.split('/')[-1] for u in s['urls']}
 feats=pd.read_csv(files['features.tsv.gz'],sep='\t',header=None);symbols=feats[1].astype(str).to_numpy()
 with gzip.open(files['matrix.mtx.gz'],'rb') as fh:x=mmread(fh).tocsr()
 assert x.shape[0]==len(symbols)
 bc=pd.read_csv(files['barcodes.tsv.gz'],sep='\t',header=None)[0].astype(str).to_numpy();assert len(bc)==x.shape[1]
 lib=np.asarray(x.sum(axis=0)).ravel();ng=np.asarray((x>0).sum(axis=0)).ravel()
 mt=np.asarray(x[np.char.startswith(symbols.astype(str),'MT-')].sum(axis=0)).ravel()/np.maximum(lib,1)
 qc=(ng>=200)&(ng<=6000)&(mt<=.2)&(lib<=40000);strict=qc&(ng>=500)&(mt<=.1)
 counts={}
 for g in genes:
  inds=np.where(symbols==g)[0];counts[g]=np.asarray(x[inds].sum(axis=0)).ravel() if len(inds) else np.zeros(x.shape[1])
 counts=pd.DataFrame(counts,index=bc)
 score={};detect={}
 for label,gg in panels.items():
  present=[g for g in gg if g in symbols]
  score[label]=np.log1p(counts[present].to_numpy()*10000/np.maximum(lib[:,None],1)).mean(axis=1) if present else np.zeros(len(lib))
  detect[label]=(counts[present]>0).sum(axis=1).to_numpy()
 sc=pd.DataFrame(score,index=bc);dd=pd.DataFrame(detect,index=bc)
 ranked=sc.to_numpy().argsort(axis=1);top=ranked[:,-1];second=ranked[:,-2];names=np.array(list(panels))
 labels=names[top].copy().astype(object);val=sc.to_numpy()[np.arange(len(lib)),top];valid=(val>=.5)&(dd.to_numpy()[np.arange(len(lib)),top]>=2)
 labels[~valid]='Unclassified'
 mixed=valid&(val-sc.to_numpy()[np.arange(len(lib)),second]<.15)&(dd.to_numpy()[np.arange(len(lib)),second]>=2)
 labels[mixed]='Mixed'
 glial=(dd[['Astrocyte','Oligodendrocyte','OPC']]>=2).any(axis=1).to_numpy()
 myeloid=(labels=='Myeloid')&(dd.Myeloid.to_numpy()>=3)&(~glial)&(dd.Platelet.to_numpy()<2)
 info=pd.DataFrame(dict(barcode=bc,library=lib,n_genes=ng,mito_fraction=mt,qc=qc,strict_qc=strict,label=labels,myeloid_clean=myeloid,glial_markers=glial,platelet_markers=dd.Platelet.to_numpy()>=2,PLXNB3_UMI=counts.PLXNB3.to_numpy()))
 for label in panels:info[label+'_score']=sc[label].to_numpy();info[label+'_detected']=dd[label].to_numpy()
 for g in ['SEMA5A','CD14','P2RY12','TMEM119','C1QA']:info[g+'_UMI']=counts[g].to_numpy()
 info.to_csv(p/'results'/(gsm+'_cell_audit.csv.gz'),index=False)
 counts[qc].to_csv(p/'results'/(gsm+'_marker_counts.csv.gz'))
 rows=[]
 for setting,mask in [('standard',qc),('strict',strict)]:
  for label in ['All','Myeloid_clean']+list(panels)+['Mixed','Unclassified']:
   ix=mask if label=='All' else mask&myeloid if label=='Myeloid_clean' else mask&(labels==label)
   n=int(ix.sum());tot=float(lib[ix].sum());y=counts.PLXNB3.to_numpy()[ix]
   rows.append(dict(series=s['series'],gsm=gsm,title=s['title'],qc=setting,celltype=label,n_cells=n,PLXNB3_total_UMI=int(y.sum()),PLXNB3_positive_cells=int((y>=1).sum()),PLXNB3_ge2_cells=int((y>=2).sum()),positive_fraction=float((y>=1).mean()) if n else None,pseudobulk_CPM=float(y.sum()/tot*1e6) if tot else None,total_UMI=tot))
 pd.DataFrame(rows).to_csv(p/'results'/(gsm+'_summary.csv'),index=False)
 print('ANALYZED',gsm,'rawcells',len(lib),'qc',int(qc.sum()),'PLXNB3+',int((counts.PLXNB3.to_numpy()[qc]>0).sum()),'myeloid+',int((counts.PLXNB3.to_numpy()[qc&myeloid]>0).sum()),flush=True)
# Native inotify waits for completed downloads rather than repeatedly inspecting processes.
libc=ctypes.CDLL('libc.so.6',use_errno=True);fd=libc.inotify_init1(0);assert fd>=0
assert libc.inotify_add_watch(fd,str(raw).encode(),0x80|0x8)>=0
done=set()
while len(done)<len(samples):
 for s in samples:
  if s['gsm'] in done:continue
  if (p/'results'/(s['gsm']+'_summary.csv')).exists():done.add(s['gsm']);continue
  if all((raw/u.split('/')[-1]).exists() for u in s['urls']):
   run(s);done.add(s['gsm'])
 if len(done)==len(samples):break
 os.read(fd,65536)
pd.concat([pd.read_csv(p/'results'/(s['gsm']+'_summary.csv')) for s in samples]).to_csv(p/'results/singlecell_sample_summary.csv',index=False)
(p/'logs/sc_analysis_complete.json').write_text(json.dumps(dict(ok=True,samples=len(done))))
print('ALL_ANALYZED',len(done),flush=True)
