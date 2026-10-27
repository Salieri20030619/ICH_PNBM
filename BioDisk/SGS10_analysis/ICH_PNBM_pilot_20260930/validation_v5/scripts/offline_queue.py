from pathlib import Path
import os,sys,json,time,hashlib,subprocess,fcntl,traceback,py_compile
P=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930/validation_v5')
L=P/'logs/offline_queue';L.mkdir(parents=True,exist_ok=True)
lock=open(L/'queue.lock','w')
try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
except BlockingIOError:raise SystemExit('An offline queue is already running')
PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'
R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
def digest(f):
 h=hashlib.sha256()
 with open(f,'rb') as x:
  for b in iter(lambda:x.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def atomic(path,obj):
 temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf-8');temp.replace(path)
state=dict(status='running',pid=os.getpid(),started=time.strftime('%Y-%m-%d %H:%M:%S'),phase='input_integrity',steps=[])
atomic(L/'status.json',state)
try:
 manifest=json.loads((P/'metadata/download_manifest.json').read_text())
 assert len(manifest)==63
 for row in manifest:
  f=P.parent/'validation_v4/raw/singlecell'/row['name']
  assert f.is_file() and f.stat().st_size==row['size'],row['name']
  assert digest(f)==row['sha256'],row['name']
 assert json.loads((P/'logs/sc_analysis_complete.json').read_text())=={'ok':True,'samples':21}
 atomic(L/'input_integrity.json',dict(ok=True,files=63,samples=21,verified_at=time.strftime('%Y-%m-%d %H:%M:%S')))
 steps=[
 ('patient_level_analysis',[PY,str(P/'scripts/04_offline_extension.py')],['results/extension/all_sample_celltype_metrics.csv','results/extension/program_correlations.csv','results/extension/brain_array_tests.csv']),
 ('alternative_annotation',[R,'--vanilla',str(P/'scripts/05_published_annotation.R')],['results/extension/published_annotation_plxnb3.csv']),
 ('report_and_package',[PY,str(P/'scripts/06_offline_report.py')],['results/extension/PLXNB3_extension_auto_report_zh.txt','results/extension/PLXNB3_extension_overview.png','PLXNB3_extension_offline_results.zip'])]
 for name,command,outputs in steps:
  script=Path(command[-1]);signature=digest(script);marker=L/(name+'.done.json')
  if script.suffix=='.py':py_compile.compile(str(script),doraise=True)
  state['phase']=name;atomic(L/'status.json',state)
  if marker.exists():
   old=json.loads(marker.read_text())
   if old.get('script_sha256')==signature and all((P/f).is_file() and digest(P/f)==old.get('outputs',{}).get(f) for f in outputs):
    state['steps'].append(dict(name=name,status='cached'));continue
  with open(L/(name+'.log'),'a',buffering=1) as log:
   log.write('\nSTART '+time.strftime('%Y-%m-%d %H:%M:%S')+'\n')
   result=subprocess.run(command,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2','MKL_NUM_THREADS':'2'})
  assert result.returncode==0,f'{name}: exit {result.returncode}, see {name}.log'
  assert all((P/f).is_file() and (P/f).stat().st_size>0 for f in outputs)
  atomic(marker,dict(ok=True,script_sha256=signature,exit_code=0,outputs={f:digest(P/f) for f in outputs}))
  state['steps'].append(dict(name=name,status='complete',exit_code=0));atomic(L/'status.json',state)
 state.update(status='complete',phase='finished',exit_code=0,finished=time.strftime('%Y-%m-%d %H:%M:%S'))
 atomic(L/'status.json',state)
 print('OFFLINE_QUEUE_COMPLETE',flush=True)
except Exception as e:
 state.update(status='failed',exit_code=1,error=str(e),finished=time.strftime('%Y-%m-%d %H:%M:%S'))
 atomic(L/'status.json',state);traceback.print_exc();sys.exit(1)
