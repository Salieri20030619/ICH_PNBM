import pathlib,gzip,re,subprocess,concurrent.futures,json,time,hashlib
p=pathlib.Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930/validation_v4')
tasks=[]; samples=[]
for series in ['266873','333172']:
 s=gzip.open(p/'raw'/('gse'+series+'.soft.gz'),'rt').read()
 for block in s.split('^SAMPLE = ')[1:]:
  gsm=block.splitlines()[0]; title=re.search(r'!Sample_title = (.+)',block)[1]
  urls=re.findall(r'!Sample_supplementary_file_\d+ = (.+)',block)
  ambiguous=series=='333172' and int(re.search(r'Patient (\d+)',title)[1]) in [2,3,4,5]
  samples.append(dict(series=series,gsm=gsm,title=title,ambiguous=ambiguous,urls=urls))
  if ambiguous:continue
  for u in urls: tasks.append((gsm,u.replace('ftp://','https://')))
(p/'metadata'/'sc_samples.json').write_text(json.dumps(samples,indent=2))
def get(t):
 gsm,u=t; out=p/'raw'/'singlecell'/u.split('/')[-1];out.parent.mkdir(exist_ok=True)
 start=time.time()
 if not out.exists():
  r=subprocess.run(['curl','-L','--fail','--retry','2','--connect-timeout','15','--max-time','1200','-sS',u,'-o',str(out)+'.part'],capture_output=True)
  if r.returncode: return dict(file=out.name,status='FAILED',error=r.stderr.decode(),seconds=time.time()-start)
  pathlib.Path(str(out)+'.part').rename(out)
 try:
  with gzip.open(out,'rb') as f:
   while f.read(1024*1024):pass
  result=dict(file=out.name,status='OK',bytes=out.stat().st_size,sha256=hashlib.sha256(out.read_bytes()).hexdigest(),seconds=round(time.time()-start,2))
 except Exception as e:result=dict(file=out.name,status='INVALID',error=str(e))
 print(json.dumps(result),flush=True)
 return result
with concurrent.futures.ThreadPoolExecutor(8) as ex:
 results=list(ex.map(get,tasks))
(p/'metadata'/'download_manifest.json').write_text(json.dumps(results,indent=2))
(p/'logs'/'download_complete.json').write_text(json.dumps({'ok':all(x['status']=='OK' for x in results),'files':len(results)}))
