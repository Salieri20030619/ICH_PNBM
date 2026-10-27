import pathlib,json,subprocess,concurrent.futures,re,gzip,hashlib,time
p=pathlib.Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930/validation_v4');raw=p/'raw'/'singlecell'
subprocess.run(['curl','-L','--fail','--max-time','40','-sS','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE266nnn/GSE266873/suppl/filelist.txt','-o',str(p/'raw/filelist266.txt')],check=True)
sizes={}
for f in [p/'raw/filelist266.txt',p/'raw/filelist333.txt']:
 for l in f.read_text().splitlines():
  a=l.split()
  if a and a[0]=='File':sizes[a[1]]=int(a[-2])
files=[]
for s in json.loads((p/'metadata/sc_samples.json').read_text()):
 if s['ambiguous']:continue
 for url in s['urls']:
  name=url.split('/')[-1]
  files.append(dict(name=name,url=url.replace('ftp://','https://'),size=sizes[name],series=s['series']))
files=[f for f in files if f['name'].split('_')[0] in ['GSM9758327','GSM9758337','GSM9758339','GSM9758345','GSM8255347']]
files.sort(key=lambda x:x['size'])
partsdir=raw/'range_parts';partsdir.mkdir(exist_ok=True)
tasks=[];remaining={};records=[]
chunk=2*1024*1024
for f in files:
 out=raw/f['name']
 if out.exists() and out.stat().st_size==f['size']:
  records.append(dict(**f,status='OK',sha256=hashlib.sha256(out.read_bytes()).hexdigest()));continue
 f['ranges']=[(i,min(i+chunk-1,f['size']-1)) for i in range(0,f['size'],chunk)]
 remaining[f['name']]=len(f['ranges'])
 tasks.extend((f,a,b) for a,b in f['ranges'])
def get(t):
 f,a,b=t;dst=partsdir/(f['name']+'.'+str(a));hdr=pathlib.Path(str(dst)+'.hdr');expected=b-a+1
 if dst.exists() and dst.stat().st_size==expected and hdr.exists() and ('content-range: bytes '+str(a)+'-'+str(b)+'/') in hdr.read_text().lower():return f,a,b,True
 for attempt in range(4):
  time.sleep(1+attempt*5)
  rr=subprocess.run(['curl','-L','--fail','--connect-timeout','15','--max-time','240','-sS','-r',str(a)+'-'+str(b),'-D',str(hdr),f['url'],'-o',str(dst)],capture_output=True)
  if rr.returncode==0 and dst.exists() and dst.stat().st_size==expected and ('content-range: bytes '+str(a)+'-'+str(b)+'/') in hdr.read_text().lower():return f,a,b,True
 return f,a,b,False
print('START',len(files),'files',len(tasks),'chunks',sum(f['size'] for f in files),flush=True)
fails=[]
with concurrent.futures.ThreadPoolExecutor(6) as ex:
 for future in concurrent.futures.as_completed([ex.submit(get,t) for t in tasks]):
  f,a,b,ok=future.result()
  if not ok:fails.append([f['name'],a,b]);continue
  remaining[f['name']]-=1
  if remaining[f['name']]==0:
   out=raw/f['name']
   with open(str(out)+'.assembled','wb') as dst:
    for lo,hi in f['ranges']:
     with open(partsdir/(f['name']+'.'+str(lo)),'rb') as src:
      while True:
       buf=src.read(1024*1024)
       if not buf:break
       dst.write(buf)
   with gzip.open(str(out)+'.assembled','rb') as z:
    while z.read(1024*1024):pass
   pathlib.Path(str(out)+'.assembled').rename(out)
   rec=dict(name=f['name'],size=f['size'],status='OK',sha256=hashlib.sha256(out.read_bytes()).hexdigest())
   records.append(rec);print('COMPLETE',f['name'],f['size'],flush=True)
   (p/'metadata/download_manifest.json').write_text(json.dumps(records,indent=2))
(p/'logs/download_complete.json').write_text(json.dumps(dict(ok=not fails,failures=fails,n=len(records))))
print('FINISHED',len(records),'failures',len(fails),flush=True)
