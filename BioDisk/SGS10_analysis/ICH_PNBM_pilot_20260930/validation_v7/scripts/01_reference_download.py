from pathlib import Path
import json,gzip,re,subprocess,concurrent.futures,hashlib,time
P=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930/validation_v7')
raw=P/'raw/reference_sc';raw.mkdir(exist_ok=True);parts=raw/'parts';parts.mkdir(exist_ok=True)
sizes={}
for l in (P/'raw/gse157827_files.txt').read_text().splitlines():
 a=l.split()
 if a and a[0]=='File':sizes[a[1]]=int(a[-2])
s=gzip.open(P/'raw/gse157827.soft.gz','rt').read();samples=[];files=[]
for block in s.split('^SAMPLE = ')[1:]:
 gsm=block.splitlines()[0];title=re.search(r'!Sample_title = (.+)',block)[1]
 diagnosis=re.search(r'!Sample_characteristics_ch1 = diagnosis: (.+)',block)[1]
 urls=re.findall(r'!Sample_supplementary_file_\d+ = (.+)',block)
 samples.append(dict(gsm=gsm,title=title,diagnosis=diagnosis,series='157827',ambiguous=False,urls=urls))
 for url in urls:
  name=url.split('/')[-1];files.append(dict(name=name,url=url.replace('ftp://','https://'),size=sizes[name]))
(P/'metadata/reference_samples.json').write_text(json.dumps(samples,indent=2))
tasks=[];remaining={};records=[];chunk=2*1024*1024
for f in sorted(files,key=lambda x:x['size']):
 out=raw/f['name']
 if out.exists() and out.stat().st_size==f['size']:records.append(dict(**f,sha256=hashlib.sha256(out.read_bytes()).hexdigest()));continue
 f['ranges']=[(a,min(a+chunk-1,f['size']-1)) for a in range(0,f['size'],chunk)]
 remaining[f['name']]=len(f['ranges'])
 tasks.extend((f,a,b) for a,b in f['ranges'])
def get(t):
 f,a,b=t;dst=parts/(f['name']+'.'+str(a));hdr=Path(str(dst)+'.hdr')
 if dst.exists() and dst.stat().st_size==b-a+1 and hdr.exists() and f'content-range: bytes {a}-{b}/' in hdr.read_text().lower():return f,True
 for attempt in range(5):
  time.sleep(.3+attempt*3)
  rr=subprocess.run(['curl','--proxy','http://127.0.0.1:7897','-L','--fail','--connect-timeout','12','--max-time','120','-sS','-r',f'{a}-{b}','-D',str(hdr),f['url'],'-o',str(dst)],capture_output=True)
  if rr.returncode==0 and dst.stat().st_size==b-a+1 and f'content-range: bytes {a}-{b}/' in hdr.read_text().lower():return f,True
 return f,False
fails=[]
print('DOWNLOAD_REFERENCE',len(samples),'samples',len(tasks),'chunks',flush=True)
with concurrent.futures.ThreadPoolExecutor(6) as ex:
 for fut in concurrent.futures.as_completed([ex.submit(get,t) for t in tasks]):
  f,ok=fut.result()
  if not ok:fails.append(f['name']);continue
  remaining[f['name']]-=1
  if remaining[f['name']]==0:
   out=raw/f['name']
   with open(str(out)+'.assembling','wb') as dst:
    for a,b in f['ranges']:dst.write((parts/(f['name']+'.'+str(a))).read_bytes())
   with gzip.open(str(out)+'.assembling','rb') as z:
    while z.read(1024*1024):pass
   Path(str(out)+'.assembling').replace(out)
   records.append(dict(name=f['name'],size=f['size'],sha256=hashlib.sha256(out.read_bytes()).hexdigest()))
   (P/'metadata/reference_download_manifest.json').write_text(json.dumps(records,indent=2))
   print('DONE',f['name'],flush=True)
(P/'logs/reference_download_complete.json').write_text(json.dumps(dict(ok=not fails,files=len(records),failed=fails)))
