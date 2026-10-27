from pathlib import Path
import subprocess,concurrent.futures,gzip,re,hashlib,json,sys
r=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930')
p=r/'raw/mrna.csv.gz.part';total=18344296;start=p.stat().st_size
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296792/suppl/GSE296792_Normalized_mRNA_counts_v1.csv.gz'
def fetch(i):
 lo=start+(total-start)*i//4;hi=start+(total-start)*(i+1)//4-1
 q=r/'raw'/('mrna_tail_'+str(i)+'.part');h=r/'logs'/('mrna_tail_'+str(i)+'.headers')
 cp=subprocess.run(['curl','-L','--fail','--connect-timeout','15','--max-time','240','-sS','-r',str(lo)+'-'+str(hi),'-D',str(h),url,'-o',str(q)],capture_output=True,text=True)
 ok=cp.returncode==0 and q.stat().st_size==hi-lo+1 and ('bytes '+str(lo)+'-'+str(hi)+'/'+str(total)) in h.read_text()
 row=dict(part=i,start=lo,end=hi,exit=cp.returncode,valid=ok,error=cp.stderr)
 print(json.dumps(row),flush=True)
 return row,q
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:parts=list(ex.map(fetch,range(4)))
(r/'metadata/tail_download_manifest.json').write_text(json.dumps([x[0] for x in parts],indent=2))
assert all(x[0]['valid'] for x in parts),'range integrity failure'
with p.open('ab') as out:
 for row,q in parts:out.write(q.read_bytes())
assert p.stat().st_size==total
with gzip.open(p,'rb') as f:
 while f.read(1024*1024):pass
p.rename(r/'raw/mrna.csv.gz')
(r/'logs/mrna_download_final.exit').write_text('0')
print('MRNA_INTEGRITY_PASS',flush=True)
exe='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'
for name in ['05_ich_pilot.py','06_package.py']:
 cp=subprocess.run([exe,str(r/'scripts'/name)],capture_output=True,text=True)
 (r/'logs'/(name+'.log')).write_text(cp.stdout+cp.stderr)
 print('STEP',name,'EXIT',cp.returncode,flush=True)
 if cp.returncode:
  print(cp.stderr[-4000:],flush=True);sys.exit(cp.returncode)
 if name.startswith('05'):print(cp.stdout,flush=True)
print('ALL_PILOTS_COMPLETE',flush=True)
