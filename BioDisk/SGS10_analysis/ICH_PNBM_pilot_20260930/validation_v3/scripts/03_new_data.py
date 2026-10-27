from pathlib import Path
import subprocess,re,concurrent.futures,gzip,json
r=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930/validation_v3')
# Preserve first proxy version. Use the common, outcome-independent 5-marker platelet panel in both datasets.
p=r/'scripts/02_proxy.R';s=p.read_text();(r/'scripts/02_proxy_six_platelet_markers.R').write_text(s)
s=s.replace('"GP1BA","GP1BB","GP9"','"GP1BA","GP9"');p.write_text(s)
if (r/'results/discovery_proxy_adjusted_all.csv.gz').exists():
 (r/'results/discovery_proxy_adjusted_all.csv.gz').rename(r/'results/discovery_proxy_adjusted_six_markers.csv.gz')
urls={
'sorted.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163256/soft/GSE163256_family.soft.gz',
'monos_counts.csv.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163256/suppl/GSE163256_monos_counts.csv.gz',
'neuts_counts.csv.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163256/suppl/GSE163256_neuts_counts.csv.gz',
'source_supp1.xlsx':'https://static-content.springer-cdn.com/esm/art%3A10.1186%2Fs12883-025-04408-w/MediaObjects/12883_2025_4408_MOESM1_ESM.xlsx'
}
def fetch(name,url):
 p=r/'raw'/name
 h=subprocess.run(['curl','-IL','--max-time','20','-sS',url],capture_output=True,text=True)
 sizes=re.findall(r'content-length:\s*(\d+)',h.stdout,re.I);size=int(sizes[-1]) if sizes else 0
 if size>1500000:
  def part(i):
   lo=size*i//4;hi=size*(i+1)//4-1;q=r/'raw'/(name+'.seg'+str(i));hd=r/'logs'/(name+'.seg'+str(i)+'.hdr')
   cp=subprocess.run(['curl','-L','--fail','--max-time','150','--connect-timeout','10','-sS','-r',f'{lo}-{hi}','-D',str(hd),url,'-o',str(q)],capture_output=True,text=True)
   assert cp.returncode==0 and q.stat().st_size==hi-lo+1 and f'bytes {lo}-{hi}/{size}' in hd.read_text(),cp.stderr
   return q
  with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:qs=list(ex.map(part,range(4)))
  with p.open('wb') as f:
   for q in qs:f.write(q.read_bytes())
 else:
  cp=subprocess.run(['curl','-L','--fail','--max-time','150','--connect-timeout','10','-sS',url,'-o',str(p)],capture_output=True,text=True)
  if cp.returncode:print(name,cp.returncode,cp.stderr,flush=True);return
 if name.endswith('.gz'):
  with gzip.open(p,'rb') as f:
   while f.read(1024*1024):pass
 print('DOWNLOADED',name,p.stat().st_size,flush=True)
for k,u in urls.items():fetch(k,u)
cp=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(r/'scripts/02_proxy.R')],capture_output=True,text=True)
(r/'logs/proxy.log').write_text(cp.stdout+cp.stderr);print('PROXY_FINAL',cp.returncode,cp.stdout,cp.stderr)
