from pathlib import Path
import subprocess,concurrent.futures,json,re
r=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930/validation_v4')
for d in ['raw','scripts','metadata','results','logs']:(r/d).mkdir(parents=True,exist_ok=True)
urls={
'gse266873_listing.html':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE266nnn/GSE266873/suppl/',
'gse266873.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE266nnn/GSE266873/soft/GSE266873_family.soft.gz',
'gse333172_listing.html':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE333nnn/GSE333172/suppl/',
'gse333172_page.html':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333172',
'hpa.xml':'https://www.proteinatlas.org/ENSG00000198753.xml',
'plxnb3_uniprot.json':'https://rest.uniprot.org/uniprotkb/Q9ULL4.json'
}
def f(k,u):
 p=r/'raw'/k;z=subprocess.run(['curl','-L','--fail','--connect-timeout','10','--max-time','45','-sS',u,'-o',str(p)],capture_output=True,text=True)
 row=dict(file=k,url=u,exit=z.returncode,error=z.stderr,bytes=p.stat().st_size if p.exists() else 0)
 print(json.dumps(row),flush=True);return row
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(lambda kv:f(*kv),urls.items()))
(r/'metadata/downloads.json').write_text(json.dumps(rows,indent=2))
for p in (r/'raw').glob('*listing.html'):print(p.name,p.read_text()[:5000])
