from pathlib import Path
import requests,concurrent.futures,subprocess,json,hashlib,time
root=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930')
for d in ['raw','metadata','scripts','results','logs']: (root/d).mkdir(parents=True,exist_ok=True)
urls={
'geo_page.html':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296792',
'geo_family.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296792/soft/GSE296792_family.soft.gz',
'mrna.csv.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296792/suppl/GSE296792_Normalized_mRNA_counts_v1.csv.gz',
'mirna.csv.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296792/suppl/GSE296792_Normalized_miRNA_counts_v1.csv.gz',
'infection_article.html':'https://www.frontiersin.org/journals/neurology/articles/10.3389/fneur.2023.1166598/full',
'infection_europepmc.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10319054/fullTextXML'
}
def fetch(k,u):
 p=root/'raw'/k
 try:
  r=requests.get(u,timeout=45);r.raise_for_status();p.write_bytes(r.content)
  return dict(file=k,url=u,bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest(),status='downloaded')
 except Exception as e:return dict(file=k,url=u,status='error',error=str(e))
rows=list(concurrent.futures.ThreadPoolExecutor(max_workers=3).map(lambda kv:fetch(*kv),urls.items()))
(root/'metadata/download_manifest.json').write_text(json.dumps(rows,indent=2))
print(json.dumps(rows,indent=2),flush=True)
base=Path('/mnt/BioDisk/wangmingyuan/miniconda3')
for p in [base/'bin/python']+list((base/'envs').glob('*/bin/python')):
 try:
  r=subprocess.run([str(p),'-c','import pandas,scipy,numpy,matplotlib; print(pandas.__version__)'],capture_output=True,text=True,timeout=12)
  if r.returncode==0: print('ANALYSIS_PYTHON',str(p),r.stdout.strip())
 except Exception:pass
