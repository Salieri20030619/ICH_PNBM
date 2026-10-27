from pathlib import Path
import os,subprocess,concurrent.futures,json,hashlib
root=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930')
print('BIODATA_OWN_DIRS',[(str(p),os.access(p,os.W_OK)) for p in Path('/mnt/BioData').iterdir() if p.stat().st_uid==os.getuid()],flush=True)
urls={
'geo_family.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296792/soft/GSE296792_family.soft.gz',
'mrna.csv.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296792/suppl/GSE296792_Normalized_mRNA_counts_v1.csv.gz',
'mirna.csv.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296792/suppl/GSE296792_Normalized_miRNA_counts_v1.csv.gz',
'infection_europepmc.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10319054/fullTextXML',
'metabolite_article.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12609026/fullTextXML',
'metabolite_correction.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12733275/fullTextXML',
'metabolite_supp.zip':'https://www.mdpi.com/1422-0067/26/21/10522/s1',
'metabolite_supp_alt.zip':'https://mdpi-res.com/d_attachment/ijms/ijms-26-10522/article_deploy/ijms-26-10522-s001.zip'
}
def fetch(k,u):
 p=root/'raw'/k
 r=subprocess.run(['curl','-L','--fail','--connect-timeout','10','--max-time','55','-sS',u,'-o',str(p)+'.part'],capture_output=True,text=True,timeout=60)
 row=dict(file=k,url=u,exit=r.returncode,error=r.stderr[:300])
 if r.returncode==0:
  Path(str(p)+'.part').rename(p);row.update(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
 print(json.dumps(row),flush=True);return row
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(lambda x:fetch(*x),urls.items()))
(root/'metadata/curl_download_manifest.json').write_text(json.dumps(rows,indent=2))
print('FETCH_COMPLETE',flush=True)
