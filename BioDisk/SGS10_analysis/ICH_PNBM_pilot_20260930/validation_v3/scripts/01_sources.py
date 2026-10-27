from pathlib import Path
import subprocess,concurrent.futures,json,re,xml.etree.ElementTree as E
b=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930');r=b/'validation_v3'
for d in ['raw','scripts','results','metadata','logs']:(r/d).mkdir(parents=True,exist_ok=True)
plan={'purpose':'Novelty and biological interpretation audit, preserving v1 and v2',
'proxy_sensitivity':'Fixed platelet panel PF4 PPBP GP1BA GP1BB GP9 TUBB1 and neutrophil panel CSF3R FCGR3B CEACAM8 S100A8 S100A9 FPR1; no ITGA2B/ITGB3 in proxies. Compare paired gene differences adjusting changes in both proxies. Same full v2 gene/time correction family. Proxies are not measured counts and can remove shared biology.',
'new_data':'Look for separate institution, comparable whole-blood time course and true day7 replication; sorted-myeloid or brain tissue only supports cellular context, not platelet temporal replication.',
'novelty':'Original source overlap must be documented, including already named ITGA2B/ITGB3 and neutrophil programs. No claim of novel mechanism from statistical reanalysis.'}
(r/'metadata/plan.json').write_text(json.dumps(plan,indent=2))
urls={x+'.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/'+x+'/fullTextXML' for x in ['PMC12465859','PMC6909779','PMC8142129','PMC10019834','PMC6727143','PMC9201407']}
urls['GSE163256_listing.html']='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163256/suppl/'
def f(k,u):
 p=r/'raw'/k;z=subprocess.run(['curl','-L','--fail','--max-time','45','--connect-timeout','10','-sS',u,'-o',str(p)],capture_output=True,text=True)
 print(k,z.returncode,flush=True)
 return dict(file=k,url=u,exit=z.returncode,error=z.stderr)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(lambda kv:f(*kv),urls.items()))
(r/'metadata/downloads.json').write_text(json.dumps(rows,indent=2))
for p in (r/'raw').glob('*.xml'):
 try:
  t=E.parse(p);s=' '.join(t.getroot().itertext())
  accessions=sorted(set(re.findall(r'GSE\d+|phs\d+|E-MTAB-\d+',s)))
  print('ACCESSIONS',p.name,accessions)
  for sec in t.iter('sec'):
   title=sec.find('title')
   if title is not None and any(k in ''.join(title.itertext()).lower() for k in ['availability','accessibility','subjects','patients']):
    st=' '.join(sec.itertext())
    if len(st)<10000:print('SECTION',p.name,st[:3500])
 except Exception as e:print(type(e).__name__)
