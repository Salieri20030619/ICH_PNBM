from pathlib import Path
import subprocess,concurrent.futures,json,re,hashlib
base=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930')
r=base/'validation_v2'
for d in ['raw','scripts','results','metadata','logs']:(r/d).mkdir(parents=True,exist_ok=True)
plan={
'status':'exploratory_extension_after_v1_not_confirmatory',
'ICH_primary':'Keep v1 panels and 12-test family unchanged; whole-transcriptome paired limma with all 3 temporal contrasts jointly BH-corrected; abundance-only filtering. Wilcoxon sensitivity with same universe.',
'ICH_replication':'GSE125512 older paired whole-blood cohort; analyze all available genes before selecting replication hits; verify 2019 vs 2025 source overlap and timing. Distinct accession is not alone proof of independent patients.',
'pathways':'Reactome full pathway universe 15-500 measured genes; paired sample scores for all three contrasts; BH across all pathways and contrasts. Exploratory.',
'infection':'All 11 metabolites plus 4 biomarkers = 15 features in the acute17 cohort; no deletion of high noninfected outlier; exact rank permutation; preserve v1 4-marker result. Chronic cohort separate sensitivity, never enlarge acute sample with chronic patients.',
'avoid':'No switching to one-sided tests, shrinking correction families, selective exclusions or repeated-sample pseudoreplication.',
'notes':'Assess downloaded external infection records for actual overlap before calling validation.'
}
(r/'metadata/analysis_plan.json').write_text(json.dumps(plan,indent=2))
urls={
'external.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125512/soft/GSE125512_family.soft.gz',
'external.xls.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125512/suppl/GSE125512_SecondSample_vs_FirstSample.all.gene.result.xls.gz',
'paired_infection.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12388189/fullTextXML',
'old_infection.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8955909/fullTextXML',
'paired_infection_supp.zip':'https://mdpi-res.com/d_attachment/metabolites/metabolites-15-00527/article_deploy/metabolites-15-00527-s001.zip',
'reactome.gmt.zip':'https://reactome.org/download/current/ReactomePathways.gmt.zip'
}
def fetch(k,u):
 p=r/'raw'/k
 z=subprocess.run(['curl','-L','--fail','--max-time','100','--connect-timeout','10','-sS',u,'-o',str(p)+'.part'],capture_output=True,text=True)
 row=dict(file=k,url=u,exit=z.returncode,error=z.stderr)
 if z.returncode==0:
  Path(str(p)+'.part').rename(p);row.update(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
 print(json.dumps(row),flush=True);return row
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(lambda kv:fetch(*kv),urls.items()))
(r/'metadata/download_manifest.json').write_text(json.dumps(rows,indent=2))
