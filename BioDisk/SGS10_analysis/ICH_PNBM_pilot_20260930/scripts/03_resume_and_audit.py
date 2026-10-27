from pathlib import Path
import subprocess,gzip,xml.etree.ElementTree as E,json,pandas as pd
r=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930')
p=r/'raw/mrna.csv.gz.part'
log=(r/'logs/mrna_resume.log').open('w')
proc=subprocess.Popen(['curl','-L','--fail','-C','-','--connect-timeout','15','--max-time','900','-sS','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296792/suppl/GSE296792_Normalized_mRNA_counts_v1.csv.gz','-o',str(p)],stdout=log,stderr=log)
print('MRNA_DOWNLOAD_PID',proc.pid,flush=True)
df=pd.read_excel(r/'raw/ijms-3877263-supplementary.xlsx',header=1)
print('COLUMNS',df.columns.tolist(),flush=True)
print('GROUPS',df.groupby(['Database','group']).agg(samples=('Sample ID','size'),patients=('Patients ID','nunique')).to_json(),flush=True)
print('DIAGNOSES',df.groupby(['Database','Type of Diagnosis','group']).agg(samples=('Sample ID','size'),patients=('Patients ID','nunique')).to_json(),flush=True)
print('PATIENTS_GROUPS',df.groupby(['Database','Patients ID'])['group'].nunique().value_counts().to_dict(),flush=True)
print('CROSS_COHORT_IDS',set(df[df.Database==1]['Patients ID'])&set(df[df.Database==2]['Patients ID']),flush=True)
print('STUDY_SUBSET',df.loc[df.Database==2,['Sample ID','Patients ID','group','Time from admission, days','Positive CSF culture','p-HPhLA','PhLA','3ILA','IL_6']].to_json(orient='records'),flush=True)
rc=proc.wait();(r/'logs/mrna_download.exit').write_text(str(rc))
if rc==0:
 with gzip.open(p,'rb') as f:
  while f.read(1024*1024): pass
 p.rename(r/'raw/mrna.csv.gz')
print('MRNA_FINISHED',rc,flush=True)
