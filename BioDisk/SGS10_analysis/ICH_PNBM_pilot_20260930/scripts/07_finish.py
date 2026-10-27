from pathlib import Path
import os,select,subprocess,sys,json
r=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930')
if not (r/'raw/mrna.csv.gz').exists():
 st=Path('/proc/3669155/stat')
 if st.exists():
  parent=int(st.read_text().split()[3]);fd=os.pidfd_open(parent)
  print('WAITING_ON_DOWNLOAD_COMPLETION',parent,flush=True)
  select.select([fd],[],[]);os.close(fd)
assert (r/'raw/mrna.csv.gz').exists(),'Download incomplete; analysis not started'
exe='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'
for name in ['05_ich_pilot.py','06_package.py']:
 p=subprocess.run([exe,str(r/'scripts'/name)],capture_output=True,text=True)
 (r/'logs'/(name+'.log')).write_text(p.stdout+p.stderr)
 print('STEP',name,'EXIT',p.returncode,flush=True)
 if p.returncode:
  print(p.stderr[-4000:],flush=True);sys.exit(p.returncode)
 if name.startswith('05'):
  j=json.loads((r/'results/ich_summary.json').read_text())
  print(json.dumps(j,indent=2),flush=True)
print('ALL_PILOTS_COMPLETE',flush=True)
