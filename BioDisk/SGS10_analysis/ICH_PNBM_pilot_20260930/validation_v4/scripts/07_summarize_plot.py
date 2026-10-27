from pathlib import Path
import pandas as pd,numpy as np,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930/validation_v4')
d=pd.read_csv(p/'results/singlecell_sample_summary.csv');a=d[d.qc=='standard']
# Require all planned pilot matrices and analyses, keep single-cell results descriptive.
assert a.gsm.nunique()==5
selected=a[a.celltype.isin(['All','Myeloid','Myeloid_clean','Oligodendrocyte','Astrocyte','OPC','Platelet'])]
selected.to_csv(p/'results/singlecell_key_results.csv',index=False)
pos=[];maps=[]
for gsm in a.gsm.unique():
 c=pd.read_csv(p/'results'/(gsm+'_cell_audit.csv.gz'));c.insert(0,'gsm',gsm)
 pos.append(c[c.qc&(c.PLXNB3_UMI>0)])
 feat=next((p/'raw/singlecell').glob(gsm+'*features.tsv.gz'))
 f=pd.read_csv(feat,sep='\t',header=None);v=f[f[1]=='PLXNB3']
 assert len(v)==1 and v.iloc[0,0]=='ENSG00000198753'
 maps.append({'gsm':gsm,'gene':'PLXNB3','ensembl':v.iloc[0,0]})
positive=pd.concat(pos);positive.to_csv(p/'results/singlecell_positive_cells.csv',index=False)
pd.DataFrame(maps).to_csv(p/'results/gene_mapping_audit.csv',index=False)
h=a[a.series==333172];t=a[a.series==266873]
metrics={}
for key,sub in [('hematoma',h),('perihematoma',t)]:
 m={}
 for ct in ['All','Myeloid','Myeloid_clean','Oligodendrocyte','Astrocyte','OPC','Platelet']:
  v=sub[sub.celltype==ct];m[ct]={'cells':int(v.n_cells.sum()),'positive':int(v.PLXNB3_positive_cells.sum()),'ge2':int(v.PLXNB3_ge2_cells.sum()),'total_UMI':int(v.PLXNB3_total_UMI.sum())}
 metrics[key]=m
(p/'results/singlecell_metrics.json').write_text(json.dumps(metrics,indent=2))
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(1,3,figsize=(15,4.5),gridspec_kw={'width_ratios':[1,1.6,1]})
blood=pd.read_csv(p/'results/sorted_patient_blood.csv',index_col=0)
hem=pd.read_csv(p/'results/sorted_patient_hematoma.csv',index_col=0)
for pid in blood.index:ax[0].plot([0,1],[blood.loc[pid,'PLXNB3'],hem.loc[pid,'PLXNB3']],color='#AEB9C5',alpha=.7,lw=.9)
ax[0].scatter(np.zeros(len(blood)),blood.PLXNB3,color='#4278A8',s=20,zorder=3);ax[0].scatter(np.ones(len(hem)),hem.PLXNB3,color='#BF684A',s=20,zorder=3)
ax[0].set_xticks([0,1],['Blood','Hematoma']);ax[0].set_ylabel('Mean TMM log2 CPM (edgeR prior.count=1)');ax[0].set_title('Sorted CD14+ samples\n18 matched patients; P = 0.0032; FDR28 = 0.0299')
cts=['Myeloid','Oligodendrocyte','Astrocyte','OPC','Platelet'];colors=['#4278A8','#BF684A','#559575','#8C71AB','#BB9A41']
xpos=np.arange(5);width=.15
order=['GSM9758327','GSM9758337','GSM9758339','GSM9758345','GSM8255347']
for j,ct in enumerate(cts):
 vals=[]
 for gsm in order:
  row=a[(a.gsm==gsm)&(a.celltype==ct)].iloc[0]
  vals.append(row.positive_fraction*100 if row.n_cells else np.nan)
 ax[1].bar(xpos+(j-2)*width,vals,width=width,label=ct,color=colors[j])
ax[1].set_xticks(xpos,['H1','H6','H7','H10','PHE28h']);ax[1].set_ylabel('PLXNB3-positive cells (%)');ax[1].set_title('Independent single-cell source pilot\n4 hematoma patients + 1 perihematomal patient')
ax[1].legend(loc='upper left',bbox_to_anchor=(0,-.15),ncol=3,frameon=False,fontsize=7)
c=pd.read_csv(p/'results/clinical_patient_inputs.csv');mcol='modified rankin score at 1 year post-ich'
ax[2].scatter(c.PLXNB3,c[mcol],s=30,color='#4278A8');ax[2].set_xlabel('Earliest hematoma PLXNB3 (log2 CPM)');ax[2].set_ylabel('1-year mRS');ax[2].set_title('Clinical association: n = 19\nSpearman rho = 0.080; P = 0.745')
fig.suptitle('PLXNB3: bulk signal, cell-source evidence, and outcome check',fontsize=13)
fig.text(.02,.012,'Single-cell labels are provisional and marker-based. No bar can mean zero detection or no cells in that category; see sample tables.',fontsize=7,color='#555555');fig.tight_layout(rect=(0,.07,1,.93));fig.savefig(p/'results/PLXNB3_evidence.png',dpi=180);plt.close(fig)
print('METRICS',json.dumps(metrics));print('KEY',selected.to_json(orient='records'));print('POSITIVE LABELS',positive.groupby(['gsm','label']).size().to_json())
