from pathlib import Path
import json,zipfile,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930/validation_v5');O=P/'results/extension'
d=pd.read_csv(O/'all_sample_celltype_metrics.csv')
time=pd.read_csv(O/'temporal_tests.csv');co=pd.read_csv(O/'program_correlations.csv');within=pd.read_csv(O/'within_ol_state_tests.csv');src=pd.read_csv(O/'source_localization_tests.csv');brain=pd.read_csv(O/'brain_array_tests.csv');pairs=pd.read_csv(O/'independent_pair_tests.csv')
ref=pd.read_csv(O/'published_annotation_plxnb3.csv')
ol=d[(d.series==266873)&(d.qc=='standard')&(d.celltype=='Oligodendrocyte')].copy()
comparison=ol[['gsm','n_cells','PLXNB3_total_UMI','pseudobulk_CPM']].merge(ref[ref.cell_type.str.lower().str.contains('oligodendro')],left_on='gsm',right_on='sample_id',how='outer')
comparison.to_csv(O/'alternative_annotation_comparison.csv',index=False)
plt.rcParams.update({'font.size':9,'axes.spines.right':False,'axes.spines.top':False})
fig,ax=plt.subplots(1,3,figsize=(15,4.5))
colors={'Myeloid':'#4278A8','Oligodendrocyte':'#bf684a','Astrocyte':'#559575'}
for ct,color in colors.items():
 a=d[(d.series==266873)&(d.qc=='standard')&(d.celltype==ct)]
 ax[0].scatter(a.hours,a.positive_fraction*100,label=ct,color=color,s=35)
ax[0].set(xlabel='Hours after ICH',ylabel='PLXNB3-positive cells (%)',title='Cell source across 9 patients');ax[0].legend(frameon=False,fontsize=8)
ax[1].scatter(ol.hours,ol.pseudobulk_CPM,color='#bf684a',s=35)
r=time[time.qc=='standard'].iloc[0]
ax[1].set(xlabel='Hours after ICH',ylabel='OL pseudobulk PLXNB3 CPM',title=f'Patient-level time association\nrho={r.rho:.3f}; exact P={r.p_raw:.4g}')
b=pd.read_csv(O/'brain_array_patient_values.csv')
for _,row in b.iterrows():ax[2].plot([0,1,2],[row.CW,row.PH,row.CG],marker='o',alpha=.7)
ax[2].set_xticks([0,1,2],['Contralateral WM','Perihematoma','Contralateral GM']);ax[2].set(ylabel='PLXNB3 log2 RMA',title='Paired brain tissue: 4 patients')
fig.text(.015,.01,'Exploratory patient-level results. Provisional cell labels; cross-sectional time groups. Grey/white matter are separate comparators.',fontsize=8)
fig.tight_layout(rect=(0,.06,1,1));fig.savefig(O/'PLXNB3_extension_overview.png',dpi=180);plt.close(fig)
def table(x):return x.to_string(index=False,float_format=lambda v:f'{v:.6g}')
text=f"""PLXNB3 扩展研究离线结果记录
本结果由服务器离线任务生成；需结合细胞注释、稀疏检测和患者差异进行人工审阅，不自动形成机制或发表结论。

数据与完成范围
GSE266873：9 位患者、9 份血肿周围组织；GSE333172：6 位组织标签一致患者的 12 份配对样本。
共 21 份单细胞样本。另外读取 GSE24265 的 4 位患者、11 份脑组织芯片样本；灰质对照 4 对，白质对照 3 对，均为患者内对照，并非健康人。
公开 Zenodo 复核包的冻结伪总体用于同一 GSE266873 数据的注释敏感性，不计为新的生物学队列。

方法
沿用上一轮固定的单细胞 QC 和经典标记大类注释；PLXNB3 不参与分类。分别保留标准和严格 QC 结果。
患者级统计，不将细胞数作为样本量。分析中所有 n<10 的结果均为探索性。
时间分析：少突胶质细胞至少 30 个的样本，以患者伪总体 PLXNB3 CPM 对发病后小时做 Spearman 精确置换；单个预设时间终点，不需 FDR。
六个预先冻结的状态程序不含 PLXNB3；其患者级相关在 6 项家族内 BH 校正。程序是人工整理的标记分数，不等同于通路活性或因果作用。
细胞内状态对比：各患者少突胶质细胞按 UMI 深度五等分，比较同层 PLXNB3 阳性和阴性细胞的程序分数，再按阳性细胞数加权形成一个患者效应；患者至少 5 个阳性、20 个阴性，层内至少 10 个阴性。对患者效应做 Wilcoxon，6 项 BH 家族。
来源比较为少突胶质细胞对髓系、对星形胶质细胞，共 2 项 BH 家族。严格 QC 是同一问题的敏感性结果，不视为额外独立证据。
血肿—血液配对按 4 个集合保留检验家族；少于 3 对，或少于 3 位患者任一组织达到 CPM≥1，则不作推断性检验。
脑组织芯片仅使用明确映射到 PLXNB3 的 205957_at 探针，使用已经 log2 RMA 标准化的值，分别计算 PH−CW 和 PH−CG；2 项 BH 家族。
筛选遵从原始 P<0.05；仅实际 FDR<0.05 才称校正后显著。全零、不可评估或缺少重复不补造 P。

患者级细胞来源
{table(src)}

时间关联
{table(time)}

六个状态程序相关
{table(co)}

控制 UMI 深度后的患者级状态对比
{table(within)}

独立血肿—血液配对
{table(pairs)}

第三队列脑组织参考区敏感性
{table(brain)}

边界
时间组来自不同患者，不是纵向随访；年龄、性别和严重度只有部分组级信息，不能用细胞量弥补。
9 人以内的秩检验存在离散 P 值下限；原始 P 或 FDR 阴性不单独否定一致的方向和效应。
人工标记分数、共表达、RNA 共检出及计算模型均不能证明细胞自主表达、吞噬来源或信号通路被激活。
全部阳性细胞含混合谱系标记的记录保留，未以严格筛选后的零检出替代原始观察。
灰质与白质的组织构成不同，不能混池来改变 PLXNB3 的结论。
公开复核包涉及相同数据集但其文本/表格检索未出现 PLXNB3；这不等于该基因或细胞谱系研究从未发表。

文件
完整输入与结果位于 results/extension；原始矩阵保留在 validation_v4/raw/singlecell。
无网络任务队列、输入校验、每步日志和退出码位于 scripts/offline_queue.py、metadata/offline_plan.json、logs/offline_queue。
"""
(O/'PLXNB3_extension_auto_report_zh.txt').write_text(text,encoding='utf-8')
with zipfile.ZipFile(P/'PLXNB3_extension_offline_results.zip','w',zipfile.ZIP_DEFLATED) as z:
 for folder in [O,P/'metadata',P/'scripts']:
  for f in folder.rglob('*'):
   if f.is_file() and '__pycache__' not in str(f):z.write(f,f.relative_to(P))
with zipfile.ZipFile(P/'PLXNB3_extension_offline_results.zip') as z:assert z.testzip() is None
print('Report and result package complete')
