from pathlib import Path
import json,pandas as pd,numpy as np,hashlib,zipfile,datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
b=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930');r=b/'validation_v3'
sources=[
('ITGA2B/ITGB3','2025 BMC Neurology','https://doi.org/10.1186/s12883-025-04408-w','Original GSE296792 manuscript explicitly names both genes and platelet-related changes.','High overlap; not a new mechanism'),
('Platelet/inflammation','2022 Frontiers Neurology','https://doi.org/10.3389/fneur.2022.886329','Public-data study already examined platelet activation and SP3 after ICH.','High thematic overlap'),
('Neutrophil/DEFA3 context','2025 BMC Neurology','https://doi.org/10.1186/s12883-025-04408-w','Day7 neutrophil degranulation is original main conclusion.','Replication alone not novel'),
('ALOX15','2025 Food Function','https://doi.org/10.1039/d5fo00425j','Experimental ICH ferroptosis work directly targets ALOX15.','Mechanistic competition exists'),
('Myeloid time/tissue programs','2021 Science Immunology','https://doi.org/10.1126/sciimmunol.abd6279','GSE163256 already mapped inflammatory to reparative stages in hematoma myeloid cells.','New target requires distinct hypothesis'),
('GRAMD1C','2022 Nature Communications','https://pmc.ncbi.nlm.nih.gov/articles/PMC9586981/','Primary experiments support cholesterol transport/autophagy function outside ICH.','Biological plausibility, not proof in ICH'),
('PLXNB3','Targeted searches 2026-09-30','https://www.proteinatlas.org/ENSG00000198753-PLXNB3/interaction','No direct ICH mechanistic study identified by this targeted search; this is not proof of absence.','Candidate novelty unresolved'),
('Day7 PBMC validation option','2026 Annals Clinical Translational Neurology','https://doi.org/10.1002/acn3.70485','7 sequenced patients: 5 BAF312 and 2 placebo; PBMC rather than whole blood.','Unsuitable as adequately replicated untreated whole-blood validation'),
('UC Davis blood/edema study','2021 Translational Stroke Research','https://doi.org/10.1007/s12975-020-00869-y','Availability upon reasonable request.','Not yet a directly downloadable validation cohort')
]
pd.DataFrame(sources,columns=['topic','study','url','evidence','decision']).to_csv(r/'results/novelty_evidence.csv',index=False)
queries=['"intracerebral hemorrhage" "ITGA2B" "ITGB3"','"intracerebral hemorrhage" "GRAMD1C"','"PLXNB3" "intracerebral hemorrhage"','"plexin-B3" hemorrhage','"PRSS33" hemorrhage stroke','"ALOX15" "intracerebral hemorrhage"','"intracerebral hemorrhage" "day 7" "whole blood" "GSE"']
(r/'metadata/search_record.json').write_text(json.dumps({'date':'2026-09-30','scope':'targeted novelty/data search, not systematic review','queries':queries},indent=2))
strict=pd.read_csv(r/'results/strict_day_matched_tissue.csv')
proxy=pd.read_csv(r/'results/discovery_proxy_adjusted_all.csv.gz');ext=pd.read_csv(r/'results/external_proxy_adjusted_all.csv.gz')
sp=pd.read_csv(r/'results/sorted_platelet_proxy_sensitivity.csv');tissue=pd.read_csv(r/'results/sorted_validation_tests.csv')
data=pd.read_csv(r/'results/sorted_patient_gene_values.csv');data['day']=data.bio_key.str.split('|').str[2]
fig,axs=plt.subplots(1,2,figsize=(9,4.2))
matched=[]
for ax,g in zip(axs,['PLXNB3','ITGB3']):
 d=data[(data.cell=='monos')&(data.stage=='early')&(data.gene==g)]
 z=d.groupby(['patientid','day','source'])[['logcpm','hrs after symptom onset']].mean().unstack('source').dropna()
 z=z[abs(z[('hrs after symptom onset','blood')]-z[('hrs after symptom onset','hematoma')])<=6]
 zz=z['logcpm'].groupby('patientid').mean()
 for pid,row in zz.iterrows():
  ax.plot([0,1],[row.blood,row.hematoma],color='#55869d',alpha=.55,marker='o',ms=3,lw=.8)
  matched.append(dict(gene=g,patientid=pid,blood=row.blood,hematoma=row.hematoma))
 a=strict[(strict.cell=='monos')&(strict.gene==g)&(strict.stage=='early')].iloc[0]
 ax.set_title(g+' | n='+str(int(a.n_pairs))+' | FDR='+format(a.fdr_28,'.4f'))
 ax.set_xticks([0,1],['Peripheral blood','Hematoma'])
 ax.set_ylabel('Mean TMM log2 CPM (prior count 1)')
 ax.spines[['top','right']].set_visible(False)
fig.suptitle('Sorted CD14+ cells: same-patient, same-day pairs (<96h)')
fig.tight_layout();fig.savefig(r/'results/strict_tissue_candidate.png',dpi=160);plt.close(fig)
pd.DataFrame(matched).to_csv(r/'results/figure_patient_values.csv',index=False)
targets=['ITGA2B','ITGB3','GRAMD1C','DEFA3','ALOX15','PRSS33','PLXNB3']
pc=proxy[proxy.gene.isin(targets)];pc.to_csv(r/'results/target_proxy_results.csv',index=False)
lines=[
'第三轮：新颖性审查与进一步验证',
'日期：2026-09-30；服务器：sgs10。',
'',
'结论与立项判断',
'1. 不建议把ITGA2B/ITGB3上调作为新机制主线直接立项：原始2025论文已在正文报告二者及血小板相关信号，2022年亦有相近公共数据分析。',
'2. 原有FDR显著事实保留。新增标记代理调整后证据减弱，说明不能把全血转录变化直接当成独立于细胞组成的激活机制。',
'3. 新候选为早期血肿CD14+单核/巨噬细胞中的PLXNB3富集，ITGB3可作为并列观察指标，不称为二者组成的调控轴。',
'4. PLXNB3属于有组织配对证据的探索性候选，尚未通过足以正式定题的细胞来源和机制验证。当前不宜承诺8–10分定位。',
'',
'一、已完成的新颖性核查',
'2025 BMC Neurology原文直接点名ITGA2B、ITGB3在72小时上调；其主结论亦包含第7天中性粒细胞脱颗粒。',
'2022 Frontiers in Neurology已做脑出血相关血小板激活与SP3的公共数据库分析。',
'ALOX15已有脑出血铁死亡实验研究，不能仅因本项目FDR显著就当成新靶点。',
'本轮针对性检索未识别到直接研究脑出血血肿髓系细胞PLXNB3的机制论文，但不是系统综述或“首次发现”证明。',
'GRAMD1C在胆固醇转运/自噬领域有实验基础，本轮独立细胞数据不支持其优先升级为脑出血髓系细胞机制主线。',
'详细文献、链接与判定见novelty_evidence.csv。',
'',
'二、全血标记代理敏感性分析',
'最终共同血小板标记：PF4、PPBP、GP1BA、GP9、TUBB1；粒细胞标记：CSF3R、FCGR3B、CEACAM8、S100A8、S100A9、FPR1。',
'外部矩阵缺少GP1BB，因此使用两队列共同可测的五基因血小板面板；保留原六基因版本供核查。',
'在原配对差值模型中加入两类标记分数的变化。发现队列仍按33,188基因×3对比共99,564项校正；外部按16,526基因校正。',
'ITGA2B、ITGB3在发现队列仍保留名义显著的早期/第7天时间效应，但完整家族FDR约0.303。',
'外部调整后ITGA2B P=0.1325、ITGB3 P=0.0549，二者FDR约0.981。',
'ITGA2B/ITGB3早期变化与血小板标记变化的Spearman相关分别约0.871和0.857。',
'这些标记只是表达代理，不是真实血小板或粒细胞计数；调整也会去除真实共享生物过程，不能据FDR消失断言“原结果是假阳性”或“完全由细胞数量解释”。',
'',
'三、新下载分选细胞数据及统计方法',
'GSE163256来自MISTIE III相关研究，包含分选CD14+单核/巨噬细胞及中性粒细胞的血液和血肿纵向样本。',
'按原研究质量标准筛选：FASTQ fragments不少于400万，检测基因数单核细胞不少于9000、中性粒细胞不少于8500；排除健康对照、缺失发病小时数及超过171小时的样本。',
'同患者、细胞类型、部位、研究日和发病小时数相同的技术重复先汇总计数，再用edgeR TMM标准化。',
'最终单核细胞193份文库合并为139份生物采样、21例患者；中性粒细胞189份文库合并为125份生物采样、17例患者。',
'固定七个候选ITGA2B、ITGB3、GRAMD1C、DEFA3、ALOX15、PRSS33、PLXNB3。',
'时间比较为发病<96小时与96–171小时；组织比较为同患者血肿与血液。两个问题各形成7基因×2细胞×2分层的28项检验家族。',
'每位患者在每个比较中只贡献一个均值或差值。使用双侧配对Wilcoxon检验和BH校正。CPM>=1的样本不足20%或缺少有效配对时不作推断检验，P/FDR标NA。',
'严格敏感性限制同采样日、发病小时数相差不超过6小时，先计算每日组织差值再按患者平均。',
'',
'四、严格组织配对结果',
strict[strict.gene.isin(['PLXNB3','ITGB3'])][['cell','gene','stage','n_pairs','detectable_fraction','median_change','fraction_up','p_raw','fdr_28']].to_string(index=False,float_format=lambda x:f'{x:.5g}'),
'早期PLXNB3：18例配对中15例血肿较高，中位TMM log2 CPM差+3.182，P=0.003205，FDR=0.02991。',
'早期ITGB3：18例中16例血肿较高，P=0.0003515，FDR=0.009843；晚期12例中11例同向，P=0.0009766，FDR=0.01367。',
'两者的分选中性粒细胞表达较低，不作可评估的组织差异宣称。',
'七候选的分选细胞早晚期变化在完整时间检验家族中均未达到FDR显著，因此上述是组织分区线索，不是第7天全血变化的直接复现。',
'',
'五、必须保留的验证限制',
'PLXNB3在严格早期分层中CPM>=1的采样比例为30%，血液端低表达较多，效应量对低表达处理和伪计数有依赖，不能直接宣称蛋白倍数变化。',
'分选单核细胞组织差值再加入血小板代理变化后，PLXNB3早期P=0.05439、FDR=0.50764。该敏感性模型的血小板差值全部为正，零代理差值的截距在观测范围之外，存在外推，不能单凭此模型断言PLXNB3信号为污染。',
'ITGB3早期在该模型下P=0.1832；晚期P=0.02794、FDR=0.50764。故不能称为已排除血小板附着/吞噬或混入的细胞内激活。',
'GRAMD1C在单核细胞中检测率过低，在可测的中性粒细胞中时间及组织比较均未支持显著变化。',
'本轮未落实到可直接复现第7天未治疗全血结果的另一充分重复队列。2026 BAF312研究虽有第1/3/7天单细胞，但仅2例安慰剂、5例治疗，且是PBMC，不能当作同等外部验证。',
'现有全血研究的37%感染定义是研究期使用抗生素/抗真菌药，属于汇总临床信息；既非明确颅内感染定义，也不是可用患者级感染配对标签。',
'不同公开项目不自动证明患者完全独立；这里新增数据支持另一临床研究的细胞/组织背景，不包装成同一结局的多中心验证。',
'',
'六、建议形成的下一版研究问题',
'候选题目：脑出血早期血肿与外周血CD14+髓系细胞的PLXNB3表达分区及其细胞来源。',
'核心需要回答：PLXNB3是否由某一真实髓系细胞状态表达，还是来自血小板相关RNA或细胞结合/吞噬；该分区是否能在独立单细胞或空间资料中重复。',
'达到这一关之前只保留候选，不将PLXNB3—ITGB3写成机制轴，不用单一FDR显著值代替新颖性与验证。',
'若客户要求现在就确定成熟发表方案，当前证据不足以支持原ITGA2B/ITGB3新机制方案；PLXNB3值得小范围后续核验，但尚未达到正式承诺阶段。',
'',
'来源',
'https://doi.org/10.1186/s12883-025-04408-w',
'https://doi.org/10.3389/fneur.2022.886329',
'https://doi.org/10.1039/d5fo00425j',
'https://doi.org/10.1126/sciimmunol.abd6279',
'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163256',
'https://doi.org/10.1002/acn3.70485',
'',
'文件与复现',
str(r),
'完整检验：sorted_validation_tests.csv、strict_day_matched_tissue.csv、sorted_platelet_proxy_sensitivity.csv。',
'全基因代理调整结果保留在discovery_proxy_adjusted_all.csv.gz和external_proxy_adjusted_all.csv.gz。',
'新增原始计数、GEO元数据、技术重复映射、脚本、质量标准和文献记录均在本目录中保存。'
]
(r/'results/novelty_validation_report_zh.txt').write_text('\n'.join(lines),encoding='utf-8')
manifest=[]
for p in (r/'raw').iterdir():
 if p.is_file() and '.seg' not in p.name and p.stat().st_size>0:
  manifest.append(dict(file=p.name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
(r/'metadata/verified_manifest.json').write_text(json.dumps(manifest,indent=2))
assert len(strict)==28 and len(tissue)==56 and len(proxy)==99564
summary={'status':'complete','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
'original_axis_novelty':'insufficient; directly discussed by source publication',
'new_candidate':'PLXNB3 early hematoma/blood CD14+ cell expression difference; exploratory',
'strict_PLXNB3':strict[(strict.gene=='PLXNB3')&(strict.cell=='monos')&(strict.stage=='early')].to_dict('records'),
'platelet_proxy_caveat':'Model extrapolates to unobserved zero proxy difference; neither proves nor rules out contamination.',
'sorted_temporal_fdr_hits':int(((tissue.test_family=='temporal_28')&(tissue.fdr<.05)).sum())}
(r/'results/final_summary.json').write_text(json.dumps(summary,indent=2))
(r/'logs/COMPLETE.json').write_text(json.dumps(summary,indent=2))
out=r/'novelty_validation_v3.zip'
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as z:
 for folder in ['results','scripts','metadata','logs']:
  for p in (r/folder).rglob('*'):
   if p.is_file():z.write(p,p.relative_to(r))
with zipfile.ZipFile(out) as z:assert z.testzip() is None
print('COMPLETE',out,out.stat().st_size)
