from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,gzip,zipfile,datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
b=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930');r=b/'validation_v2'
g=pd.read_csv(r/'results/ich_all_gene_tests.csv.gz');e=pd.read_csv(r/'results/gene_replication_all.csv.gz')
loo=pd.read_csv(r/'results/gene_leave_one_patient_out.csv');paths=pd.read_csv(r/'results/pathway_discovery_all.csv');ep=pd.read_csv(r/'results/pathway_external_all.csv')
inf=pd.read_csv(r/'results/infection_all_analytes.csv');old=pd.read_csv(r/'results/original_panels_external.csv')
x=pd.read_csv(b/'raw/mrna.csv.gz',index_col=0);x.columns=[c.replace('_bas','_base') if c.endswith('_bas') else c.replace('_7da','_7day') if c.endswith('_7da') else c for c in x]
ext=pd.read_csv(r/'raw/external.xls.gz',sep='\t');ext=ext[~ext.symbol.duplicated(False)&ext.symbol.notna()].set_index('symbol')
fig,axs=plt.subplots(1,2,figsize=(9,4))
for pid in sorted(set(c.split('_')[0] for c in x)):
 vals=[np.log2(x.loc['ITGA2B',pid+'_'+t]+1)-np.log2(x.loc['ITGA2B',pid+'_base']+1) if pid+'_'+t in x else np.nan for t in ['base','72h','7day']]
 axs[0].plot([0,1,2],vals,color='#507f9d',alpha=.4,lw=1)
axs[0].set_xticks([0,1,2],['Baseline','72h','Day 7']);axs[0].set_title('GSE296792: 27 patients')
for pid in range(1,12):
 v=np.log2(ext.loc['ITGA2B',f'ICH{2*pid}']+1)-np.log2(ext.loc['ITGA2B',f'ICH{2*pid-1}']+1)
 axs[1].plot([0,1],[0,v],color='#b85a3d',alpha=.6,lw=1)
axs[1].set_xticks([0,1],['First draw','~72h after first']);axs[1].set_title('GSE125512: 11 patients')
for ax in axs:
 ax.axhline(0,color='gray',lw=.6);ax.spines[['top','right']].set_visible(False);ax.set_ylabel('Within-patient change in log2(value + 1)')
fig.suptitle('ITGA2B expression: discovery FDR, external nominal evidence')
fig.tight_layout();fig.savefig(r/'results/ITGA2B_validation.png',dpi=160);plt.close(fig)
early=e[(e.contrast=='72h_vs_base')&(e.fdr_limma_joint<.05)]
late=e[(e.contrast=='7day_vs_base')&(e.gene.isin(['ITGA2B','ITGB3','DEFA3','ALOX15','PRSS33']))]
robust=loo[loo.loo_fraction_fdr_significant==1]
selected=e[e.gene.isin(['ITGA2B','ITGB3','CD151','GRAMD1C'])].merge(loo[['gene','contrast','loo_fdr_max','loo_fraction_fdr_significant']],on=['gene','contrast'],how='left')
selected.to_csv(r/'results/recommended_gene_evidence.csv',index=False)
summary={
'gene_universe':g.gene.nunique(),'gene_tests':len(g),'gene_fdr_comparisons':int((g.fdr_limma_joint<.05).sum()),
'gene_fdr_unique':g.loc[g.fdr_limma_joint<.05,'gene'].nunique(),'both_methods_fdr_comparisons':int(((g.fdr_limma_joint<.05)&(g.fdr_wilcoxon_joint<.05)).sum()),
'all_loo_significant':len(robust),'pathway_tests':len(paths),'pathway_fdr_comparisons':int((paths.fdr_rank_joint<.05).sum()),
'external_gene_fdr':int((pd.read_csv(r/'results/external_gene_tests.csv.gz').fdr_limma<.05).sum()),
'external_pathway_fdr':int((ep.fdr_rank_all<.05).sum()),
'acute_infection_fdr':int(((inf.cohort=='acute_postoperative')&(inf.fdr_within_cohort_15<.05)).sum()),
'chronic_infection_fdr':int(((inf.cohort=='chronic_brain_injury')&(inf.fdr_within_cohort_15<.05)).sum()),
'overlap':json.loads((r/'results/infection_v2_summary.json').read_text()),
'status':'completed_exploratory_validation'
}
(r/'results/final_summary.json').write_text(json.dumps(summary,indent=2))
lines=[
'第二轮验证：FDR显著方向与可重复性',
'2026-09-30；服务器sgs10；探索性延伸，非预注册确认性研究。',
'',
'结论',
'脑出血可以找到FDR显著方向。较适合继续的是ITGA2B/ITGB3相关血小板黏附表达变化；FDR证据目前在基因层面，不能写成整条通路已经校正后显著。',
'另有第7天DEFA3等粒细胞抗菌相关基因信号，但没有同一时间窗的外部第7天队列。',
'现有术后感染17例未达到全指标FDR显著；不建议为达到阈值缩小检验家族或删除高值非感染患者。',
'',
'一、数据与方法',
'保留第一轮全部结果。新分析固定33,188个通过表达量与非零变异过滤的基因，三个时间对比共99,564项检验联合BH校正。',
'以患者为独立单位，用log2(标准化表达+1)的配对差值做limma经验贝叶斯稳健方差估计，配对Wilcoxon作为敏感性分析，分别报告对应完整家族FDR。',
'逐一剔除27例患者，每次对全部99,564项重新拟合和校正，保留所有43项初始FDR候选的稳定性结果。',
'外部GSE125512包含11例患者22份FPKM表达数据，按GEO注释确认ICH奇数/偶数列为同一患者的前后采样。剔除缺失/重复基因符号，按表达量过滤后对16,526个基因统一校正。',
'两个来源共同覆盖14,513个基因。Reactome取共同可测基因数15至500的全部1,231条通路，计算样本平均log2表达分数。',
'发现队列通路的三个配对时间对比共3,693项联合BH校正；外部1,231条通路另作完整家族校正。秩检验为通路主分析，t检验为敏感性分析。',
'原四面板在外部仅有一个时间对比，四项完整校正；第一轮12项发现队列FDR不变。',
'感染扩展到11个代谢物及4个生物标志物。急性术后17例按患者中位数；慢性脑损伤35个可识别患者单位另取最早采样，分别报告。',
'缺失或恒定指标不进行检验，原始P和FDR记NA；BH计算为保守保留计划家族大小15而使用非显著占位，不把占位值报告为检验P值。',
'两个感染队列各15项构成不同人群的检验家族，并额外提供全部30项合并校正作为更保守敏感性结果。',
'',
'二、脑出血：发现队列存在稳定FDR信号',
f"共{summary['gene_fdr_comparisons']}个基因×时间对比通过limma FDR<0.05，涉及{summary['gene_fdr_unique']}个不同基因；43项均同时通过秩检验完整家族FDR。",
'其中72小时对基线7项、第7天对72小时1项、第7天对基线35项。',
'8个第7天信号在逐患者剔除后仍全部保持FDR<0.05：',
robust[['gene','contrast','median_log2_change','fraction_increased','fdr_limma_joint','loo_fdr_max']].to_string(index=False,float_format=lambda z:f'{z:.5g}'),
'',
'较可比早期时间窗的跨数据集证据如下（外部仍仅名义显著）：',
early[['gene','median_log2_change_discovery','fraction_increased_discovery','p_limma_discovery','fdr_limma_joint','median_log2_change_external','p_limma_external','fdr_limma']].to_string(index=False,float_format=lambda z:f'{z:.5g}'),
'',
'ITGA2B第7天对基线：24对，23例升高，中位log2(value+1)变化+1.230，limma FDR=0.00503；逐患者剔除后的最差FDR=0.01511。',
'ITGA2B早期26对中24例升高，FDR=0.02308。外部11例中方向一致且P=0.00577，但全基因FDR=0.39161。',
'ITGB3早期FDR=0.03368，第7天FDR=0.03214；外部P=0.00565，FDR=0.39161。早期基因信号的留一法FDR会跨越0.05，详见完整表。',
'',
'三、必须保留的验证阴性结果',
'外部11例队列在全基因和完整通路家族下均无FDR<0.05结果，不能称为两个独立队列均校正后显著。',
'原抗原呈递面板在外部P=0.27832，四面板FDR=0.55664，未复现名义显著结果。因此撤下其优先主线地位。',
'发现队列Reactome通路有16个通路×时间对比通过秩检验FDR，但均未在外部复现名义显著；其t检验敏感性也未给出对应FDR支持。',
'血小板聚集通路早期FDR约0.05212，未达0.05。ITGA2B/ITGB3的基因级阳性不等于整条通路校正后显著。',
'',
'四、术后感染：当前数据尚不能达到稳定FDR证据',
inf[['cohort','feature','n_infected','n_noninfected','median_ratio','p_raw','fdr_within_cohort_15','fdr_all_30_sensitivity']].to_string(index=False,float_format=lambda z:f'{z:.5g}'),
'急性术后17例全指标FDR阴性；保留第一轮PhLA、3ILA的名义显著候选定位。',
'慢性脑损伤35个患者单位（最早采样感染7、未感染28）中，p-HPhLA和PhLA的中位数比为5.63和5.46，15项家族FDR均为0.03547。',
'将两人群全部30项合并校正时两者FDR为0.07094。逐一剔除慢性队列患者后最差FDR约0.10361，提示稳定性有限。',
'慢性脑损伤人群不能替代急性术后队列。该结果仅作为继发性中枢感染的探索性线索。',
'新下载脑脊液/血清配对补充表共提取29份样本；其中23份样本的11项脑脊液代谢物可与现有资料完全匹配，涉及14个现有患者单位，不能作为独立验证。',
'',
'五、立项建议与边界',
'建议脑出血优先：拟题可聚焦脑出血后血小板黏附相关表达程序的纵向变化，以ITGA2B/ITGB3为候选基因，不宣称已经证实细胞激活或因果机制。',
'目前FDR问题已有答案，下一关是新颖性与临床解释。两个原始研究已经讨论血小板相关基因和粒细胞反应，不能把复现已有结论直接作为新发现。',
'发现队列第7天结果的稳定性较好，但外部仅提供早期窗口；晚期结果不构成同时间窗外部验证。',
'两脑出血来源来自同一研究团队/机构，项目号不同但公开信息不足以完全排除患者重叠，不作跨中心独立验证宣传。',
'发现队列标准化counts与外部FPKM处理方法不同，只比较各自队列内方向与证据，不合并绝对表达量。',
'本轮是看到第一轮结果后的探索性延伸；保留全部分析，不能包装成预先注册的确认性检验。',
'感染方向不建议用当前数据单独承诺完整发表方案；若保持急性术后人群，需要新的可核查、患者独立的数据。',
'',
'来源与位置',
str(r),
'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296792',
'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125512',
'https://reactome.org/download/current/ReactomePathways.gmt.zip',
'https://doi.org/10.1186/s12883-025-04408-w',
'https://doi.org/10.3390/ijms262110522',
'https://doi.org/10.3390/metabo15080527',
'核心表：recommended_gene_evidence.csv、gene_leave_one_patient_out.csv、pathway_replication.csv、infection_all_analytes.csv、infection_external_overlap.csv。'
]
(r/'results/validation_report_zh.txt').write_text('\n'.join(lines),encoding='utf-8')
manifest=[]
for p in (r/'raw').iterdir():
 if p.is_file() and not p.name.endswith('.part') and not p.name.startswith('ext_tail'):
  manifest.append(dict(file=p.name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
(r/'metadata/verified_manifest.json').write_text(json.dumps(manifest,indent=2))
assert len(g)==99564 and len(paths)==3693 and len(ep)==1231 and len(inf)==30
assert not (inf.loc[inf.n_infected<2,'p_raw'].notna().any())
(r/'logs/COMPLETE.json').write_text(json.dumps({'status':'complete','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'exploratory FDR and replication validation'},indent=2))
out=r/'validation_v2_results.zip'
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as z:
 for folder in ['results','scripts','metadata','logs']:
  for p in (r/folder).rglob('*'):
   if p.is_file():z.write(p,p.relative_to(r))
with zipfile.ZipFile(out) as z:assert z.testzip() is None
print(json.dumps({k:v for k,v in summary.items() if k!='overlap'},indent=2))
print('PACKAGE',out,out.stat().st_size)
