from pathlib import Path
import json,hashlib,gzip,zipfile,datetime,pandas as pd
r=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930')
a=json.loads((r/'results/ich_summary.json').read_text())
b=json.loads((r/'results/infection_summary.json').read_text())
ich=pd.read_csv(r/'results/ich_module_pilot.csv')
inf=pd.read_csv(r/'results/infection_marker_pilot.csv')
lines=[
'脑出血与术后颅内感染：公共数据小规模验证',
'日期：2026-09-30；执行服务器：sgs10',
'',
'一、选定方向',
'1. 脑出血：发病后第一周抗原呈递的时序变化及其与中性粒细胞相关表达的关系。',
'   数据：GSE296792的纵向外周血mRNA、miRNA及GEO样本元数据。',
'2. 术后颅内感染：神经外科术后脑脊液芳香族代谢物与感染状态的关联。',
'   数据：IJMS 2025论文10.3390/ijms262110522公开Table S1中的急性神经外科术后队列。',
'   该队列全部为肿瘤术后患者；不将慢性脑损伤患者混入同一个术后主队列。',
'',
'二、数据与方法',
'脑出血：以患者为独立单位。处理后标准化表达量转为log2(value+1)，每个预设基因组取平均表达分数。',
'基因表达以所有基线样本的对应基因中位数中心化；该常数不改变患者内时间差值。',
'预设抗原呈递、干扰素反应、中性粒细胞相关三个六基因面板，以及抗原呈递减中性粒细胞分数。',
'三种时间对比均按患者ID配对，采用双侧Wilcoxon符号秩检验。12个面板对比为一个BH校正家族。',
'四个既有miRNA候选成熟臂的12个复现对比独立作为另一个BH家族。',
'感染：仅取Database=2，以队列号和患者号建立唯一标识，重复样本按患者取中位数。',
'四项文献预选指标为p-HPhLA、PhLA、3ILA、IL-6。用全部标签排列的双侧秩和置换检验，四项P值进行BH校正。',
'报告组间中位数比值、描述性AUC及2000次组内bootstrap区间；用删除单个患者和首个来源记录进行敏感性检查。',
'所有有效筛选以原始P<0.05为名义显著门槛，同时完整保留FDR、效应量与样本一致性。',
'',
'三、脑出血结果',
f"完整下载mRNA矩阵：{a['mrna_rows']}行 × {a['mrna_samples']}份样本；miRNA矩阵：{a['mirna_rows']}行 × {a['mirna_samples']}份样本。",
f"独立患者{a['n_patients']}例，三时间点齐全者{a['complete_three_time_patients']}例。基线27份、72小时26份、7天24份。",
f"患者内跨测序批次人数：{a['patients_with_multiple_batches']}。",
f"12项预设面板对比中，原始P<0.05有{a['nominal_panel_tests']}项，FDR<0.05有{a['fdr_panel_tests']}项。",
'全部预设面板结果如下；变化单位为平均log2表达分数：',
ich[['feature','contrast','n_pairs','median_change','fraction_increased','p_raw','fdr_bh']].to_string(index=False,float_format=lambda x:f'{x:.5g}'),
'',
'抗原呈递72小时对基线：26对，20例升高，中位变化+0.349，95% bootstrap区间0.143至0.628，P=0.00729，FDR=0.06910。',
'第7天对基线：24对，17例升高，中位变化+0.534，区间0.0849至0.6964，P=0.01152，FDR=0.06910。',
'抗原呈递减中性粒细胞分数仅72小时对基线名义显著（P=0.04633，FDR=0.18532）；未建立持续不同步变化。',
'当前预设干扰素与中性粒细胞六基因面板的所有对比均未达原始P<0.05，不列入显著候选。该小面板不等同于原论文的完整中性粒细胞脱颗粒通路。',
'',
'四、术后感染结果',
'公开补充表为96行记录；论文报告53例患者，但表内按队列与患者号仅能识别52个单位。',
'急性术后队列19份样本、17例患者，其中感染10例、未感染7例。四项指标均无缺失。',
'原始病名为Tumors，故主分析人群定义为肿瘤神经外科术后。按原文临床标准保留感染标签。',
inf[['marker','n_infected','n_noninfected','median_ratio','p_raw','fdr_bh','auc_descriptive','auc_bootstrap_lo','auc_bootstrap_hi']].to_string(index=False,float_format=lambda x:f'{x:.5g}'),
'PhLA与3ILA达到原始P<0.05，为名义显著信号；四项检验均未达到FDR<0.05。',
'p-HPhLA与IL-6在本次患者级检验中未达原始P<0.05，不列入名义显著候选。',
'上述为同一已发表来源的再分析和可计算性验证，AUC为本队列描述性区分度，未经独立验证。',
'',
'五、联合与配对数据核查',
'出血原发病标签Haemorrage仅见于慢性队列：8例患者，其中感染3例、未感染5例。',
'该标签未细分脑实质出血、蛛网膜下腔出血，不能直接将8例归为自发性脑实质出血。',
'存在感染标签随采样变化的1例患者，其原发病为肿瘤；本表未发现出血患者感染前后配对。',
'因此本次保留两个独立方向；不将这8例作为联合主队列，不作联合组推断性检验。',
'',
'六、证据边界与立项判断',
'脑出血已具备纵向配对和可复算数据基础，可继续作为优先方向。自定义面板仅作探索性表达摘要。',
'全血组成变化可影响面板分数，不能据此确认细胞内机制；当前下载元数据没有核实感染或预后标签。',
'后续需要独立数据检验方向与细胞来源，且必须与原研究的中性粒细胞脱颗粒/miRNA结论区分。',
'术后感染的代谢物信号可以继续作为候选方向，但17例规模不足以独立支持复杂预测模型。',
'术后19份样本的采样时间全部缺失，不能称为感染前预警；培养阳性仅1例，不能做可靠的培养确证敏感性推断。',
'IL-6有8份样本恰为5000，存在检测上限截断疑虑，不能简单解释为真实相等浓度。',
'慢性队列存在患者数与论文不一致、两个样本号重复用于不同日期、一个患者号跨队列复用等问题，已输出完整记录，未擅自修正原始值。',
'原论文更正涉及参考浓度换算及相应文字；本次不使用这些参考阈值作诊断或分组。',
'数据中大部分候选指标由原论文选出，当前复算不构成新发现；独立验证和新问题是正式立项前提。',
'本次验证不支持直接承诺3分或8–10分发表。',
'',
'七、文件位置与复现',
str(r),
'原始数据raw/；样本映射metadata/；脚本scripts/；统计表与图results/；退出码logs/。',
'核心结果：ich_module_pilot.csv、ich_mirna_reproduction.csv、infection_marker_pilot.csv、infection_sensitivity.csv。',
'执行环境：/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python。',
'原始矩阵下载齐全后执行scripts/04_infection_pilot.py与scripts/05_ich_pilot.py。',
'本项目仅约几十MB。BioData约5.6TiB空闲，但当前账号无自己的可写目录，故未在那里创建或更改权限。',
'',
'八、来源',
'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296792',
'https://doi.org/10.1186/s12883-025-04408-w',
'https://doi.org/10.3390/ijms262110522',
'https://doi.org/10.3390/ijms262411856',
'https://mdpi-res.com/d_attachment/ijms/ijms-26-10522/article_deploy/ijms-26-10522-s001.zip',
'此前Olink论文10.3389/fneur.2023.1166598的EuropePMC补充文件清单仅列相关系数图片，未提供可用患者级矩阵，因此本次改用上述代谢物来源。'
]
(r/'results/pilot_report_zh.txt').write_text('\n'.join(lines),encoding='utf-8')
manifest=[]
for name in ['mrna.csv.gz','mirna.csv.gz','geo_family.soft.gz','metabolite_supp_alt.zip','ijms-3877263-supplementary.xlsx','metabolite_article.xml','metabolite_correction.xml','infection_europepmc.xml']:
 p=r/'raw'/name
 if name.endswith('.gz'):
  with gzip.open(p,'rb') as f:
   while f.read(1024*1024):pass
 if name=='mrna.csv.gz':assert p.stat().st_size==18344296
 manifest.append(dict(file=str(p.relative_to(r)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
(r/'metadata/verified_files.json').write_text(json.dumps(manifest,indent=2))
assert len(ich)==12 and len(inf)==4
assert (r/'logs/ich_pilot.exit').read_text()=='0'
assert (r/'logs/infection_pilot.exit').read_text()=='0'
zipout=r/'pilot_results_20260930.zip'
with zipfile.ZipFile(zipout,'w',zipfile.ZIP_DEFLATED) as z:
 for dirname in ['results','metadata','scripts','logs']:
  for p in (r/dirname).rglob('*'):
   if p.is_file():z.write(p,p.relative_to(r))
(r/'logs/validation_complete.json').write_text(json.dumps({'exit_status':0,'status':'complete','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'archive':str(zipout),'scope':'small-scale pilots only'},indent=2))
print('VALIDATION_COMPLETE',zipout,zipout.stat().st_size)
