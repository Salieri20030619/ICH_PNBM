from pathlib import Path
import pandas as pd,json,zipfile,hashlib,platform
p=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930/validation_v4')
met=json.loads((p/'results/singlecell_metrics.json').read_text())
sc=pd.read_csv(p/'results/singlecell_sample_summary.csv')
ss=sc[sc.qc=='standard']
h=met['hematoma'];t=met['perihematoma']
tab=ss[ss.celltype=='All'][['gsm','title','n_cells','PLXNB3_positive_cells','PLXNB3_total_UMI']]
positive=pd.read_csv(p/'results/singlecell_positive_cells.csv')
lo=pd.read_csv(p/'results/sorted_leave_one_patient_out.csv')
strict=sc[(sc.qc=='strict')&(sc.series==333172)&(sc.celltype=='Myeloid')]
text=f"""PLXNB3 新颖性及进一步验证报告
2026-09-30 | SGS10 | 公共数据库调查及小规模验证

结论
PLXNB3 在脑出血后的全血时间变化及血肿—血液分选样本对比中存在统计信号，值得保留为候选。当前证据不支持把它直接确定为“单核/微胶质细胞内 PLXNB3 机制轴”，也不支持预后标志物主线。本轮独立单细胞小规模核查提示优先追查少突胶质细胞来源。建议将“机制主靶点”降为“细胞来源待确证的候选”，暂不据此形成 8–10 分级别机制方案或发表承诺。

方法
既有全血发现队列 GSE296792 使用 log2(x+1) 患者内差值和 limma；33,188 个基因 × 3 个对比共 99,564 项为同一 BH 家族。本轮读取原始结果，未缩小校正池。
GSE163256 分选数据按患者、组织、试验日及发病后小时合并技术重复，edgeR TMM 标准化。保留原早期 <96 h、同一患者同一试验日、两组织时间差 ≤6 h 的配对。每个患者先在可配对日期内平均差值；18 人参与主要组织对比。原 7 基因 × 2 类细胞 × 2 时段构成 28 项检验家族。
本轮另做 log2(CPM+c) 的 c=0.1、1、5 敏感性、留一患者法、4 类来源标记在血液/血肿/配对差值中的 Spearman 相关（12 项 BH 家族）。留一法保持原可评估假设及 28 项家族。
临床探索使用 19 人各自最早的早期血肿样本，计算一年 mRS、初始血肿体积、采样时间的 3 项相关并作 BH 校正；另以 mRS 秩为结局，PLXNB3、初始血肿体积、年龄的秩为自变量作敏感性回归，报告 PLXNB3 目标效应的 P。
单细胞选用 GSE333172 的 4 份标签一致的血肿样本及 GSE266873 的 1 份 28 h 血肿周围组织样本，按下载可用性在查看表达结果前选定。下载核对 HTTP 区间、文件大小、gzip 完整性并保存 SHA256。
细胞 QC：检测基因 200–6,000、线粒体比例 ≤20%、总 UMI ≤40,000；严格敏感性改为 ≥500 基因、线粒体 ≤10%。根据固定经典标记的检测数及每万 UMI 归一化 log1p 均值作初步大类标注，分类不使用 PLXNB3。额外定义 ≥3 个髓系标记且无胶质/血小板标记组合共现的严格髓系集合。具体标记及阈值保存在 sc_analysis_plan.json。
PLXNB3 一致映射至 ENSG00000198753。单细胞部分只报告每样本阳性细胞、UMI、检测率及伪总体 CPM，不作细胞级推断性检验。

结果一：全血及分选样本中的差异保留，但不能把两种对比称为同一现象的独立复现
全血 72 h 相对基线：26 对，23 人增加；平均 log2(x+1) 差 1.042，中位差 0.816，P=8.72e-6，FDR=0.03214。
全血第 7 天相对基线：24 对，22 人增加；平均差 1.249，中位差 1.371，P=8.81e-7，FDR=0.01253。
第 7 天相对 72 h：P=0.668，FDR=0.968，当前数据不支持此阶段继续增加。
加入血小板及中性粒细胞转录代理量后的全血时间效应仍名义显著：72 h P=0.000544，第 7 天 P=0.001595；两者全家族 FDR=0.3031。效应方向和量级保留，校正后证据减弱。
早期 CD14+ 分选血肿相对外周血：18 人中 15 人增加，中位 TMM logCPM 差 3.182，P=0.003205，原 28 项家族 FDR=0.02991。
留一法原始 P 均 <0.006，但 FDR<0.05 仅 {int((lo.fdr_28<.05).sum())}/18 次；最差 FDR={lo.fdr_28.max():.6f}。因此差异方向有支持，FDR 阈值附近的稳定性有限。
c=0.1/1/5 时，中位差分别 3.591/1.365/0.477，P 分别 0.003205/0.002282/0.001919，均 15/18 人增加。这些是敏感性结果，未重新包装成新的 FDR 阳性发现。方向稳定，但低表达变换明显影响倍数表述。
GSE125512 外部全血矩阵没有 PLXNB3，不能进行同靶点时间变化验证。

结果二：独立单细胞小规模核查将来源线索指向胶质细胞
GSE333172 共 4 位患者、{h['All']['cells']} 个 QC 后细胞，{h['All']['positive']} 个 PLXNB3 阳性细胞，合计 {h['All']['total_UMI']} UMI。
初步髓系类别共 {h['Myeloid']['cells']} 个细胞，PLXNB3 阳性 {h['Myeloid']['positive']} 个；进一步限制来源共现后为 {h['Myeloid_clean']['cells']} 个细胞，阳性 {h['Myeloid_clean']['positive']} 个。
少突胶质细胞类别共 {h['Oligodendrocyte']['cells']} 个细胞，PLXNB3 阳性 {h['Oligodendrocyte']['positive']} 个。
严格 QC 下，髓系类别共 {int(strict.n_cells.sum())} 个细胞，阳性 {int(strict.PLXNB3_positive_cells.sum())} 个。
GSE266873 的单份 28 h 血肿周围组织：{t['All']['cells']} 个 QC 后细胞，PLXNB3 阳性 {t['All']['positive']} 个，合计 {t['All']['total_UMI']} UMI；其初步分类如下：
{positive[positive.gsm=='GSM8255347'].groupby('label').size().to_string()}

其中 2 个髓系标注阳性细胞各只有 1 个 PLXNB3 UMI；两者同时检出 MBP（7/6 UMI）和 PLP1（均 1 UMI），因此在排除胶质标记共现的严格髓系集合中不再保留。这一敏感性筛选不能证明其为污染或双细胞；真实髓系自身表达、吞噬携带的组织 RNA 和环境 RNA 仍无法区分。全部 46 个阳性细胞在严格 QC 下仍保留。4 份血肿样本的 12 个阳性细胞均只有 1 UMI，需强调检出稀疏。\n\n逐样本总览：
{tab.to_string(index=False)}

这些结果显示本轮检出稀疏，支持把胶质来源放在后续核查首位；当前数据不支持已确立的髓系细胞 PLXNB3 表达机制。

结果三：来源关联与临床结局没有形成新的校正后证据
18 人血肿样本中，PLXNB3 与 OPC 标记代理分数相关 rho=0.528，P=0.02417，FDR12=0.28998，属于名义显著线索。星形胶质、成熟少突及血小板代理相关不满足 P<0.05。标记共变不能鉴定细胞来源或证明污染。
19 人一年 mRS：rho=0.07996，P=0.74487，FDR3=0.74487；调整初始血肿体积和年龄的秩回归中 PLXNB3 系数=-0.1507，95% CI [-0.7043,0.4029]，P=0.57043。当前数据不支持预后主线。
初始血肿体积 rho=0.450，P=0.05324；采样时间 rho=-0.399，P=0.09016；均未达预设 P<0.05 筛选标准。

结果四：已知机制只能提供假设，尚未在本数据中建立通路
早期最早血肿样本中，PLXNB3 的患者中位 CPM=0.991，9/19 人达到 CPM≥1；SEMA5A 中位 CPM=0.174，0/19 人达到 CPM≥1。RAC1、ARHGDIA、FSCN1、MET 可测，MST1R 未检出。
现有数据不支持建立 SEMA5A–PLXNB3 自分泌机制。其他细胞提供配体的旁分泌假设本轮未验证。RAC1 的 RNA 不能替代 RAC1-GTP 活性，MET 的 RNA 不能替代磷酸化状态。

结果五：ICH 特异研究有探索空间，但脑损伤和胶质细胞背景已有文献
2026-09-30 Europe PMC 按记录的标题/摘要检索式，PLXNB3/plexin B3/plexin-B3 与 intracerebral hemorrhage/haemorrhage 组合为 0 条；与颅内感染、脑膜炎、脑室炎组合也为 0 条。该结果限定于检索字段和收录范围，不能写成“首次发现”。
2021 Communications Biology 已研究 plexin-B3 阳性成年 OPC 及皮层刺伤/KCl 损伤反应（10.1038/s42003-021-02404-7）。因此“PLXNB3 与脑损伤/少突胶质细胞有关”本身不新。
正常 HPA 图谱显示胶质细胞群富集，并包含巨核细胞增强表达；HPA 免疫细胞未检出。Haematologica 原始血小板制备蛋白组附件也列有 PLXNB3。
SEMA5A–PLXNB3–RhoGDIalpha/RAC1 抑制迁移来自胶质瘤研究；HIF1–PLXNB3–MET/SRC 来自乳腺癌研究，不能直接迁移成脑出血机制。
另检索到 2024 年 AD 微胶质研究预印本提及计算推断的 SEMA3A/PLXNB3 通讯；本轮不把该推断作为已验证配体结合证据。详细来源及边界见 literature_evidence_matrix.csv。

讨论与下一步判断
1. 本轮为明确界定的小规模来源核查，未完成原计划 9 份血肿周围组织及 6 对血肿—血液的全队列验证。NCBI 大文件初轮出现大量 503，故按下载可用性收敛为 5 个样本；全部范围变更、部分文件和下载清单保留。不得将便利样本描述为整个队列代表或独立配对复现。
2. 髓系为经典标记的初步大类，不能区分单核细胞与活化微胶质细胞；未取得作者逐细胞注释。低 RNA 含量、掉落、环境 RNA、双细胞及吞噬后的组织 RNA 均需考虑。无检出不等于不存在表达；胶质来源线索也不等于已证明原分选样本污染。当前质控对低 RNA 血小板不敏感，不能借此排除全血血小板来源。
3. GSE333172 患者 2–5 的 GEO 组织字段与文件名冲突，未自行改标签；原作者 GitHub 公开链接返回 404。文件、样本映射及检索失败均记录。原文说明此队列与 MISTIE 分选验证队列独立招募；本轮文件患者编号未与原分选 21 人重合。
4. 全血时间变化、分选组织差异和单细胞来源是不同问题。GSE266873 本轮无健康对照、无足够生物学重复，不能据此声称 PLXNB3 在 ICH 中上调或推动损伤/修复。
5. 若继续纯公共数据路线，优先完整核查 GSE266873 的少突/OPC 状态及 GSE333172 的标签一致配对样本，以患者为单位分析，再评估是否形成“脑出血后 PLXNB3 的细胞来源及少突胶质谱系状态”主线。只有这一证据链得到支持，才考虑扩大关联通路分析。
6. 当前保留候选，暂停把单核/微胶质机制轴或预后模型作为主线。单基因 FDR 阳性不能单独支撑期刊档次判断；不能基于本轮结果承诺 3 分以上或 8–10 分发表。若独立来源与时间变化仍不能稳定验证，应更换主靶点，PLXNB3 保留为次要结果。

关键参考
HPA: https://www.proteinatlas.org/ENSG00000198753-PLXNB3/single+cell
OPC/脑损伤: https://doi.org/10.1038/s42003-021-02404-7
胶质瘤机制: https://pmc.ncbi.nlm.nih.gov/articles/PMC2952245/
缺氧/乳腺癌机制: https://doi.org/10.1016/j.celrep.2023.112164
独立血肿队列: https://doi.org/10.1177/0271678X261471282
血肿周围组织: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266873

可复核文件
results/: 完整统计表、逐患者输入、单细胞逐样本汇总、阳性细胞审计、标记计数、结果图。
metadata/: 原/修订分析口径、样本清单、下载 SHA256 和范围记录。
scripts/: 核心分析及下载脚本。公开矩阵保留在服务器 raw/singlecell/，压缩交付包不包含这些大文件。
"""
(p/'results/PLXNB3_further_investigation_zh.txt').write_text(text,encoding='utf-8')
assert json.loads((p/'logs/sc_analysis_complete.json').read_text())['ok']
assert json.loads((p/'logs/download_complete.json').read_text())['ok']
files=[]
for folder in ['results','metadata','scripts','logs']:
 files.extend(f for f in (p/folder).rglob('*') if f.is_file())
with zipfile.ZipFile(p/'PLXNB3_investigation_v4.zip','w',zipfile.ZIP_DEFLATED) as z:
 for f in files:z.write(f,f.relative_to(p))
print('PACKAGE', (p/'PLXNB3_investigation_v4.zip').stat().st_size,'REPORT',len(text),'FILES',len(files))
