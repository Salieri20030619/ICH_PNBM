from pathlib import Path
import pandas as pd,numpy as np,json,zipfile,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930');P=B/'validation_v6';E=B/'validation_v6_external';V=B/'validation_v5'
tests=pd.read_csv(P/'results/matched_state_tests.csv');ef=pd.read_csv(P/'results/matched_patient_effects.csv');ex=pd.read_csv(E/'results/matched_state_tests.csv')
marker=pd.read_csv(P/'results/leave_one_marker_out.csv');loo=pd.read_csv(P/'results/leave_one_patient_out.csv');genes=pd.read_csv(P/'results/all_marker_tests.csv')
# Summarize key claims with exact statistics, including sensitivity failures.
claims=[]
for name in ['Differentiating_OL','Mature_myelin']:
 a=tests[(tests.setting=='standard')&(tests.program==name)].iloc[0].to_dict()
 a['strict_qc_fdr']=float(tests[(tests.setting=='strict_QC')&(tests.program==name)].fdr_6.iloc[0])
 a['lineage_specific_fdr']=float(tests[(tests.setting=='lineage_specific')&(tests.program==name)].fdr_6.iloc[0])
 a['worst_leave_patient_fdr']=float(loo[loo.program==name].fdr_6.max())
 a['worst_leave_marker_fdr']=float(marker[marker.program==name].fdr_6.max())
 a['external_n']=int(ex[(ex.setting=='standard')&(ex.program==name)].n.iloc[0])
 a['external_n_up']=int(ex[(ex.setting=='standard')&(ex.program==name)].n_up.iloc[0])
 claims.append(a)
pd.DataFrame(claims).to_csv(P/'results/key_evidence_summary.csv',index=False)
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(1,3,figsize=(15,4.8),gridspec_kw={'width_ratios':[1.2,1,1.3]})
colors={'Differentiating_OL':'#487fb1','Mature_myelin':'#bd674c'}
for j,name in enumerate(colors):
 y=ef[(ef.setting=='standard')&(ef.program==name)].delta.to_numpy()
 ax[0].scatter(np.linspace(j-.08,j+.08,len(y)),y,color=colors[name],s=35)
 ax[0].hlines(np.median(y),j-.18,j+.18,color='black',lw=2)
 ee=pd.read_csv(E/'results/matched_patient_effects.csv');yy=ee[(ee.setting=='standard')&(ee.program==name)].delta.to_numpy()
 ax[1].scatter(np.linspace(j-.05,j+.05,len(yy)),yy,color=colors[name],s=42)
ax[0].axhline(0,color='gray',ls='--',lw=1);ax[0].set_xticks([0,1],['Differentiation','Mature myelin']);ax[0].set_ylabel('Matched score difference: PLXNB3+ minus PLXNB3-')
ax[0].set_title('Discovery: 9 patients\nBoth: P = 0.00391; FDR6 = 0.01172')
ax[1].axhline(0,color='gray',ls='--',lw=1);ax[1].set_xticks([0,1],['Differentiation','Mature myelin']);ax[1].set_title('Independent hematoma: 3 eligible patients\nBoth: P = 0.25; FDR6 = 0.50');ax[1].set_ylabel('Matched score difference')
for j,name in enumerate(colors):
 dd=tests[tests.program==name].set_index('setting').loc[['standard','strict_QC','lineage_specific']]
 ax[2].plot([0,1,2],dd.fdr_6,marker='o',label=name.replace('_',' '),color=colors[name])
ax[2].axhline(.05,color='gray',ls='--',lw=1);ax[2].set_xticks([0,1,2],['Main match','Strict QC','Lineage filter']);ax[2].set_ylabel('FDR within six programs');ax[2].set_title('Sensitivity matters');ax[2].legend(frameon=False,fontsize=8)
fig.text(.015,.012,'Each dot in the first two panels represents a patient. Small-N exploratory evidence; state association does not establish repair or causality.',fontsize=8)
fig.tight_layout(rect=(0,.07,1,1));fig.savefig(P/'results/PLXNB3_reviewed_state_evidence.png',dpi=180);plt.close(fig)
report="""PLXNB3 扩展研究：显式复核与方向判断
2026-09-30 | 版本 v6

结论
当前数据支持继续保留 PLXNB3，并将主线收敛为“脑出血患者中 PLXNB3 相关的少突胶质细胞髓鞘状态”。本轮完成了全队列结果解读、独立注释核对、双变量匹配、严格 QC/谱系限定、留一患者、逐标记删除和独立血肿队列的同规则验证。
与此前相比，证据新增在于：控制测序深度和检测基因数后，分化及成熟髓鞘标记分数仍在发现队列 9 人中全部同向升高，且留一患者后维持 FDR 显著。
推荐优先保留“成熟髓鞘相关状态”，将“分化状态”作为支持性线索。当前不能写“PLXNB3 驱动髓鞘修复”，也未证明其为脑出血特异性变化。

方法
数据：GSE266873 的 9 位患者；GSE333172 的 6 位标签一致患者、12 份血肿/血液样本；GSE24265 的 4 位患者、11 份配对脑组织芯片样本。
沿用 v5 固定的经典标记大类注释及标准 QC；PLXNB3 不参与分类。v6 为对已有状态信号的敏感性复核，不是独立患者重复。
在每位患者的少突胶质细胞中，按标准化后的 log1p(总 UMI−PLXNB3 UMI) 和 log1p(检测基因数−PLXNB3 是否检出) 进行二维近邻匹配。每个阳性细胞匹配 5 个阴性细胞，允许阴性细胞复用，最远距离不超过 0.5；患者至少 5 个匹配阳性、20 个候选阴性。负对照复用次数、独立控制细胞数和匹配诊断均保存。
每个程序的基因表达为 log1p(每万 UMI 归一化计数)，逐基因计算阳性与匹配阴性的均值差，再对程序基因求均值；每位患者提供一个效应。Wilcoxon 检验以患者为单位，6 个固定状态程序组成 BH 家族。
严格 QC 和严格谱系分别为敏感性分析；后者在原少突胶质类别内要求 SOX10 检出，或 MOG/MAG/MOBP 至少 2 个检出。逐标记删除固定原匹配细胞，只修改对应程序，然后重算完整 6 项检验家族。留一患者法也保持 6 项家族。
另外对全部 37 个程序标记做同类患者内对比，以完整 37 项为 BH 家族，不只挑选阳性程序内基因。
独立血肿队列使用相同匹配和可评估标准。所有样本量以患者为单位，n<10 按探索性证据呈现。
公开 Zenodo 冻结伪总体来自同一 GSE266873 的另一套聚类/注释方法，仅用于分析方法敏感性，未计作新的生物学重复。

结果一：两个状态信号在精细匹配后仍有 FDR 证据
主匹配保留每位患者全部合格阳性细胞；匹配后的 UMI 和检测基因数绝对标准化差异均小于 0.004。
分化程序：9/9 人阳性细胞更高，中位程序分数差 +0.02989，P=0.003906，FDR6=0.011719。
成熟髓鞘程序：9/9 人更高，中位差 +0.08745，P=0.003906，FDR6=0.011719。
分数差是平均 log1p(每万 UMI 归一化表达) 的差值，不能解释成表达倍数、修复幅度或临床效应。
严格 QC 后：分化 9/9 同向，成熟髓鞘 8/9 同向，两者 FDR6=0.023438。
留一患者后：两者全部 9 次 FDR6=0.023438，未出现由单个患者支撑显著性的情况。
OPC 程序在精细匹配后 P=0.05469、FDR=0.08203，不再达到原始 P<0.05 标准。应纠正仅依据上一版分箱匹配而强调 OPC 富集的叙事。
应激、干扰素及吞噬溶酶体整体程序均不满足主匹配 P<0.05；不围绕这些阴性程序拼接机制。

结果二：标记依赖和谱系限定划定了证据边界
严格谱系限定后：分化 8/9、成熟髓鞘 7/9 同向；两者 P=0.019531、FDR6=0.058594，属于名义显著。
分化程序去掉 MYRF 后 P=0.128906、FDR=0.193359；去掉 BCAS1 后 P=0.039063、FDR=0.082031。分化叙事较依赖 MYRF/BCAS1，不能称为广泛、稳定的再髓鞘化程序。
成熟髓鞘程序逐个去掉 MBP、PLP1、MOG、MAG 或 MOBP 后仍保留 FDR<0.05；去掉 OPALIN 后 P=0.019531、FDR=0.058594，仍有方向及名义证据，但不再校正后显著。
37 个标记完整家族中，MYRF、MBP、MOBP、OPALIN 较高的 FDR 均为 0.048177；OAS1、TREM2 较低也通过该阈值。后两者为少突胶质类别内稀疏标记差异，不能外推成髓系抗炎或吞噬机制。所有 37 项结果保留。

结果三：独立血肿队列支持方向，但尚未完成有统计效能的独立验证
6 位患者中，只有 3 位满足至少 5 个 PLXNB3 阳性少突胶质细胞及匹配要求。
主匹配中，分化及成熟髓鞘分数均 3/3 同向，中位差分别 +0.04667 和 +0.46060；两者 P=0.25、FDR6=0.50。
n=3 的双侧精确 Wilcoxon 最小可达 P 为 0.25，因此本结果用于方向一致性，不能因未达 FDR 否定方向，也不能称为独立显著复现。
严格 QC 下成熟髓鞘仍 3/3 同向，分化为 2/3。独立证据更偏向成熟髓鞘相关状态。

结果四：细胞来源已比上一轮清楚，时间和预后主线仍缺证据
9 份血肿周围样本共 40,584 个 QC 后细胞，PLXNB3 阳性 1,496 个；其中 1,411 个初步标为少突胶质细胞，髓系类别 13 个。
患者级少突胶质相对髓系的 PLXNB3 表达 9/9 更高：中位 log2(CPM+1) 差 2.484，P=0.003906、来源家族 FDR2=0.007813。与星形胶质比较仅 6 位患者可评估，P=0.0625，不能称少突特异。
公开替代聚类/注释下少突与合并髓系比较同样 9/9 更高，中位差 2.430、P=0.003906；这是同一数据的敏感性复核。
独立 6 份血肿中共 15,642 个 QC 后细胞，PLXNB3 阳性 49 个，47 个位于初步少突胶质类别；2,148 个髓系细胞未检出。6 份血液仅 10 个阳性细胞，检出过低，不补造配对推断 P。
发现队列中患者级 PLXNB3 少突伪总体 CPM 对采样时间的 rho=-0.485、精确 P=0.18955；6 个患者间状态分数相关也均不满足 P<0.05。患者间总量关系与患者内细胞状态差异是不同问题，不能混写成时间演变。
既有 19 人一年 mRS 关联 rho=0.080、P=0.745 的阴性结论不变，当前不支持预后主线。

结果五：脑组织参照选择改变方向
GSE24265 中，血肿周围组织相对对侧白质为 3/3 更低，中位 log2 RMA 差 -0.59968、P=0.25。
相对对侧灰质则为 4/4 更高，中位差 +0.37798、P=0.125。两者 FDR2=0.25。
这是具有方向一致性的小样本描述证据，提示组织构成和参照区域十分重要；不能合并灰白质后宣称脑出血统一上调。所有样本均来自脑出血患者，对侧组织不是健康人群对照。

讨论与研究安排
1. 当前最合适的表述为“PLXNB3 阳性少突胶质细胞与成熟髓鞘标记状态相关”。成熟髓鞘程序的逐标记稳定性和外部方向一致性优于分化程序。
2. MYRF 是优先跟进的共表达线索；这里没有调控方向、直接结合、干预或蛋白活性证据，不能命名为已经建立的 PLXNB3–MYRF 调控轴。
3. 经典细胞标记、RNA 深度、稀疏检出、环境 RNA、双细胞、吞噬携带的 RNA 仍可能影响标注。严格谱系筛选的敏感性结果已保留，不以筛选后零检出替代原始阳性例外。
4. 全部关键状态发现来自 9 位患者；大量细胞不会增加独立样本量。匹配后的患者效应使用 9 个生物学重复，阴性细胞复用不构成额外重复。分化/成熟状态程序是人工整理标记，不能替代经过验证的功能读出。
5. 本数据没有健康、同脑区、同流程的匹配参考，也缺少干预。需优先区分“正常少突成熟状态”与“脑出血特异状态”；当前不能把关联升级为损伤后新生少突细胞或髓鞘修复。
6. 下一阶段以更完整细胞亚群注释及外部疾病/参考背景核查为优先，维持患者级比较，保留 MYRF、OPALIN 作为候选节点；暂不扩展预后模型、细胞通讯因果链或药物推荐。
7. 结论是“值得继续、主线已收敛”，不是“8–10 分文章已经具备条件”。当前也不能基于单个候选及小样本 FDR 承诺发表。

相关文献边界
既有成年 OPC/脑损伤研究：10.1038/s42003-021-02404-7。
同数据的公开分析复核包：10.5281/zenodo.22080648；其中未见 PLXNB3 文本命中，但这不是首次发现证明。
本轮 MYRF/PLXNB3 定向检索主要出现共同细胞标记资料，未据此认定直接调控。完整原始方法、数据和此前新颖性审计保留在 v4/v5。

交付
key_evidence_summary.csv：关键效应及敏感性边界。
matched_state_tests.csv、matching_diagnostics.csv：全部程序统计与匹配质量。
leave_one_patient_out.csv、leave_one_marker_out.csv、all_marker_tests.csv：全部敏感性及37标记统计。
external/：独立队列同规则验证。
PLXNB3_reviewed_state_evidence.png：患者级效应及 FDR 敏感性图。
"""
(P/'results/PLXNB3_reviewed_report_zh.txt').write_text(report,encoding='utf-8')
with zipfile.ZipFile(P/'PLXNB3_reviewed_v6.zip','w',zipfile.ZIP_DEFLATED) as z:
 for folder in ['results','scripts','metadata','logs']:
  for f in (P/folder).rglob('*'):
   if f.is_file():z.write(f,f.relative_to(P))
 for f in E.rglob('*'):
  if f.is_file():z.write(f,Path('external')/f.relative_to(E))
with zipfile.ZipFile(P/'PLXNB3_reviewed_v6.zip') as z:assert z.testzip() is None
print('REPORT',len(report),'ZIP', (P/'PLXNB3_reviewed_v6.zip').stat().st_size)
print('SHA256',hashlib.sha256((P/'PLXNB3_reviewed_v6.zip').read_bytes()).hexdigest())
