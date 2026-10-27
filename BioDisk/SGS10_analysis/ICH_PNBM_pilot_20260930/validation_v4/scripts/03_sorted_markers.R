suppressPackageStartupMessages(library(edgeR))
base <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930"
x <- as.matrix(read.csv(file.path(base,"validation_v3/raw/monos_collapsed_counts.csv.gz"),row.names=1,check.names=FALSE))
d <- calcNormFactors(DGEList(x))
genes <- c("PLXNB3","SEMA5A","ARHGDIA","RAC1","FSCN1","MET","MST1R","GFAP","AQP4","SLC1A3","ALDH1L1","HEPACAM","PLP1","MBP","MOG","SOX10","OLIG1","OLIG2","PDGFRA","CSPG4","GPR17","PF4","PPBP","GP1BA","GP9","TUBB1","LST1","TYROBP","FCER1G","CTSS","AIF1","C1QA","CD14","CSF1R")
genes<-intersect(genes,rownames(x))
write.csv(t(cpm(d,log=FALSE)[genes,,drop=FALSE]),file.path(base,"validation_v4/results/sorted_marker_cpm.csv"))
write.csv(t(x[genes,,drop=FALSE]),file.path(base,"validation_v4/results/sorted_marker_counts.csv"))
write.csv(t(cpm(d,log=TRUE,prior.count=1)[genes,,drop=FALSE]),file.path(base,"validation_v4/results/sorted_marker_logcpm.csv"))
