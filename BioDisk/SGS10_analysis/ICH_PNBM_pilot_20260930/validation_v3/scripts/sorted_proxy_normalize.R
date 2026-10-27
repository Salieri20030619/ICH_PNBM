suppressPackageStartupMessages(library(edgeR))
a<-commandArgs(trailingOnly=TRUE);x<-as.matrix(read.csv(gzfile(a[1]),row.names=1,check.names=FALSE));d<-calcNormFactors(DGEList(x))
gg<-c("ITGA2B","ITGB3","GRAMD1C","DEFA3","ALOX15","PRSS33","PLXNB3","PF4","PPBP","GP1BA","GP9","TUBB1")
write.csv(cpm(d,log=TRUE,prior.count=1)[intersect(gg,rownames(x)),,drop=FALSE],a[2])
