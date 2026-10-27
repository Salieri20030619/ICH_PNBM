suppressPackageStartupMessages(library(limma))
r<-"/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930/validation_v2"
d<-read.delim(gzfile(file.path(r,"raw/external.xls.gz")),check.names=FALSE)
dup<-duplicated(d$symbol)|duplicated(d$symbol,fromLast=TRUE)|is.na(d$symbol)|d$symbol==""
d<-d[!dup,]
x<-as.matrix(d[,paste0("ICH",1:22)]);rownames(x)<-d$symbol
keep<-rowSums(x>=.1)>=ceiling(.2*ncol(x)) & apply(x,1,var)>0
x<-x[keep,,drop=FALSE];y<-log2(x+1);delta<-y[,seq(2,22,2)]-y[,seq(1,21,2)]
fit<-eBayes(lmFit(delta,matrix(1,11,1)),robust=TRUE,trend=FALSE)
wp<-apply(delta,1,function(z) if(all(z==0)) 1 else suppressWarnings(wilcox.test(z,exact=FALSE)$p.value))
out<-data.frame(gene=rownames(y),n_pairs=11,mean_log2_change=rowMeans(delta),median_log2_change=apply(delta,1,median),
fraction_increased=rowMeans(delta>0),p_limma=fit$p.value[,1],p_wilcoxon=wp)
out$fdr_limma<-p.adjust(out$p_limma,"BH");out$fdr_wilcoxon<-p.adjust(out$p_wilcoxon,"BH")
write.csv(out,gzfile(file.path(r,"results/external_gene_tests.csv.gz")),row.names=FALSE)
writeLines(rownames(x),file.path(r,"metadata/external_gene_universe.txt"))
cat("External genes",nrow(out),"limma FDR",sum(out$fdr_limma<.05),"rank FDR",sum(out$fdr_wilcoxon<.05),"\n")
