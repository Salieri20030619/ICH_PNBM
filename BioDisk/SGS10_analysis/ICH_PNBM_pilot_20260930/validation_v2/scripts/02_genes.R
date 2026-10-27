suppressPackageStartupMessages(library(limma))
base <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930"
r <- file.path(base,"validation_v2")
x <- as.matrix(read.csv(gzfile(file.path(base,"raw/mrna.csv.gz")),row.names=1,check.names=FALSE))
colnames(x) <- sub("_bas$","_base",sub("_7da$","_7day",colnames(x)))
keep <- rowSums(x>=1)>=ceiling(.2*ncol(x)) & apply(x,1,var)>0
x <- x[keep,,drop=FALSE]; y<-log2(x+1)
cons<-list(c("base","72h"),c("72h","7day"),c("base","7day"))
out<-list()
for(k in seq_along(cons)){
 tt<-cons[[k]]
 ids0<-sub(paste0("_",tt[1],"$"),"",grep(paste0("_",tt[1],"$"),colnames(y),value=TRUE))
 ids1<-sub(paste0("_",tt[2],"$"),"",grep(paste0("_",tt[2],"$"),colnames(y),value=TRUE))
 ids<-intersect(ids0,ids1)
 delta<-y[,paste0(ids,"_",tt[2]),drop=FALSE]-y[,paste0(ids,"_",tt[1]),drop=FALSE]
 fit<-eBayes(lmFit(delta,matrix(1,length(ids),1)),robust=TRUE,trend=FALSE)
 wp<-apply(delta,1,function(z) if(all(z==0)) 1 else suppressWarnings(wilcox.test(z,exact=FALSE)$p.value))
 out[[k]]<-data.frame(gene=rownames(y),contrast=paste0(tt[2],"_vs_",tt[1]),n_pairs=length(ids),
 mean_log2_change=rowMeans(delta),median_log2_change=apply(delta,1,median),fraction_increased=rowMeans(delta>0),
 p_limma=fit$p.value[,1],p_wilcoxon=wp)
}
z<-do.call(rbind,out);z$fdr_limma_joint<-p.adjust(z$p_limma,"BH");z$fdr_wilcoxon_joint<-p.adjust(z$p_wilcoxon,"BH")
write.csv(z,gzfile(file.path(r,"results/ich_all_gene_tests.csv.gz")),row.names=FALSE)
writeLines(rownames(y),file.path(r,"metadata/ich_gene_universe.txt"))
writeLines(capture.output(sessionInfo()),file.path(r,"metadata/R_session.txt"))
print(as.data.frame(aggregate(cbind(limma=z$fdr_limma_joint<.05,rank=z$fdr_wilcoxon_joint<.05),list(contrast=z$contrast),sum)))
writeLines("0",file.path(r,"logs/ich_genes.exit"))
