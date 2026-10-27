suppressPackageStartupMessages(library(limma))
base<-"/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930";r<-file.path(base,"validation_v2")
x<-as.matrix(read.csv(gzfile(file.path(base,"raw/mrna.csv.gz")),row.names=1,check.names=FALSE))
colnames(x)<-sub("_bas$","_base",sub("_7da$","_7day",colnames(x)));x<-x[readLines(file.path(r,"metadata/ich_gene_universe.txt")),,drop=FALSE];y<-log2(x+1)
cons<-list(c("base","72h"),c("72h","7day"),c("base","7day"));Ds<-list();keys<-c()
for(k in 1:3){
 tt<-cons[[k]];i0<-sub(paste0("_",tt[1],"$"),"",grep(paste0("_",tt[1],"$"),colnames(y),value=TRUE));i1<-sub(paste0("_",tt[2],"$"),"",grep(paste0("_",tt[2],"$"),colnames(y),value=TRUE));ids<-intersect(i0,i1)
 dd<-y[,paste0(ids,"_",tt[2])]-y[,paste0(ids,"_",tt[1])];colnames(dd)<-ids;Ds[[k]]<-dd;keys<-c(keys,paste(rownames(y),paste0(tt[2],"_vs_",tt[1]),sep="::"))
}
z<-read.csv(gzfile(file.path(r,"results/ich_all_gene_tests.csv.gz")))
sel<-which(z$fdr_limma_joint<.05);selected_keys<-paste(z$gene[sel],z$contrast[sel],sep="::")
at<-match(selected_keys,keys);pids<-unique(sub("_.*","",colnames(y)));qq<-matrix(NA,length(sel),length(pids));effects<-qq
for(j in seq_along(pids)){
 ps<-c();ms<-c()
 for(dd in Ds){
  dd<-dd[,colnames(dd)!=pids[j],drop=FALSE];ff<-eBayes(lmFit(dd,matrix(1,ncol(dd),1)),robust=TRUE,trend=FALSE)
  ps<-c(ps,ff$p.value[,1]);ms<-c(ms,rowMeans(dd))
 }
 qq[,j]<-p.adjust(ps,"BH")[at];effects[,j]<-ms[at]
}
out<-z[sel,];out$loo_fdr_min<-apply(qq,1,min);out$loo_fdr_max<-apply(qq,1,max);out$loo_fraction_fdr_significant<-rowMeans(qq<.05)
out$loo_fraction_same_direction<-rowMeans(sign(effects)==sign(out$mean_log2_change))
write.csv(out,file.path(r,"results/gene_leave_one_patient_out.csv"),row.names=FALSE)
cat("LOO done",nrow(out),"candidates;",sum(out$loo_fraction_fdr_significant==1),"all deletions FDR\n")
writeLines("0",file.path(r,"logs/loo.exit"))
