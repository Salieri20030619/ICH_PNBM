suppressPackageStartupMessages(library(limma))
b<-"/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930";r<-file.path(b,"validation_v3")
panels<-list(platelet=c("PF4","PPBP","GP1BA","GP1BB","GP9","TUBB1"),neutrophil=c("CSF3R","FCGR3B","CEACAM8","S100A8","S100A9","FPR1"))
x<-as.matrix(read.csv(gzfile(file.path(b,"raw/mrna.csv.gz")),row.names=1,check.names=FALSE));colnames(x)<-sub("_bas$","_base",sub("_7da$","_7day",colnames(x)))
y<-log2(x+1);univ<-readLines(file.path(b,"validation_v2/metadata/ich_gene_universe.txt"))
stopifnot(all(unlist(panels)%in%rownames(y)));proxy<-sapply(panels,function(g)colMeans(y[g,,drop=FALSE]))
write.csv(proxy,file.path(r,"metadata/discovery_proxy_scores.csv"))
out<-list();corout<-list();k<-0
for(tt in list(c("base","72h"),c("72h","7day"),c("base","7day"))){
 k<-k+1;i0<-sub(paste0("_",tt[1],"$"),"",grep(paste0("_",tt[1],"$"),colnames(y),value=TRUE));i1<-sub(paste0("_",tt[2],"$"),"",grep(paste0("_",tt[2],"$"),colnames(y),value=TRUE));ids<-intersect(i0,i1)
 c0<-paste0(ids,"_",tt[1]);c1<-paste0(ids,"_",tt[2]);delta<-y[univ,c1]-y[univ,c0]
 dp<-proxy[c1,,drop=FALSE]-proxy[c0,,drop=FALSE];design<-cbind(time_intercept=1,dp)
 fit<-eBayes(lmFit(delta,design),robust=TRUE,trend=FALSE)
 out[[k]]<-data.frame(gene=univ,contrast=paste0(tt[2],"_vs_",tt[1]),n_pairs=length(ids),adjusted_time_effect=fit$coefficients[,1],p_raw=fit$p.value[,1],proxy_gene=univ%in%unlist(panels),design_condition=kappa(design))
 for(g in c("ITGA2B","ITGB3","GRAMD1C","DEFA3","ALOX15","PRSS33","PLXNB3")){
  fit0<-lm(delta[g,]~dp[,1]+dp[,2])
  corout[[length(corout)+1]]<-data.frame(gene=g,contrast=paste0(tt[2],"_vs_",tt[1]),platelet_rho=cor(delta[g,],dp[,1],method="spearman"),neutrophil_rho=cor(delta[g,],dp[,2],method="spearman"),proxy_R2=summary(fit0)$r.squared)
 }
}
z<-do.call(rbind,out);z$fdr_joint<-p.adjust(z$p_raw,"BH");write.csv(z,gzfile(file.path(r,"results/discovery_proxy_adjusted_all.csv.gz")),row.names=FALSE)
write.csv(do.call(rbind,corout),file.path(r,"results/target_proxy_correlations.csv"),row.names=FALSE)
d<-read.delim(gzfile(file.path(b,"validation_v2/raw/external.xls.gz")),check.names=FALSE);dup<-duplicated(d$symbol)|duplicated(d$symbol,fromLast=TRUE)|is.na(d$symbol)|d$symbol=="";d<-d[!dup,]
xx<-as.matrix(d[,paste0("ICH",1:22)]);rownames(xx)<-d$symbol;yy<-log2(xx+1);uu<-readLines(file.path(b,"validation_v2/metadata/external_gene_universe.txt"))
pp<-sapply(panels,function(g)colMeans(yy[g,,drop=FALSE]));dd<-yy[uu,seq(2,22,2)]-yy[uu,seq(1,21,2)];dp<-pp[seq(2,22,2),]-pp[seq(1,21,2),]
fit<-eBayes(lmFit(dd,cbind(time_intercept=1,dp)),robust=TRUE,trend=FALSE)
zz<-data.frame(gene=uu,n_pairs=11,adjusted_time_effect=fit$coefficients[,1],p_raw=fit$p.value[,1],proxy_gene=uu%in%unlist(panels));zz$fdr<-p.adjust(zz$p_raw,"BH")
write.csv(zz,gzfile(file.path(r,"results/external_proxy_adjusted_all.csv.gz")),row.names=FALSE)
writeLines("0",file.path(r,"logs/proxy.exit"))
cat("Discovery adjusted FDR:",sum(z$fdr_joint<.05 & !z$proxy_gene)," External:",sum(zz$fdr<.05 & !zz$proxy_gene),"\n")
