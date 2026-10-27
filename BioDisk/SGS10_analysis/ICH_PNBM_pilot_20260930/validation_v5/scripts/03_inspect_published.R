p <- '/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930/validation_v5'
fs <- list.files(file.path(p,'raw'),pattern='published_.*counts.rds$',full.names=TRUE)
for(f in fs){x<-readRDS(f);str(x,max.level=1); if(is.matrix(x)||inherits(x,'Matrix')){cat('PLXNB3_ROW',which(rownames(x)=='PLXNB3'),'\n'); print(head(colnames(x)));}}
