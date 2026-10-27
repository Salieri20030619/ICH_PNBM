p <- '/mnt/BioDisk/wangmingyuan/SGS10_analysis/ICH_PNBM_pilot_20260930/validation_v5'
x <- readRDS(file.path(p,'raw/published_GSE266873_de_novo_pseudobulk_counts.rds'))
stopifnot(sum(rownames(x)=='PLXNB3')==1)
md <- read.csv(file.path(p,'raw/published_GSE266873_de_novo_pseudobulk_profiles.csv'),check.names=FALSE)
stopifnot(all(md$profile_id %in% colnames(x)))
md$PLXNB3_counts <- as.numeric(x['PLXNB3',md$profile_id])
md$total_UMI <- as.numeric(colSums(x)[md$profile_id])
md$PLXNB3_CPM <- md$PLXNB3_counts/pmax(md$total_UMI,1)*1e6
write.csv(md,file.path(p,'results/extension/published_annotation_plxnb3.csv'),row.names=FALSE)
cat('Published alternative annotation extracted; this is the same biological cohort, not a new cohort.\n')
