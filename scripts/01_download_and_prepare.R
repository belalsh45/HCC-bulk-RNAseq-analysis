# Reproducible subset query: GDC TCGA-LIHC, 10 matched pairs
suppressPackageStartupMessages({library(TCGAbiolinks);library(SummarizedExperiment)})
dir.create('data',showWarnings=FALSE);dir.create('results',showWarnings=FALSE)
q <- GDCquery(project='TCGA-LIHC',data.category='Transcriptome Profiling',
 data.type='Gene Expression Quantification',workflow.type='STAR - Counts',
 sample.type=c('Primary Tumor','Solid Tissue Normal'),access='open')
m <- getResults(q)
stopifnot(all(c('cases','sample_type') %in% names(m)))
m$patient_id <- substr(as.character(m$cases),1,12)
stopifnot(all(grepl('^TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}$',m$patient_id)))
pair_counts <- table(m$patient_id,m$sample_type)
eligible <- rownames(pair_counts)[pair_counts[,'Primary Tumor']==1 & pair_counts[,'Solid Tissue Normal']==1]
stopifnot(length(eligible)>=10)
ids <- head(sort(eligible),10)
sel <- m[m$patient_id %in% ids,,drop=FALSE]
stopifnot(nrow(sel)==20)
write.csv(sel,'results/selected_paired_samples.csv',row.names=FALSE)
q2 <- GDCquery(project='TCGA-LIHC',data.category='Transcriptome Profiling',
 data.type='Gene Expression Quantification',workflow.type='STAR - Counts',
 barcode=as.character(sel$cases),sample.type=c('Primary Tumor','Solid Tissue Normal'),access='open')
stopifnot(nrow(getResults(q2))==20)
GDCdownload(q2,method='api',directory='data/GDCdata',files.per.chunk=5)
se <- GDCprepare(q2,directory='data/GDCdata',summarizedExperiment=TRUE)
stopifnot(ncol(se)==20,'unstranded' %in% assayNames(se))
barcodes <- colnames(se); patients <- substr(barcodes,1,12)
condition <- ifelse(substr(barcodes,14,15)=='01','Tumor',ifelse(substr(barcodes,14,15)=='11','Normal',NA_character_))
stopifnot(!anyNA(condition),all(rowSums(table(patients,condition)==1)==2))
cnt <- assay(se,'unstranded')
stopifnot(!anyNA(cnt),all(cnt>=0),all(cnt==round(cnt)))
saveRDS(se,'data/TCGA_LIHC_paired_SE.rds')
writeLines(capture.output(sessionInfo()),'results/session_info.txt')
cat('Prepared',nrow(se),'features and',ncol(se),'samples from',length(unique(patients)),'pairs\n')
