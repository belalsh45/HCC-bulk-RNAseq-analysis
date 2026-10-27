# Gene-symbol-labelled paired VST change heatmap; requires model and prepared SE
suppressPackageStartupMessages({library(DESeq2);library(SummarizedExperiment);library(pheatmap)})
dds<-readRDS('results/HCC_DESeq2_model.rds');se<-readRDS('data/TCGA_LIHC_paired_SE.rds')
r<-read.csv('results/HCC_DESeq2_Tumor_vs_Normal.csv',stringsAsFactors=FALSE)
sig<-subset(r,!is.na(padj)&padj<0.05&abs(log2FoldChange)>1)
up<-head(sig[sig$log2FoldChange>1,][order(sig[sig$log2FoldChange>1,'padj']),],20)
dn<-head(sig[sig$log2FoldChange<(-1),][order(sig[sig$log2FoldChange<(-1),'padj']),],20)
genes<-c(up$gene_id,dn$gene_id)
stopifnot(length(genes)==40,all(genes %in% rownames(dds)))
gdat<-as.data.frame(rowData(se));possible<-c('external_gene_name','gene_name','hgnc_symbol','gene_symbol')
field<-intersect(possible,colnames(gdat));stopifnot(length(field)>0)
strip<-function(x)sub('\\.[0-9]+$','',x)
ix<-match(strip(genes),strip(rownames(gdat)))
lab<-as.character(gdat[[field[1]]][ix]);bad<-is.na(lab)|!nzchar(trimws(lab))
lab[bad]<-strip(genes)[bad];lab<-make.unique(lab)
vsd<-varianceStabilizingTransformation(dds,blind=FALSE)
expr<-assay(vsd)[genes,,drop=FALSE];info<-as.data.frame(colData(dds))
stopifnot(identical(colnames(expr),rownames(info)))
patients<-sort(unique(as.character(info$patient_id)))
delta<-matrix(NA_real_,nrow=length(genes),ncol=length(patients),dimnames=list(lab,patients))
for(i in seq_along(patients)){
 n<-which(as.character(info$patient_id)==patients[i]&as.character(info$condition)=='Normal')
 t<-which(as.character(info$patient_id)==patients[i]&as.character(info$condition)=='Tumor')
 stopifnot(length(n)==1,length(t)==1)
 delta[,i]<-expr[,t]-expr[,n]
}
stopifnot(all(is.finite(delta)))
co<-cor(delta,method='pearson');stopifnot(all(is.finite(co)))
dm<-1-co;dm[dm<0]<-0;diag(dm)<-0
hc<-hclust(as.dist(dm),method='average')
lim<-max(abs(as.numeric(quantile(delta,c(.02,.98),na.rm=TRUE))))
stopifnot(lim>0)
dir.create('figures',showWarnings=FALSE);dir.create('results',showWarnings=FALSE)
write.csv(delta,'results/HCC_Paired_VST_Differences.csv')
for(ext in c('png','pdf')) {
 pheatmap(delta,cluster_rows=FALSE,cluster_cols=hc,gaps_row=20,
 color=colorRampPalette(c('#2166AC','white','#B2182B'))(100),
 breaks=seq(-lim,lim,length.out=101),angle_col=45,fontsize_row=9,fontsize_col=8,
 main='TCGA-LIHC: Paired Tumor-Normal Expression Changes',border_color=NA,
 filename=paste0('figures/HCC_Paired_Delta_Heatmap.',ext),width=11,height=10)
}
cat('Heatmap generated for',nrow(delta),'genes ×',ncol(delta),'pairs\n')
