# Paired DESeq2 and volcano; requires data/TCGA_LIHC_paired_SE.rds
suppressPackageStartupMessages({library(DESeq2);library(SummarizedExperiment);library(ggplot2)})
dir.create('results',showWarnings=FALSE);dir.create('figures',showWarnings=FALSE)
se <- readRDS('data/TCGA_LIHC_paired_SE.rds')
cnt <- assay(se,'unstranded'); b <- colnames(cnt)
code <- substr(b,14,15)
stopifnot(all(code %in% c('01','11')))
info <- data.frame(patient_id=factor(substr(b,1,12)),
 condition=factor(ifelse(code=='01','Tumor','Normal'),levels=c('Normal','Tumor')),row.names=b)
stopifnot(all(rowSums(table(info$patient_id,info$condition)==1)==2))
dds <- DESeqDataSetFromMatrix(countData=cnt,colData=info,design=~patient_id+condition)
dds <- dds[rowSums(counts(dds)>=10)>=5,]
dds <- DESeq(dds)
r <- as.data.frame(results(dds,contrast=c('condition','Tumor','Normal'),alpha=0.05))
r$gene_id <- rownames(r);r <- r[order(r$padj,na.last=TRUE),]
write.csv(r,'results/HCC_DESeq2_Tumor_vs_Normal.csv',row.names=FALSE)
saveRDS(dds,'results/HCC_DESeq2_model.rds')
r$Group <- 'Not Significant'
r$Group[!is.na(r$padj)&r$padj<0.05&r$log2FoldChange>1] <- 'Upregulated'
r$Group[!is.na(r$padj)&r$padj<0.05&r$log2FoldChange<(-1)] <- 'Downregulated'
r$Group <- factor(r$Group,levels=c('Downregulated','Not Significant','Upregulated'))
r$plot_padj <- -log10(pmax(r$padj,.Machine$double.xmin))
p <- ggplot(r,aes(x=log2FoldChange,y=plot_padj,color=Group))+
 geom_point(alpha=.6,size=1.1,na.rm=TRUE)+
 scale_color_manual(values=c(Downregulated='#4575B4','Not Significant'='grey75',Upregulated='#D73027'))+
 geom_vline(xintercept=c(-1,1),linetype='dashed',color='grey40')+
 geom_hline(yintercept=-log10(.05),linetype='dashed',color='grey40')+
 labs(title='TCGA-LIHC: Tumor vs Normal',subtitle='Paired DESeq2 analysis (10 patients)',
 x='log2 Fold Change (Tumor / Normal)',y='-log10 Adjusted P-value',color='Gene status')+
 theme_classic(base_size=13)
ggsave('figures/Volcano_Plot.png',p,width=9,height=7,dpi=300,bg='white')
cat('Tested',nrow(r),'features; threshold hits',sum(r$Group!='Not Significant'),'\n')
