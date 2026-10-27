# Hallmark preranked GSEA from DESeq2 Wald statistics
# Requires local CRAN msigdbr + Bioconductor clusterProfiler/enrichplot
suppressPackageStartupMessages({library(msigdbr);library(clusterProfiler);library(enrichplot);library(ggplot2)})
dir.create('results',showWarnings=FALSE);dir.create('figures',showWarnings=FALSE)
r<-read.csv('results/HCC_DESeq2_Tumor_vs_Normal.csv',stringsAsFactors=FALSE)
r<-r[is.finite(r$stat)&!is.na(r$gene_id),]
r$ensembl<-sub('\\.[0-9]+$','',r$gene_id)
r<-r[order(-abs(r$stat)),];r<-r[!duplicated(r$ensembl),]
ranked<-setNames(r$stat,r$ensembl);ranked<-sort(ranked,decreasing=TRUE)
h<-msigdbr(species='Homo sapiens',collection='H')
stopifnot(all(c('gs_name','ensembl_gene') %in% names(h)))
t2g<-unique(data.frame(term=h$gs_name,gene=h$ensembl_gene))
t2g<-t2g[!is.na(t2g$gene)&nzchar(t2g$gene),]
stopifnot(length(unique(t2g$term))==50,length(intersect(names(ranked),t2g$gene))>=1000)
set.seed(2026)
gs<-clusterProfiler::GSEA(geneList=ranked,TERM2GENE=t2g,exponent=1,
 minGSSize=10,maxGSSize=500,pvalueCutoff=1,pAdjustMethod='BH',eps=0,
 by='fgsea',seed=TRUE,verbose=FALSE)
tab<-as.data.frame(gs);stopifnot(nrow(tab)>0)
tab<-tab[order(tab$p.adjust,na.last=TRUE),]
write.csv(tab,'results/Hallmark_GSEA_All_Results.csv',row.names=FALSE)
write.csv(tab[!is.na(tab$p.adjust)&tab$p.adjust<.05,],
 'results/Hallmark_GSEA_Significant.csv',row.names=FALSE)
saveRDS(gs,'results/Hallmark_GSEA_Model.rds')
version<-if('db_version' %in% names(h))paste(unique(h$db_version),collapse=', ') else 'unknown'
qc<-c(paste('R:',R.version.string),paste('MSigDB:',version),
 paste('clusterProfiler:',as.character(packageVersion('clusterProfiler'))),
 paste('msigdbr:',as.character(packageVersion('msigdbr'))),
 paste('Ranked genes:',length(ranked)),paste('Pathways tested:',nrow(tab)),
 paste('Significant FDR<0.05:',sum(tab$p.adjust<.05,na.rm=TRUE)))
writeLines(qc,'results/Hallmark_GSEA_QC.txt')
writeLines(capture.output(sessionInfo()),'results/GSEA_session_info.txt')
show<-tab[!is.na(tab$p.adjust)&tab$p.adjust<.05,]
if(nrow(show)==0)show<-head(tab,15) else show<-head(show,20)
show$Label<-factor(gsub('_',' ',sub('^HALLMARK_','',show$Description)),levels=rev(gsub('_',' ',sub('^HALLMARK_','',show$Description))))
show$Direction<-ifelse(show$NES>0,'Tumor-enriched','Normal-enriched')
p<-ggplot(show,aes(NES,Label,fill=Direction))+geom_col()+geom_vline(xintercept=0)+
 scale_fill_manual(values=c('Tumor-enriched'='#B2182B','Normal-enriched'='#2166AC'))+
 labs(title='Hallmark GSEA in TCGA-LIHC',subtitle='Paired Tumor vs Normal | 10 patients',x='NES',y=NULL)+theme_classic(base_size=12)
ggsave('figures/Hallmark_GSEA_NES_Barplot.png',p,width=11,height=8,dpi=300,bg='white')
ggsave('figures/Hallmark_GSEA_NES_Barplot.pdf',p,width=11,height=8)
# Prefer E2F targets if tested; otherwise best adjusted-P pathway
id<-if('HALLMARK_E2F_TARGETS' %in% tab$ID)'HALLMARK_E2F_TARGETS' else tab$ID[1]
curve<-enrichplot::gseaplot2(gs,geneSetID=id,title=gsub('_',' ',id),pvalue_table=TRUE)
ggsave('figures/Hallmark_GSEA_Enrichment_Curve.png',curve,width=10,height=8,dpi=300,bg='white')
cat('GSEA exported:',nrow(tab),'pathways, FDR-significant',sum(tab$p.adjust<.05,na.rm=TRUE),'\n')
