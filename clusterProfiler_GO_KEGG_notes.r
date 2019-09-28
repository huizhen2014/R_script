##clusterprofiler analysis
library(clusterProfiler)
DE_28vs21_h_keys <- as.character(DE_28vs21_h$GeneNames)
##enrichGO analysis / turn gene SYMBOL to ENTREZID
DE_28vs21_h_keys <- bitr(DE_28vs21_h_keys,fromType = "SYMBOL",
                         toType=c("ENTREZID","REFSEQ"),OrgDb = kp_db)
##Do not recommend using SYMBOL directly, use setReadable function to 
##translate geneID to gene symbol
library(AnnotationDbi)
kp_db <- loadDb("/Users/carlos/.AnnotationHub/74476/org.Klebsiella_pneumoniae_subsp._pneumoniae_HS11286.eg.sqlite")
##or user makeOrgpackage() with AnnotationForge package
ego_mf <- enrichGO(gene=DE_28vs21_h_keys$ENTREZID,OrgDb = kp_db,
                   ont="MF",pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,readable = T)
##setReadable(ego_mf,OrgDb = kp_db)
head(summary(ego_mf))
dotplot(ego_mf,showCategory=30)
cnetplot(ego_mf,foldChange = geneList_28vs21_h)
emapplot(ego_mf)
upsetplot(ego_mf)
plotGOgraph(ego_mf)

##clusterProfiler KEGG module analysis"
kpc <- search_kegg_organism("Klebsiella pneumoniae",by="scientific_name",
                            ignore.case = T)
kkm_h <- enrichMKEGG(gene=as.character(DE_28vs21_h_keys$SYMBOL),organism = "kpm",
                     keyType="kegg",pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",minGSSize = 1,maxGSSize = 500,
                     qvalueCutoff = 0.05)
##clusterProfiler KEGG analysis
## This parameter is by default setting to FALSE,
##and enrichKEGG function will download the latest KEGG data for enrichment analysis. 
##If the parameter use_internal_data is explicitly setting to TRUE, 
##it will use the KEGG.db which is still supported but not recommended.
##since KEGG.db is not updated since 2012. 
kk_h <- enrichKEGG(gene=as.character(DE_28vs21_h_keys$SYMBOL),organism = "kpm",
                   pvalueCutoff = 0.05,pAdjustMethod = "BH",
                   minGSSize = 1, maxGSSize = 500, qvalueCutoff = 0.05,
                   use_internal_data = FALSE)
summary(kk_h)
##clusterProfiler analyze data witt unsupported organisms, GO and KEGG enrichment
##GO via blastgo / KEGG via KAAS
##TERM2GENE, data.frame first column term ID, second cloumn mapped gene
##TERM2NAME, data.frame first column term ID, second column corresponding term name
##universal enrichment analyzer, this example not proper. Usually, it could be used to 
##analyze the phenotype with corresponding gene expression which are annotated by user
term2gene <- select(kp_db,keys=as.character(rownames(cts)),columns=c("GO"),
                    keytype = "SYMBOL")
term2gene <- term2gene[,c(2,1)]
gene_vector <- as.character(DE_28vs21_h_keys$SYMBOL)
gene_universal <- as.character(rownames(cts))
go_enricher <- enricher(gene=gene_vector,pvalueCutoff = 0.05,pAdjustMethod = "BH",
                        universe = gene_universal,minGSSize = 1,maxGSSize = 500,
                        qvalueCutoff = 0.05,TERM2GENE = term2gene)

