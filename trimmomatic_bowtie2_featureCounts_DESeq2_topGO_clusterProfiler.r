###part 1
##trimmomatic
#trimmomatic PE -summary kp_21_trimmed.log kp_21_R1.fastq.gz kp_21_R2.fastq.gz -baseout kp_21_trimmed.fastq.gz \
#ILLUMINACLIP:Sabgon_adapter.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:35 
##sangon_adapter.fa
##>PrefixPE/1
##AGATCGGAAGAGCACACGTCTGAAC
##>PrefixPE/2
##AGATCGGAAGAGCGTCGTGTAGGGA

###align bowtie2 with shell scripts
#bowtie2-build GCF_000240185.1_ASM24018v2_genomic.fna hs11286_index
#bowtie2 -x hs11286_index -1 kp_21_sangon_1P.fastq.gz -2 kp_21_sangon_2P.fastq.gz -S kp_21_sangon_mapped.sam
#bowtie2 -x hs11286_index -1 kp_28_sangon_1P.fastq.gz -2 kp_28_sangon_2P.fastq.gz -S kp_28_sangon_mapped.sam
#samtools view -bS kp_21_sangon_mapped.sam > kp_21_sangon_mapped.bam
#samtools view -bS kp_28_sangon_mapped.sam > kp_28_sangon_mapped.bam
#samtools sort -n kp_21_sangon_mapped.bam -o kp_21_sangon_mapped_sorted.bam
#samtools sort -n kp_28_sangon_mapped.bam -o kp_28_sangon_mapped_sorted.bam

##part 2
##考虑到基因重组，仅计算primary比对即可
####featureCounts ignoreDup=F,primaryOnly=TRUE
##featureCounts
library(Rsubread)
samples <- c("ab_1_c507_0440","ab_2_c507_0440","ab_3_c507_0440",
             "ab_1", "ab_2")
gtf_file <- "ab030_gff.gtf"
featuretype <- "transcript"
attrtype <- "locus_tag"
anno_file <- "ab030_annotation_extraction.txt"
go_file <- "ab030_Sangon_go.txt"
output <- "ab_c507_0440_vs_ab_c"

Results <- list()
for(sample in samples){
  tmp <- c()
  tmp <- paste0(sample,"_sangon_mapped_sorted.bam")
  Results[[sample]] <- featureCounts(tmp,annot.ext = gtf_file,isGTFAnnotationFile = TRUE,
                                     GTF.featureType = featuretype,
                                     GTF.attrType = attrtype,isPairedEnd = TRUE,
                                     allowMultiOverlap=TRUE,primaryOnly=TRUE,
                                     strandSpecific=2,nthreads=4)
}

###collect counts
fc_counts <- list()
for(sample in samples){
  tmp <- data.frame()
  fc <- list()
  fc <- Results[[sample]]
  tmp <- data.frame(Gene_ID=rownames(fc$counts),Pos=paste0(
    fc$annotation$Chr,"[",fc$annotation$Strand,"]",
    fc$annotation$Start,"-",fc$annotation$End),
    Length=fc$annotation$Length,Count=fc$counts[,1],
    TPM=((fc$counts[,1]/fc$annotation$Length)*1000000)/
      (sum(fc$counts[,1]/fc$annotation$Length))
  )
  fc_counts[[sample]] <- tmp
}

##构建DESeq2 countData/colData
countData <- sapply(fc_counts,function(x)x$Count)
rownames(countData) <- fc_counts[[1]]$Gene_ID
#idx <- apply(countData,1,function(x)length(x[x>0])>1)
#countData <- countData[idx,]
colData <- matrix(666,nrow=length(fc_counts))
rownames(colData) <- colnames(countData)
colData[,1] <- c(rep("ab_c507_0440",3),rep("ab_c",2)) ##根据实际修改
colnames(colData) <- "condition"

##DESeq2差异分析
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,
                              design = ~ condition)

##设置factor levels
dds$condition <- factor(dds$condition,levels=c("ab_c507_0440","ab_c"))
##dds$condition <- relevel(dds$condition,ref="MDR")
##dds$condition <- droplevels(dds$level)

##DE analysis
##1, estimate of size factors: estimateSizeFactors
##2, estimate of dispersion: esitmatedispersions
##3, Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
##4, results函数生成log2倍数改变及对应p值
dds <- DESeq(dds)

##默认为last level vs. ref level
##resultsNames(dds) 查看coefficient名称可知
##这里通过contrast指定 MDR/AS，指定adjusted p-value cutoff (FDR)阈值为0.05
res <- results(dds,contrast=c("condition","ab_c507_0440","ab_c"))
#res <- results(dds,contrast = c("condition","28s","21s"),alpha=0.05)
##removeResults函数返回至DESeqDataSet对象

##collection and print to file
##因为重名，所有.1后缀为对应的TPM值,.2后缀为sizeFactor矫正后counts，
##然后通过sub替换对应名称
library(xlsx)
res_data_frame <- NULL
Total_counts <- NULL
res_data_frame <- as.data.frame(res)
Anno <- read.delim(anno_file,sep="\t",header=FALSE,
                   row.names = 3)
Total_counts <- data.frame(row.names=rownames(fc_counts[[1]]),
                           Gene_ID=fc_counts[[1]]$Gene_ID,
                           Pos=fc_counts[[1]]$Pos,
                           Length=fc_counts[[1]]$Length,
                           sapply(fc_counts,function(x)x$Count),
                           counts(dds,normalized=TRUE),
                           sapply(fc_counts,function(x)x$TPM)
)
colnames(Total_counts) <- sub("(*)\\.1","\\1_Nor_count",colnames(Total_counts))
colnames(Total_counts) <- sub("(*)\\.2","\\1_TPM",colnames(Total_counts))
Total_counts <- cbind(Total_counts,res_data_frame[rownames(Total_counts),])
Total_counts$Anno_location <- Anno[rownames(Total_counts),]$V2
Total_counts$Anno_geneid <- Anno[rownames(Total_counts),]$V4
Total_counts$Anno_protein <- Anno[rownames(Total_counts),]$V5
Total_counts$Anno_product <- Anno[rownames(Total_counts),]$V6
write.table(Total_counts,file=paste0("ab_c507_0440_vs_ab_c","_DE_total_3vs2.txt"),sep="\t",quote=FALSE)
write.xlsx(Total_counts,file =paste0("ab_c507_0440_vs_ab_c","_DE_total_3vs2.xlsx"),row.names=FALSE)

##res结果根据padj排序
resOrdered <- res[order(res$padj),]
##summary(resOrdered)

##lfcshrink, 仅为绘图和排序使用，最终差异结果和results无异
res_shrunken_normal <- lfcShrink(dds,contrast = c("condition","ab_c507_0440","ab_c"),
res=res,type="normal",alpha=0.05)

##绘制MA-plot
#plotMA(res,main="DESeq2",ylim=c(-2,2))
plotMA(res_shrunken_normal,main="DESeq2",ylim=c(-2,2))
##idx <- identity(res$baseMean,res$log2FoldChange)
##rownames(res)[idex] 交互式获得对应点坐标信息

##count plots,检查单个基因的count read
plotCounts(dds,gene=which.min(res$padj),intgroup = "condition")

##查看rsults返回结果描述
##p-values == NA
##1，一行中，所有样本counts都为0,baseMean为0，其他都为NA
##2，若一行中包含极端count值，prl和adjusted p都为NA
##outlier counts的检出是根据Cook's distance
##3，若一行由于低于mean normalized count
##而被automatic independent filtering, 那么只有adjusted p为NA
##mcols(resOrdered)$description
##详细查看metadata(res)
##sizeFacots(dds)
##coef(dds)
##

##为检测差异性表达，对原始counts数据使用离散分布分析。但是为了
##下游其他类型分析，例如可视化和聚类，counts数据对转换会更有用
##rlog/varianceStabilizingformation，都有参数blind，默认为TRUE,
##当期待许多或大部分基因会根据实验设计而出现大的counts差异时，
##此时blind dispersion estimation会导致偏大的离散度评估。通过设置
##blind=FALSE，将已得到dispesion用于数据转弯，若不存在，函数将根据
##design fromula重新评估
rld <- rlog(dds,blind=FALSE)
##绘制数据转换和变异关系图(row stand deviations vs row means, 
##both are transfered)
library(vsn)
notAllZero <- (rowSums(counts(dds))>0)
meanSdPlot(assay(rld[notAllZero,]))
dev.off()

##Heatmap ,使用颜色展示矩阵数据的个体数值，这里的展示是经过转化后的
##数据的热图
#library(pheatmap)
#select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing = TRUE)[1:20]
#nt <- normTransform(dds) ##默认为logs(x+1),x为sizefactor处理后的counts
#log2.norm.counts <- assay(nt)[select,] ##选择前20个基因
#df <- as.data.frame(colData(dds)[,c("condition")])
#pheatmap(log2.norm.counts,cluster_rows = FALSE,show_rownames = FALSE,
#         cluster_cols = FALSE,annotation_col = df)
#pheatmap(assay(rld)[select,],cluster_rows = FALSE,show_rownames = FALSE,
#         cluster_cols = FALSE, annotation_col = df)
#dev.off()

##绘制样本间的距离热图
##dist(x, method="euclidean",diag=FALSE,upper=FALSE,p=2)默认参数
sampleDists <- dist(t(assay(rld))) 
##直接使用plot绘制距离树状图
pdf(file=paste0(output,"_hclust_dist.pdf"))
plot(hclust(sampleDists))
dev.off()

library(RColorBrewer)
library(pheatmap)
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pheatmap(sampleDistMatrix,clustering_distance_cols = sampleDists,
         clustering_distance_rows = sampleDists,col=colors)
dev.off()

##绘制PCA，距离相关矩阵就是样本的PCA图
##plotPCA(rld,intgroup="condition")
##使用gglot2绘制
library(ggplot2)
library(ggrepel)
data <- plotPCA(rld,intgroup="condition",returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"),digits = 3)
##aes尽量在需要用时才指定，ggplot里指定影响面太大
ggplot(data,aes(PC1,PC2))+
  geom_point(aes(color=condition,shape=condition),size=3)+
  geom_text_repel(aes(label=name))+labs(x=paste0("PC1(",percentVar[1],"%)"),
                                        y=paste0("PC2(",percentVar[2],"%)"),
                                        color="Type",shape="Type")
ggsave(paste0(output,"_PCA.pdf"),width=9.5,height=7)
dev.off()

##挑选出差异基因做GO/KEGG分析
DE_up <- Total_counts[Total_counts$log2FoldChange > 1 & 
                        (! is.na(Total_counts$padj)) & Total_counts$padj < 0.05,]

DE_down <- Total_counts[Total_counts$log2FoldChange < -1 & 
                                        (! is.na(Total_counts$padj)) & Total_counts$padj < 0.05,]

write.table(DE_up,file=paste0("ab_c507_0440_vs_ab_c","_DE_up_3vs2.txt"),sep="\t",quote=FALSE)
write.xlsx(DE_up,file =paste0("ab_c507_0440_vs_ab_c","_DE_up_3vs2.xlsx"),row.names=FALSE)
write.table(DE_down,file=paste0("ab_c507_0440_vs_ab_c","_DE_down_3vs2.txt"),sep="\t",quote=FALSE)
write.xlsx(DE_down,file =paste0("ab_c507_0440_vs_ab_c","_DE_down_3vs2.xlsx"),row.names=FALSE)

###Part 3
##topGO enrichment
library(topGO)
## get the GO annotation file and filter the duplicated GO terms in Genes
geneID2GO <- readMappings(go_file)
geneID2GO <- sapply(geneID2GO,function(var)unique(var))

geneList_up <- as.factor(as.integer(Total_counts$Gene_ID %in%
                                             DE_up$Gene_ID))
names(geneList_up) <- Total_counts$Gene_ID
geneList_down <- as.factor(as.integer(Total_counts$Gene_ID %in% 
                                        DE_down$Gene_ID))
names(geneList_down) <- Total_counts$Gene_ID 
go_type <- c("MF","BP","CC")

##DE up 构建富集GOdata数据
up_go <- list()
for(i in 1:length(go_type)){
  type=go_type[i]
  godata <- new("topGOdata",ontology=type,allGenes=geneList_up,
                description=paste("GOdata_up",type,sep="\t"),annot=annFUN.gene2GO,
                gene2GO=geneID2GO,nodeSize=1)
  ##renew the GOdata
  .geneList_up <- as.factor(as.integer(genes(godata) %in% sigGenes(godata)))
  names(.geneList_up) <- genes(godata)
  godata <- new("topGOdata", ontology=type,allGenes=.geneList_up,
                description=paste("GOdata_up",type,sep="\t"),annot=annFUN.gene2GO,
                gene2GO=geneID2GO,nodeSize=1)
  up_go[[i]] <- godata
}

#statistic <- c("classic","weight01","elim"), 富集检测
up_go_results <- list()
up_go_results_table <- list()
up_go_results_gentable <- list()
for(i in 1:length(up_go)){
  godata <- up_go[[i]]
  result <- runTest(godata,algorithm = "classic",statistic = "fisher")
  up_go_results_gentable[[i]] <- GenTable(godata,classic=result,orderBy="classic",ranksOf="classic",
                                                   topNodes=length(usedGO(godata)),numChar=50)
  up_go_results_gentable[[i]]$qvalue <- p.adjust(sort(score(result)),"BH")
  up_go_results_gentable[[i]] <- up_go_results_gentable[[i]][order(up_go_results_gentable[[i]]$qvalue),]
  up_go_results_table[[i]] <- up_go_results_gentable[[i]][1:30,]
  up_go_results_gentable[[i]]$Term <- ifelse(is.na(Definition(up_go_results_gentable[[i]]$GO.ID)),
                                         up_go_results_gentable[[i]]$Term,Definition(up_go_results_table[[i]]$GO.ID))
  up_go_results_gentable[[i]]$Sig_Genes <- sapply(
    sapply(genesInTerm(godata,up_go_results_gentable[[i]]$GO.ID),function(x){
      sigGenes(godata)[sigGenes(godata) %in% x]
    }),function(y){paste(y,collapse = ",")})
  up_go_results_table[[i]]$Sig_Genes <- up_go_results_gentable[[i]][1:30,]$Sig_Genes
  
  up_go_results_gentable[[i]]$All_Genes <- sapply(
    genesInTerm(godata,up_go_results_gentable[[i]]$GO.ID),
    function(x){paste(x,collapse = ",")})
  
  up_go_results[[i]] <- result
}

##DE down 构建GOdata数据
down_go <- list()
for(i in 1:length(go_type)){
  type=go_type[i]
  godata <- new("topGOdata",ontology=type,allGenes=geneList_down,
                description=paste("GOdata_down",type,sep="\t"),annot=annFUN.gene2GO,
                gene2GO=geneID2GO,nodeSize=1)
  ##renew the genelist
  .geneList_down <- as.factor(as.integer(genes(godata) %in% sigGenes(godata)))
  names(.geneList_down) <- genes(godata)
  godata <- new("topGOdata",ontology=type,allGenes=.geneList_down,
                description=paste("GOdata_down",type,sep="\t"),annot=annFUN.gene2GO,
                gene2GO=geneID2GO,nodeSize=1)
  down_go[[i]] <- godata
}

down_go_results <- list()
down_go_results_table <- list()
down_go_results_gentable <- list()                
for(i in 1:length(down_go)){
  godata <- down_go[[i]]
  result <- runTest(godata,algorithm = "classic",statistic = "fisher")
  down_go_results_gentable[[i]] <- GenTable(godata,classic=result,orderBy="classic",ranksOf="classic",
                                                   topNodes=length(usedGO(godata)),numChar=50)
  down_go_results_gentable[[i]]$qvalue <- p.adjust(sort(score(result)),"BH")
  down_go_results_gentable[[i]] <- down_go_results_gentable[[i]][order(down_go_results_gentable[[i]]$qvalue),]
  down_go_results_table[[i]] <- down_go_results_gentable[[i]][1:30,]
  down_go_results_gentable[[i]]$Term <- ifelse(is.na(Definition(down_go_results_gentable[[i]]$GO.ID)),
                                               down_go_results_gentable[[i]]$Term,Definition(down_go_results_table[[i]]$GO.ID))
  down_go_results_gentable[[i]]$Sig_Genes <- sapply(
    sapply(genesInTerm(godata,down_go_results_gentable[[i]]$GO.ID),function(x){
      sigGenes(godata)[sigGenes(godata) %in% x]
      }),function(y){paste(y,collapse = ",")})
  down_go_results_table[[i]]$Sig_Genes <- down_go_results_gentable[[i]][1:30,]$Sig_Genes
    
  down_go_results_gentable[[i]]$All_Genes <- sapply(
    genesInTerm(godata,down_go_results_gentable[[i]]$GO.ID),
    function(x){paste(x,collapse = ",")})

  down_go_results[[i]] <- result
}

##绘制GO 富集散点图,ggsave根据后缀判断类型(device)
library(ggplot2)
##创建富集结果输出目录
dir.create("./GO_enrichment_results")
for(i in 1:3){
  tmp=up_go_results_table[[i]]
  ##命名图片名称
  name=paste0(output,"_Up_",go_type[i],"_","Enrichment_Map")
  tmp$Annot_comb <- paste(tmp$GO.ID,tmp$Term,sep=" : ")
  tmp$qvalue <- as.numeric(tmp$qvalue)
  tmp$Significant <- as.numeric(tmp$Significant)
  tmp$Annot_comb <- factor(tmp$Annot_comb,levels = rev(tmp$Annot_comb))
  p<- ggplot(tmp,aes(qvalue,Annot_comb))+geom_point(aes(size=Significant,color=qvalue))+
    scale_color_gradient(low="red",high="green",limits=c(0,1))+
    labs(color="P.adjust",size="Significant Count",x="P.adjust",
         y="GO Terms")+theme(plot.title=element_text(hjust = 0.5))+theme_bw()
  ggsave(paste0("./GO_enrichment_results/",name,".pdf"),plot=p,width=25,height=15,units = "cm")
  printGraph(up_go[[i]],up_go_results[[i]],
                 firstSigNodes = 10,useInfo = "all",pdfSW=F,
                 fn.prefix=paste0("./GO_enrichment_results/",name,"DAG"))
  write.table(up_go_results_gentable[[i]],file=paste0(
    "./GO_enrichment_results/",name,".xls"),sep="\t",quote=F,row.names = F)
} 

for(i in 1:3){
  tmp=down_go_results_table[[i]]
  name=paste0(output,"_Down_",go_type[i],"_","Enrichment_Map")
  tmp$Annot_comb <- paste(tmp$GO.ID,tmp$Term,sep=" : ")
  tmp$qvalue <- as.numeric(tmp$qvalue)
  tmp$Significant <- as.numeric(tmp$Significant)
  tmp$Annot_comb <- factor(tmp$Annot_comb,levels = rev(tmp$Annot_comb))
  p<- ggplot(tmp,aes(qvalue,Annot_comb))+geom_point(aes(size=Significant,color=qvalue))+
    scale_color_gradient(low="red",high="green",limits=c(0,1))+
    labs(color="P.adjust",size="Significant Count",x="P.adjust",
         y="GO Terms")+theme(plot.title=element_text(hjust = 0.5))+theme_bw()
  ggsave(paste0("./GO_enrichment_results/",name,".pdf"),plot=p,width=25,height=15,units = "cm")
  printGraph(down_go[[i]],down_go_results[[i]],
             firstSigNodes = 10,useInfo = "all",pdfSW=F,
             fn.prefix=paste0("./GO_enrichment_results/",name,"DAG"))
  write.table(down_go_results_gentable[[i]],file=paste0(
    "./GO_enrichment_results/",name,".xls"),sep="\t",quote=F,row.names = F)
} 

##绘制netwrok和heatmap图
library(igraph)
library(fields)
library(reshape2)
library(ggplot2)
library(scales)
up_igraph_results_table <- list()
down_igraph_results_table <- list()

##DE up, 仅选取前15个显著性富集GO terms绘制
for(i in 1:3){
  num <- 0
  links <- data.frame()
  vertices <- data.frame()
  vertices <- up_go_results_table[[i]][1:15,c(1,2,4,7)]  ##前15个
  colnames(vertices) <- c("GO.ID","Term","Significant","qvalue")
  for(j in 1:(nrow(up_go_results_table[[i]][1:15,])-1)){  ##前15个
    from <- up_go_results_table[[i]][j,]$GO.ID
    from_genes <- unlist( strsplit( up_go_results_table[[i]][j,]$Sig_Genes,","))
    for(k in (j+1):nrow(up_go_results_table[[i]][1:15,])){  ##前15个
      end <- up_go_results_table[[i]][k,]$GO.ID
      end_genes <- unlist( strsplit( up_go_results_table[[i]][k,]$Sig_Genes,","))
      from_end_num <- sum(from_genes %in% end_genes)
      num <- num+1
      links[num,"from"] <- from
      links[num,"end"] <- end
      links[num,"count"] <- from_end_num
    }
  }
  links <- links[links$count>0,]
  up_igraph_results_table[[i]] <- list(vertices,links)
}
##draw up the netwroks
for(n in 1:3){
  tmp <- list()
  vertices <- data.frame()
  d <- data.frame()
  tmp <- up_igraph_results_table[[n]]
  vertices <- tmp[[1]]
  d <- tmp[[2]]
  net <- graph_from_data_frame(
    d=d,vertices=vertices,directed = FALSE)
  
  rbPal <- colorRampPalette(c("red","yellow"))
  V(net)$color <- rbPal(1000)[cut(seq(0,1,0.001),breaks = 1000)][
    round(vertices$qvalue,digits = 3)*1000+1] ##防止 qvalue=0
  
  V(net)$label <- vertices$Term
  V(net)$label.family <- "Times"
  V(net)$label.cex <- 0.6
  V(net)$label.dist <- 0.5
  E(net)$width <- d$count*0.15
  
  pdf(paste0("./GO_enrichment_results/",output,"_Up_",go_type[n],"_network.pdf"))
  plot(net,layout=layout_nicely)
  image.plot(zlim=seq(0,1),legend.only = TRUE,
             col=rbPal(1000)[cut(seq(0,1,by=0.001),breaks=1000)],
             horizontal = TRUE,legend.shrink=0.2,legend.width = 1,
             legend.args=list(text="P.adjust",cex=0.7,line=0.1),
             axis.args=list(at=c(0,0.5,1),labels=c(0,0.5,1),cex.axis=0.7,tck=-0.15,lwd=0,line=-0.95),
             smallplot = c(0.85,0.95,0.9,0.95))
  dev.off()
}

##DE down 同样绘制前15个GO terms
for(i in 1:3){
  num <- 0
  vertices <- data.frame()
  links <- data.frame()
  vertices <- down_go_results_table[[i]][1:15,c(1,2,4,7)]  ##前15个
  colnames(vertices) <- c("GO.ID","Term","Significant","qvalue")
  for(j in 1:(nrow(down_go_results_table[[i]][1:15,])-1)){  ##前15个
    from <- down_go_results_table[[i]][j,]$GO.ID
    from_genes <- unlist( strsplit( down_go_results_table[[i]][j,]$Sig_Genes,","))
    for(k in (j+1):nrow(down_go_results_table[[i]][1:15,])){  ##前15个
      end <- down_go_results_table[[i]][k,]$GO.ID
      end_genes <- unlist( strsplit( down_go_results_table[[i]][k,]$Sig_Genes,","))
      from_end_num <- sum(from_genes %in% end_genes)
      num <- num+1
      links[num,"from"] <- from
      links[num,"end"] <- end
      links[num,"count"] <- from_end_num
    }
  }
  links <- links[links$count>0,]
  down_igraph_results_table[[i]] <- list(vertices,links)
}

##draw down the netwroks
for(n in 1:3){
  tmp <- list()
  vertices <- data.frame()
  d <- data.frame()
  tmp <- down_igraph_results_table[[n]]
  vertices <- tmp[[1]]
  d <- tmp[[2]]
  net <- graph_from_data_frame(
    d=d,vertices=vertices,directed = F)
  
  rbPal <- colorRampPalette(c("red","yellow"))
  V(net)$color <- rbPal(1000)[cut(seq(0,1,by=0.001),breaks = 1000)][
    round(vertices$qvalue,digits=3)*1000+1] ##防止 qvalue=0
  
  V(net)$label <- vertices$Term
  V(net)$label.family <- "Times"
  V(net)$label.cex <- 0.6
  V(net)$label.dist <- 0.5
  E(net)$width <- d$count*0.15
  
  pdf(paste0("./GO_enrichment_results/",output,"_Down_",go_type[n],"_network.pdf"))
  plot(net,layout=layout_nicely)
  image.plot(zlim=seq(0,1),legend.only = TRUE,
             col=rbPal(1000)[cut(seq(0,1,by=0.001),breaks=1000)],
             horizontal = TRUE,legend.shrink=0.2,legend.width = 1,
             legend.args=list(text="P.adjust",cex=0.7,line=0.1),
             axis.args=list(at=c(0,0.5,1),labels=c(0,0.5,1),cex.axis=0.7,tck=-0.15,lwd=0,line=-0.95),
             smallplot = c(0.85,0.95,0.9,0.95))
  dev.off()
}

##up heatmap of GO terms with genes by names
for(m in 1:3){
  tmp <- data.frame()
  genes_up <- vector()
  Data <- data.frame()
  tmp=up_go_results_table[[m]]
  for(i in 1:nrow(tmp)){
    genes_up <- append(genes_up, unlist(strsplit(tmp[i,]$Sig_Genes,",")))
  }
  genes_up <- sort(unique(genes_up))
  
  Data <- data.frame(matrix(1:length(genes_up),nrow=1))
  for(j in 1:nrow(tmp)){
    Data[j,] <- as.integer(genes_up %in% unlist(strsplit(tmp[j,]$Sig_Genes,",")))
  }
  colnames(Data) <- factor(genes_up,levels=genes_up)
  rownames(Data) <- factor(tmp$GO.ID,levels=rev(tmp$GO.ID))
  x1 <- vector()
  x2 <- vector()
  y1 <- vector()
  y2 <- vector()
  q <- vector()
  d <- data.frame()
  p <- NULL
  for(k in 1:nrow(Data)){
    for(n in 1:ncol(Data)){
      x1 <- append(x1,Data[k,n]*(n-0.45))
      x2 <- append(x2,Data[k,n]*(n+0.45))
      y1 <- append(y1, Data[k,n]*(k-0.45))
      y2 <- append(y2, Data[k,n]*(k+0.45))
    }
    q <- append(q,rep(tmp[k,]$qvalue,length(genes_up)))
  }
  d <- data.frame(x1=x1,x2=x2,y1=y1,y2=y2,q=q)
  
  p <- ggplot() + theme_bw()+ geom_rect(
    data=d,mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=q))+
    scale_fill_gradient(low="red",high="blue",limits=c(0,1))+
    scale_y_continuous(breaks=seq(1,length(tmp$GO.ID)),labels=tmp$GO.ID,expand = c(0,0))+
    scale_x_continuous(breaks=seq(1,length(genes_up)),labels=genes_up,expand = c(0,0))+
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1,size=6.5),
          plot.title=element_text(hjust = 0.5))+
    labs(y="GO Terms",x="DE Genes",fill="P.adjust")
  ggsave(paste0("./GO_enrichment_results/",output,"_Up_",go_type[m],"_by_Name_heapmap.pdf"),
         plot=p,width = 28,height=18,units = "cm")
}

##up heatmap of GO terms with genes sorted by fold
for(m in 1:3){
  tmp <- data.frame()
  genes_up <- vector()
  Data <- data.frame()
  tmp=up_go_results_table[[m]]
  for(i in 1:nrow(tmp)){
    genes_up <- append(genes_up, unlist(strsplit(tmp[i,]$Sig_Genes,",")))
  }
  genes_up <- sort(unique(genes_up))
  genes_up <- genes_up[order(DE_up[genes_up,"log2FoldChange"],
                             decreasing = TRUE)]
  
  Data <- data.frame(matrix(1:length(genes_up),nrow=1))
  for(j in 1:nrow(tmp)){
    Data[j,] <- as.integer(genes_up %in% unlist(strsplit(tmp[j,]$Sig_Genes,",")))
  }
  colnames(Data) <- factor(genes_up,levels=genes_up)
  rownames(Data) <- factor(tmp$GO.ID,levels=rev(tmp$GO.ID))
  x1 <- vector()
  x2 <- vector()
  y1 <- vector()
  y2 <- vector()
  q <- vector()
  d <- data.frame()
  p <- NULL
  for(k in 1:nrow(Data)){
    for(n in 1:ncol(Data)){
      x1 <- append(x1,Data[k,n]*(n-0.45))
      x2 <- append(x2,Data[k,n]*(n+0.45))
      y1 <- append(y1, Data[k,n]*(k-0.45))
      y2 <- append(y2, Data[k,n]*(k+0.45))
    }
    q <- append(q,rep(tmp[k,]$qvalue,length(genes_up)))
  }
  d <- data.frame(x1=x1,x2=x2,y1=y1,y2=y2,q=q)
  
  p <- ggplot() + theme_bw()+ geom_rect(
    data=d,mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=q))+
    scale_fill_gradient(low="red",high="blue",limits=c(0,1))+
    scale_y_continuous(breaks=seq(1,length(tmp$GO.ID)),labels=tmp$GO.ID,expand = c(0,0))+
    scale_x_continuous(breaks=seq(1,length(genes_up)),labels=genes_up,expand = c(0,0))+
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1,size=6.5),
          plot.title=element_text(hjust = 0.5))+
    labs(y="GO Terms",x="DE Genes",fill="P.adjust")
  ggsave(paste0("./GO_enrichment_results/",output,"_Up_",go_type[m],"_by_Fold_heapmap.pdf"),
         plot=p,width = 28,height=18,units = "cm")
}


##down heatmap of GO terms with genes by names
for(m in 1:3){
  tmp <- data.frame()
  genes_down <- vector()
  Data <- data.frame()
  tmp=down_go_results_table[[m]]
  for(i in 1:nrow(tmp)){
    genes_down <- append(genes_down, unlist(strsplit(tmp[i,]$Sig_Genes,",")))
  }
  genes_down <- sort(unique(genes_down))
  
  Data <- data.frame(matrix(1:length(genes_down),nrow=1))
  for(j in 1:nrow(tmp)){
    Data[j,] <- as.integer(genes_down %in% unlist(strsplit(tmp[j,]$Sig_Genes,",")))
  }
  colnames(Data) <- factor(genes_down,levels=genes_down)
  rownames(Data) <- factor(tmp$GO.ID,levels=rev(tmp$GO.ID))
  x1 <- vector()
  x2 <- vector()
  y1 <- vector()
  y2 <- vector()
  q <- vector()
  p <- NULL
  for(k in 1:nrow(Data)){
    for(n in 1:ncol(Data)){
      x1 <- append(x1,Data[k,n]*(n-0.45))
      x2 <- append(x2,Data[k,n]*(n+0.45))
      y1 <- append(y1, Data[k,n]*(k-0.45))
      y2 <- append(y2, Data[k,n]*(k+0.45))
    }
    q <- append(q,rep(tmp[k,]$qvalue,length(genes_down)))
  }
  d <- data.frame(x1=x1,x2=x2,y1=y1,y2=y2,q=q)
  
  p <- ggplot() + theme_bw()+ geom_rect(
    data=d,mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=q))+
    scale_fill_gradient(low="red",high="blue",limits=c(0,1))+
    scale_y_continuous(breaks=seq(1,length(tmp$GO.ID)),labels=tmp$GO.ID,expand = c(0,0))+
    scale_x_continuous(breaks=seq(1,length(genes_down)),labels=genes_down,expand = c(0,0))+
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1,size=6.5),
          plot.title=element_text(hjust = 0.5))+
    labs(y="GO Terms",x="DE Genes",fill="P.adjust")
  ggsave(paste0("./GO_enrichment_results/",output,"_Down_",go_type[m],"_by_Name_heapmap.pdf"),
         plot=p,width = 28,height=18,units = "cm")
}

##down heatmap of GO terms with genes by fold
for(m in 1:3){
  tmp <- data.frame()
  genes_down <- vector()
  Data <- data.frame()
  tmp=down_go_results_table[[m]]
  for(i in 1:nrow(tmp)){
    genes_down <- append(genes_down, unlist(strsplit(tmp[i,]$Sig_Genes,",")))
  }
  genes_down <- sort(unique(genes_down))
  genes_down <- genes_down[order(DE_down[genes_down,"log2FoldChange"])]
  
  Data <- data.frame(matrix(1:length(genes_down),nrow=1))
  for(j in 1:nrow(tmp)){
    Data[j,] <- as.integer(genes_down %in% unlist(strsplit(tmp[j,]$Sig_Genes,",")))
  }
  colnames(Data) <- factor(genes_down,levels=genes_down)
  rownames(Data) <- factor(tmp$GO.ID,levels=rev(tmp$GO.ID))
  x1 <- vector()
  x2 <- vector()
  y1 <- vector()
  y2 <- vector()
  q <- vector()
  p <- NULL
  for(k in 1:nrow(Data)){
    for(n in 1:ncol(Data)){
      x1 <- append(x1,Data[k,n]*(n-0.45))
      x2 <- append(x2,Data[k,n]*(n+0.45))
      y1 <- append(y1, Data[k,n]*(k-0.45))
      y2 <- append(y2, Data[k,n]*(k+0.45))
    }
    q <- append(q,rep(tmp[k,]$qvalue,length(genes_down)))
  }
  d <- data.frame(x1=x1,x2=x2,y1=y1,y2=y2,q=q)
  
  p <- ggplot() + theme_bw()+ geom_rect(
    data=d,mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=q))+
    scale_fill_gradient(low="red",high="blue",limits=c(0,1))+
    scale_y_continuous(breaks=seq(1,length(tmp$GO.ID)),labels=tmp$GO.ID,expand = c(0,0))+
    scale_x_continuous(breaks=seq(1,length(genes_down)),labels=genes_down,expand = c(0,0))+
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1,size=6.5),
          plot.title=element_text(hjust = 0.5))+
    labs(y="GO Terms",x="DE Genes",fill="P.adjust")
  ggsave(paste0("./GO_enrichment_results/",output,"_Down_",go_type[m],"_by_Fold_heapmap.pdf"),
         plot=p,width = 28,height=18,units = "cm")
}

##part 4 
##clusterProfiler kegg analysis
library(AnnotationDbi)
library(AnnotationForge)
library(AnnotationHub)
library(clusterProfiler)
library(enrichplot)

##abau for Acinetobacter baumannii AB0303
##KEGG Module是人工审核定义的功能单元，在一些情况下，KEGG Module具有明确直接的解释
##https://guangchuangyu.github.io/2016/04/kegg-module-enrichment-analysis/
##https://www.genome.jp/kegg/catalog/org_list.html
dir.create("./KEGG_enrichment_results")
locus_exchange <- read.delim("ab030_locustag_exchange.txt",header=F,sep="\t",quote="",row.names = 1)
DE_up$old_Gene_ID <- locus_exchange[DE_up$Gene_ID,]
DE_down$old_Gene_ID <- locus_exchange[DE_down$Gene_ID,]

kegg_DE_up <- enrichMKEGG(gene=DE_up$old_Gene_ID,organism = "abau",keyType ="kegg",minGSSize = 1,
                          pAdjustMethod="BH",pvalueCutoff =2,qvalueCutoff = 2)
kegg_DE_down <- enrichMKEGG(gene=DE_down$old_Gene_ID,organism = "abau",keyType="kegg",minGSSize= 1,
                            pAdjustMethod="BH",pvalueCutoff = 2,qvalueCutoff = 2)

kegg_DE_up_dataframe <- data.frame(kegg_DE_up)
kegg_DE_down_dataframe <- data.frame(kegg_DE_down)

write.table(kegg_DE_up_dataframe,file=paste0("./KEGG_enrichment_results/",output,"_KEGG_up_3vs2.txt"),
            sep="\t",quote=FALSE)
write.xlsx(kegg_DE_up_dataframe,file =paste0("./KEGG_enrichment_results/",output,"_KEGG_up_3vs2.xlsx"),
           row.names=FALSE)
write.table(kegg_DE_down_dataframe,file=paste0("./KEGG_enrichment_results/",output,"_KEGG_down_3vs2.txt"),
            sep="\t",quote=FALSE)
write.xlsx(kegg_DE_down_dataframe,file =paste0("./KEGG_enrichment_results/",output,"_KEGG_down_3vs2.xlsx"),
           row.names=FALSE)

#kegg_DE_up_simple <- enrichKEGG(gene=DE_up$old_Gene_ID,organism = "abau",keyType ="kegg",minGSSize = 1,
#                          pAdjustMethod="BH",pvalueCutoff =2,qvalueCutoff = 2)
#kegg_DE_down_simple <- enrichKEGG(gene=DE_down$old_Gene_ID,organism = "abau",keyType="kegg",minGSSize= 1,
#                            pAdjustMethod="BH",pvalueCutoff = 2,qvalueCutoff = 2)

##browseKEGG
##browseKEGG(kk, 'hsa04110')
##library("pathview")
#hsa04110 <- pathview(gene.data  = geneList,
#                     pathway.id = "hsa04110",
#                     species    = "hsa",
#                     limit      = list(gene=max(abs(geneList)), cpd=1))
