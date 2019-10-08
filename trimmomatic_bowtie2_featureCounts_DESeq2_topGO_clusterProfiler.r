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
####featureCounts ignoreDup=F,countMultiMappingReads=F
##featureCounts
library(Rsubread)
Results <- list()
for(sample in c("PAO1","Y89","Y71","Y82","Y31")){
  tmp <- c()
  tmp <- paste0(sample,"_sorted.bam")
  Results[[sample]] <- featureCounts(tmp,annot.ext="GCF_000006765.1_ASM676v1_genomic.gtf",
                                     isGTFAnnotationFile = TRUE,GTF.featureType = "gene",
                                     GTF.attrType = "gene_id",isPairedEnd = FALSE,ignoreDup = FALSE,
                                     nthreads = 4,countChimericFragments = FALSE,
                                     countMultiMappingReads = FALSE)
}
##paired-end reads featureCounts
#sample <- featureCounts("sample.bam",annot.ext = "HS11286.gtf",
#                        isGTFAnnotationFile = TRUE,GTF.featureType = "transcript",
#                        GTF.attrType = "Name",isPairedEnd = TRUE,ignoreDup=FALSE,
#                        requireBothEndsMapped=TRUE,nthreads=4,countChimericFragments=FALSE,
#                        countMultiMappingReads=FALSE)

###collect counts
fc_counts <- list()
for(sample in c("PAO1","Y89","Y71","Y82","Y31")){
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
##collection and print to file
##因为重名，所有.1后缀为对应的TPM值
Total_counts <- data.frame(row.names=rownames(fc_counts[[1]]),
                           Gene_ID=fc_counts[[1]]$Gene_ID,
                           Pos=fc_counts[[1]]$Pos,
                           Length=fc_counts[[1]]$Length,
                           sapply(fc_counts,function(x)x$Count),
                           sapply(fc_counts,function(x)x$TPM)
)
write.table(Total_counts,file="All_samples_count_statistic.txt",sep="\t",quote=F)
write.csv(Total_counts,file="All_samples_count_statistic.csv")

##构建DESeq2 countData/colData
countData <- sapply(fc_counts,function(x)x$Count)
rownames(countData) <- fc_counts[[1]]$Gene_ID
colData <- matrix(666,nrow=length(fc_counts))
rownames(colData) <- colnames(countData)
colData[,1] <- c("AS","MDR","MDR","MDR","AS")
colnames(colData) <- "condition"

##DESeq2差异分析
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,
                              design = ~ condition)
##pre-filtering，过滤掉low count的genes(没有read或几乎没有read)
##此过程非必要，后续会自动过滤，为实现普遍结论，而非个别差异，过滤read count
##edgeR 建议5-10个方为表达
dds <- dds[rowSums(counts(dds)>=5)==5,]

##设置factor levels
dds$condition <- factor(dds$condition,levels=c("MDR","AS"))
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
res <- results(dds,contrast = c("condition","MDR","AS"),alpha=0.05)
##removeResults函数返回至DESeqDataSet对象

##结果根据padj排序
resOrdered <- res[order(res$padj),]

##lfcshrink, 仅为绘图和排序使用，最终差异结果和results无异
#res_shrunken_normal <- lfcShrink(dds,contrast = c("condition","MDR","AS"),
#res=res,type="normal",alpha=0.05)

##绘制MA-plot
plotMA(res,main="DESeq2",ylim=c(-2,2))
##idx <- identity(res$baseMean,res$log2FoldChange)
##rownames(res)[idex] 交互式获得对应点坐标信息

##count plots,检查单个基因的count read
plotCounts(dds,gene=which.min(res$padj),intgroup = "condition")

##查看rsults返回结果描述
##p-values == NA
##1，一行中，所有样本counts都为0,baseMean为0，其他都为NA
##2，若一行中包含极端count至，p和adjusted p都为NA
##outlier counts的检出是根据Cook's distance
##3，若一行由于低于mean normalized count
##而被automatic independent filtering, 那么只有adjusted p为NA
##mcols(resOrdered)$description
##详细查看metadata(res)
##sizeFacots(dds)
##coef(dds)
##

##整理结果，绘制PCA/heatmap,volcano,DEgenes





###Part 3
##topGO enrichment
library(topGO)
geneID2GO <- readMappings("HS11286_sangon_go.txt")
geneList_up <- as.factor(as.integer(rownames(cts) %in%
                                             DE_28vs21_up$GeneNames))
names(geneList_up) <- rownames(cts)
geneList_down <- as.factor(as.integer(rownames(cts) %in% 
                                        DE_28vs21_down$GeneNames))
names(geneList_down) <- rownames(cts)  
go_type <- c("MF","BP","CC")

##28 vs 21 high 构建富集GOdata数据
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

##28 vs 21 down 构建GOdata数据
down_go <- list()
for(i in 1:length(go_type)){
  type=go_type[i]
  godata <- new("topGOdata",ontology=type,allGenes=geneList_down,
                description=paste("GOdata_28vs21_l",type,sep="\t"),annot=annFUN.gene2GO,
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
  name=paste0("Kp_28vs21_Up","_",go_type[i],"_","Enrichment_Map")
  tmp$Annot_comb <- paste(tmp$GO.ID,tmp$Term,sep=" : ")
  tmp$qvalue <- as.numeric(tmp$qvalue)
  tmp$Significant <- as.numeric(tmp$Significant)
  tmp$Annot_comb <- factor(tmp$Annot_comb,levels = rev(tmp$Annot_comb))
  p<- ggplot(tmp,aes(qvalue,Annot_comb))+geom_point(aes(size=Significant,color=qvalue))+
    scale_color_gradient(low="red",high="green")+scale_x_reverse()+
    labs(color="Classic Fisher Qvalue",size="Significant Count",x="Classic Fisher Qvalue",
         y="GO Terms",title=name)+theme(plot.title=element_text(hjust = 0.5))+theme_bw()
  ggsave(paste0("./GO_enrichment_results/",name,".pdf"),plot=p,width=25,height=15,units = "cm")
  printGraph(up_go[[i]],up_go_results[[i]],
                 firstSigNodes = 10,useInfo = "all",pdfSW=F,
                 fn.prefix=paste0("./GO_enrichment_results/",name,"DAG"))
  write.table(up_go_results_gentable[[i]],file=paste0(
    "./GO_enrichment_results/",name,".xls"),sep="\t",quote=F,row.names = F)
} 

for(i in 1:3){
  tmp=down_go_results_table[[i]]
  name=paste0("Kp_28vs21_Down","_",go_type[i],"_","Enrichment_Map")
  tmp$Annot_comb <- paste(tmp$GO.ID,tmp$Term,sep=" : ")
  tmp$qvalue <- as.numeric(tmp$qvalue)
  tmp$Significant <- as.numeric(tmp$Significant)
  tmp$Annot_comb <- factor(tmp$Annot_comb,levels = rev(tmp$Annot_comb))
  p<- ggplot(tmp,aes(qvalue,Annot_comb))+geom_point(aes(size=Significant,color=qvalue))+
    scale_color_gradient(low="red",high="green")+scale_x_reverse()+
    labs(color="Classic Fisher Qvalue",size="Significant Count",x="Classic Fisher Qvalue",
         y="GO Terms",title=name)+theme(plot.title=element_text(hjust = 0.5))+theme_bw()
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
  rbPal <- colorRampPalette(c("yellow","red"))
  V(net)$color <- rbPal(100)[as.numeric(cut(rescale(-log10(vertices$qvalue)),breaks = 100))]
  V(net)$label <- vertices$Term
  V(net)$label.family <- "Times"
  V(net)$label.cex <- 0.6
  V(net)$label.dist <- 0.5
  E(net)$width <- d$count*0.15
  
  #pdf(paste0("./GO_enrichment_results/","28vs21_up_",go_type[n],"_star_network.pdf"))
  #plot(net,layout=layout_as_star,main=paste0(
  #  "28vs21_up_",go_type[n],"_star_network"))
  #dev.off()
  
  pdf(paste0("./GO_enrichment_results/","28vs21_up_",go_type[n],"_network.pdf"))
  plot(net,layout=layout_nicely,main=paste0(
    "28vs21_up_",go_type[n],"_network"))
  image.plot(legend.only = TRUE,zlim=range(rescale(-log10(vertices$qvalue))),
             col=rbPal(100)[cut(seq(0,1,by=0.001),breaks=100)],
             horizontal = TRUE,legend.shrink=0.2,legend.width = 1,
             legend.args=list(text=expression(-log[10](Qvalue)),cex=0.7,line=0.1),
             axis.args=list(at=c(0,0.5,1),labels=c(0,0.5,1),line=0.05),
             smallplot = c(0.8,0.9,0.85,0.9))
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
  rbPal <- colorRampPalette(c("yellow","red"))
  V(net)$color <- rbPal(100)[as.numeric(cut(rescale(-log10(vertices$qvalue)),breaks = 100))]
  V(net)$label <- vertices$Term
  V(net)$label.family <- "Times"
  V(net)$label.cex <- 0.6
  V(net)$label.dist <- 0.5
  E(net)$width <- d$count*0.15
  
  #pdf(paste0("./GO_enrichment_results/","28vs21_down_",go_type[n],"_star_network.pdf"))
  #plot(net,layout=layout_as_star,main=paste0(
  #  "28vs21_down_",go_type[n],"_star_network"))
  #dev.off()
  
  pdf(paste0("./GO_enrichment_results/","28vs21_down_",go_type[n],"_network.pdf"))
  plot(net,layout=layout_nicely,main=paste0(
    "28vs21_down_",go_type[n],"_network"))
  image.plot(legend.only = TRUE,zlim=range(rescale(-log10(vertices$qvalue))),
             col=rbPal(100)[cut(seq(0,1,by=0.001),breaks=100)],
             horizontal = TRUE,legend.shrink=0.2,legend.width = 1,
             legend.args=list(text=expression(-log[10](Qvalue)),cex=0.7,line=0.1),
             axis.args=list(at=c(0,0.5,1),labels=c(0,0.5,1),line=0.05),
             smallplot = c(0.8,0.9,0.85,0.9))
  dev.off()
}

##up heatmap of GO terms with genes
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
    scale_fill_gradient(low="red",high="blue")+
    scale_y_continuous(breaks=seq(1,length(tmp$GO.ID)),labels=tmp$GO.ID,expand = c(0,0))+
    scale_x_continuous(breaks=seq(1,length(genes_up)),labels=genes_up,expand = c(0,0))+
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1,size=6.5),
          plot.title=element_text(hjust = 0.5))+
    labs(title=paste0(go_type[m],"_Heatmap"),y="GO Terms",x="DE Genes",fill="Qvalue")
  ggsave(paste0("./GO_enrichment_results/","Kp_28vs21_Up_",go_type[m],"_heapmap.pdf"),
         plot=p,width = 28,height=18,units = "cm")
}

##down heatmap of GO terms with genes
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
    scale_fill_gradient(low="red",high="blue")+
    scale_y_continuous(breaks=seq(1,length(tmp$GO.ID)),labels=tmp$GO.ID,expand = c(0,0))+
    scale_x_continuous(breaks=seq(1,length(genes_down)),labels=genes_down,expand = c(0,0))+
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1,size=6.5),
          plot.title=element_text(hjust = 0.5))+
    labs(title=paste0(go_type[m],"_Heatmap"),y="GO Terms",x="DE Genes",fill="Qvalue")
  ggsave(paste0("./GO_enrichment_results/","Kp_28vs21_Down_",go_type[m],"_heapmap.pdf"),
         plot=p,width = 28,height=18,units = "cm")
}

##part 4 
##clusterProfiler kegg analysis
library(AnnotationDbi)
library(AnnotationForge)
library(AnnotationHub)
library(clusterProfiler)
library(enrichplot)
##kegg_database <- search_kegg_organism("Klebsiella pneumoniae",by="scientific_name",
##                                      ignore.case = T)
##kpm
##KEGG Module是人工审核定义的功能单元，在一些情况下，KEGG Module具有明确直接的解释
kegg_28vs21_up <- enrichMKEGG(gene=DE_28vs21_up$GeneNames,organism = "kpm",
                            pvalueCutoff = 1,minGSSize= 1,qvalueCutoff = 1)
kegg_28vs21_down <- enrichMKEGG(gene=DE_28vs21_down$GeneNames,organism = "kpm",
                             pvalueCutoff = 1,minGSSize = 1,qvalueCutoff = 1)

##dotplot pic / 
kegg_28vs21_up <- as.data.frame(kegg_28vs21_up)
kegg_28vs21_down <- as.data.frame(kegg_28vs21_down)
kegg_28vs21_up$Terms <- paste0(kegg_28vs21_up$ID,":",kegg_28vs21_up$Description)
kegg_28vs21_down$Terms <- paste0(kegg_28vs21_down$ID,":",kegg_28vs21_down$Description)
kegg_28vs21_results <- list(kegg_28vs21_up,kegg_28vs21_down)
name_kegg <- c("KEGG_28vs21_Up","KEGG_28vs21_Down")

for(i in 1:2){
  tmp <- data.frame()
  name <- c()
  if(nrow(kegg_28vs21_results[[i]])>30){
    tmp <- kegg_28vs21_results[[i]][1:30,]
    }else{
      tmp <- kegg_28vs21_results[[i]]
    }
  name <- name_kegg[i]
  tmp$p.adjust <- as.numeric(tmp$p.adjust)
  tmp$Count <- as.numeric(tmp$Count)
  tmp$Terms <- factor(tmp$Terms,levels = rev(tmp$Terms))
  p<- ggplot(tmp,aes(p.adjust,Terms))+geom_point(aes(size=Count,color=p.adjust))+
    scale_color_gradient(low="red",high="green")+scale_x_reverse()+
    labs(color="p.adjust",size="Significant Count",x="p.adjust",
         y="KO Terms",title=paste0(name,"_Enrichment"))+theme(plot.title=element_text(hjust = 0.5))+theme_bw()
  ggsave(paste0("./GO_enrichment_results/",name,"_Enrichment",".pdf"),plot=p,width=25,height=15,units = "cm")
} 

##heatmap for kegg enrichment 
for(m in 1:2){
  tmp <- data.frame()
  genes <- vector()
  Data <- data.frame()
  if(nrow(kegg_28vs21_results[[m]])>30){
    tmp <- kegg_28vs21_results[[m]][1:30,]
  }else{
    tmp <- kegg_28vs21_results[[m]]
  }
  
  for(i in 1:nrow(tmp)){
    genes <- append(genes, unlist(strsplit(tmp$geneID,"/")))
  }
  genes <- sort(unique(genes))
  
  Data <- data.frame(matrix(1:length(genes),nrow=1))
  for(j in 1:nrow(tmp)){
    Data[j,] <- as.integer(genes %in% unlist(strsplit(tmp[j,]$geneID,"/")))
  }
  colnames(Data) <- factor(genes,levels=genes)
  rownames(Data) <- factor(tmp$ID,levels=rev(tmp$ID))
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
    q <- append(q,rep(tmp[k,]$p.adjust,length(genes)))
  }
  d <- data.frame(x1=x1,x2=x2,y1=y1,y2=y2,q=q)
  
  p <- ggplot() + theme_bw()+ geom_rect(
    data=d,mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=q))+
    scale_fill_gradient(low="red",high="blue")+
    scale_y_continuous(breaks=seq(1,length(tmp$ID)),labels=tmp$ID,expand = c(0,0))+
    scale_x_continuous(breaks=seq(1,length(genes)),labels=genes,expand = c(0,0))+
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1,size=6.5),
          plot.title=element_text(hjust = 0.5))+
    labs(title=paste0(name_kegg[m],"_Heatmap"),y="KO Terms",x="DE Genes",fill="p.adjust")
  ggsave(paste0("./GO_enrichment_results/",name_kegg[m],"_heapmap.pdf"),
         plot=p,width = 28,height=18,units = "cm") 
}

##browseKEGG
##browseKEGG(kk, 'hsa04110')
##library("pathview")
#hsa04110 <- pathview(gene.data  = geneList,
#                     pathway.id = "hsa04110",
#                     species    = "hsa",
#                     limit      = list(gene=max(abs(geneList)), cpd=1))
