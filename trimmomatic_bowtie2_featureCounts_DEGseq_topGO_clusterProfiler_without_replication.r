###part 1
##trimmomatic
#trimmomatic PE -summary C_trimmed.log C_R1.fastq.gz C_R2.fastq.gz -baseout C_trimmed.fastq.gz \
#ILLUMINACLIP:Sabgon_adapter.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:35 
##sangon_adapter.fa
##>PrefixPE/1
##AGATCGGAAGAGCACACGTCTGAAC
##>PrefixPE/2
##AGATCGGAAGAGCGTCGTGTAGGGA

###align bowtie2 with shell scripts
#bowtie2-build GCF_000240185.1_ASM24018v2_genomic.fna hs11S6_index
#bowtie2 -x hs11S6_index -1 C_sangon_1P.fastq.gz -2 C_sangon_2P.fastq.gz -S C_sangon_mapped.sam
#bowtie2 -x hs11S6_index -1 S_sangon_1P.fastq.gz -2 S_sangon_2P.fastq.gz -S S_sangon_mapped.sam
#samtools view -bS C_sangon_mapped.sam > C_sangon_mapped.bam
#samtools view -bS S_sangon_mapped.sam > S_sangon_mapped.bam
#samtools sort -n C_sangon_mapped.bam -o C_sangon_mapped_sorted.bam
#samtools sort -n S_sangon_mapped.bam -o S_sangon_mapped_sorted.bam

##part 2
##考虑到基因重组，仅计算primary 比对即可
####featureCounts ignoreDup=F,primaryOnly=TRUE
##featureCounts
library(Rsubread)
##modify the input according to the reality, strandSpecific=2
sample_s <- "KP_28_sangon_mapped_sorted.bam"
sample_c <- "KP_21_sangon_mapped_sorted.bam"
gtf_file <- "HS11286.gtf"
featuretype <- "transcript"
attrtype <- "locus_tag"
anno_file <- "HS11286_annotation_extraction.txt"
#go_file <- "HS11286_sangon_go.txt"
go_file <- "atcc13883_uniprot_go.txt"
output <- "KP28vsKP21"
##
S_exp <- featureCounts(sample_s,annot.ext = gtf_file,
                           isGTFAnnotationFile = TRUE,GTF.featureType = featuretype,
                           GTF.attrType = attrtype,isPairedEnd = TRUE,allowMultiOverlap=TRUE,
                           primaryOnly=TRUE,strandSpecific=2,nthreads=4)
C_exp <- featureCounts(sample_c,annot.ext = gtf_file,
                           isGTFAnnotationFile = TRUE,GTF.featureType = featuretype,
                           GTF.attrType = attrtype,isPairedEnd = TRUE,allowMultiOverlap=TRUE,
                           primaryOnly=TRUE,strandSpecific=2,nthreads=4)

###collect counts
fc_S <- data.frame(Gene_ID=rownames(S_exp$counts),Pos=paste0(
  S_exp$annotation$Chr,"[",S_exp$annotation$Strand,"]",
  S_exp$annotation$Start,"-",S_exp$annotation$End),
  Length=S_exp$annotation$Length,Count=S_exp$counts[,1])
fc_C <- data.frame(Gene_ID=rownames(C_exp$counts),Pos=paste0(
  C_exp$annotation$Chr,"[",C_exp$annotation$Strand,"]",
  C_exp$annotation$Start,"-",C_exp$annotation$End),
  Length=C_exp$annotation$Length,Count=C_exp$counts[,1])

##选择不过滤，最后得出差异再过滤
library(edgeR)
library(xlsx)
cts <- data.frame(S=fc_S$Count,C=fc_C[rownames(fc_S),]$Count,
                  row.names=rownames(fc_S))
group=factor(c("S","C"),levels=c("S","C"))
y <- DGEList(cts,group=group)
#y$samples$lib.size=colSums(y$counts)
y <- calcNormFactors(y,method="TMM")

#norm.count
fc_S$Norm.Count <- fc_S$Count/y$samples$norm.factors[1]
fc_C$Norm.Count <- fc_C$Count/y$samples$norm.factors[2]
##TPM
fc_S$TPM <- ((fc_S$Norm.Coun/fc_S$Length)*1000000)/
  (sum(fc_S$Norm.Coun/fc_S$Length))
fc_C$TPM <- ((fc_C$Norm.Coun/fc_C$Length)*1000000)/
  (sum(fc_C$Norm.Coun/fc_C$Length))
##updated fc_sangon cpm rpkm normalized count
cpm_sangon <- data.frame(cpm(y))
fc_S$CPM <- cpm_sangon[rownames(fc_S),]$S
fc_C$CPM <- cpm_sangon[rownames(fc_C),]$C
#length
rpkm_sangon <- data.frame(rpkm(y,gene.length =
  c(fc_S[rownames(y$counts),]$Length,
    fc_C[rownames(y$counts),]$Length)))
fc_S$RPKM <- rpkm_sangon[rownames(fc_S),]$S
fc_C$RPKM <- rpkm_sangon[rownames(fc_C),]$C

##collection and print to file, add the annotation contents
Anno <- read.delim(anno_file,sep="\t",header=FALSE,
                   row.names = 3)
Total_counts <- data.frame(row.names=rownames(fc_S),
                           Geneid=rownames(fc_S),
                           Pos=fc_S$Pos,
                           Length=fc_S$Length,
                           Counts_S=fc_S$Count,
                           Counts_C=fc_C[rownames(fc_S),]$Count,
                           CPM_S=fc_S$CPM,
                           CPM_C=fc_C[rownames(fc_S),]$CPM,
                           RPKM_S=fc_S$RPKM,
                           RPKM_C=fc_C[rownames(fc_S),]$RPKM,
                           TPM_S=fc_S$TPM,
                           TPM_C=fc_C[rownames(fc_S),]$TPM,
                           Norm.Count_S=fc_S$Norm.Count,
                           Norm.Count_C=fc_C[rownames(fc_S),]$Norm.Count,
                           Anno=Anno[rownames(fc_S),]$V5
)

##手动修改样本名称或直接全部设置变量提前修改
##手动修改样本名称或直接全部设置变量提前修改
colnames(Total_counts) <- sub("(*\\_)S","\\1KP28",colnames(Total_counts))
colnames(Total_counts) <- sub("(*\\_)C","\\1KP21",colnames(Total_counts))
#write.table(Total_counts,file=paste0(output,"_count_statistic.txt"),sep="\t",quote=FALSE,
#            row.names=FALSE)
#write.xlsx(Total_counts,file =paste0(output,"_count_statistic.xlsx"),row.names=FALSE)

### cts count normalized by norm.factors and is used for DEGseq package analysis
library(DEGseq)
cts$S_norm <- cts$S/(y$samples$norm.factors[1])
cts$C_norm <- cts$C/(y$samples$norm.factors[2])

write.table(cts,file="cts_normalized_by_norm_factor.txt",
            sep="\t",row.names = TRUE)
S_cts <- readGeneExp(file="cts_normalized_by_norm_factor.txt",geneCol=1,valCol = 4)
C_cts <- readGeneExp(file="cts_normalized_by_norm_factor.txt",geneCol=1,valCol = 5)
layout(matrix(c(1,2,3,4,5,6),3,2,byrow=T))
par(mar=c(2,2,2,2))
DEGexp(geneExpMatrix1 = S_cts,geneCol1 = 1,expCol1 = 2,groupLabel1 = "S_cts",
       geneExpMatrix2 = C_cts,geneCol2 = 1,expCol2 = 2,groupLabel2 = "C_cts",
       method="MARS",rawCount = F,thresholdKind=3,qValue=0.05,
       outputDir="./sangon_cts_norm_degseq_results")

##topGO analysis 
score <- read.delim("./sangon_cts_norm_degseq_results/output_score.txt",
                      header=T,sep="\t")
rownames(score) <- score$GeneNames

##选择后过滤TMP>=5
score <- cbind(score,Total_counts[rownames(score),2:ncol(Total_counts)])

##log2.Fold_change. 对应修改TPM名称
DE_SvsC_up <- score[score$log2.Fold_change. > 1 & 
                           score$q.value.Benjamini.et.al..1995. < 0.05 &
                           score$TPM_KP28 >= 5,]
DE_SvsC_down <- score[score$log2.Fold_change. < -1 &
                             score$q.value.Benjamini.et.al..1995. < 0.05 &
                             score$TPM_KP21 >= 5,]

DE_SvsC_up <- DE_SvsC_up[,c(1,11,12,13,14,21,22,19,20,4,7,8,10,23)]
DE_SvsC_down <- DE_SvsC_down[,c(1,11,12,13,14,21,22,19,20,4,7,8,10,23)]
DE_SvsC_total <- score[,c(1,11,12,13,14,21,22,19,20,4,7,8,10,23)]
write.table(DE_SvsC_down,file=paste0(output,"_DE_Down.txt"),sep="\t",quote = FALSE,row.names=FALSE)
write.table(DE_SvsC_up,file=paste0(output,"_DE_Up.txt"),sep="\t",quote = FALSE,row.names=FALSE)
write.table(DE_SvsC_total,file=paste0(output,"_DE_Total.txt"),sep="\t",quote=FALSE,row.names=FALSE)
write.xlsx(DE_SvsC_down,file=paste0(output,"_DE_Down.xlsx"),row.names=FALSE)
write.xlsx(DE_SvsC_up,file=paste0(output,"_DE_Up.xlsx"),row.names=FALSE)
write.xlsx(DE_SvsC_total,file=paste0(output,"_DE_Total.xlsx"),row.names=FALSE)

###Part 3
##topGO enrichment
library(topGO)
geneID2GO <- readMappings(go_file)
geneID2GO_uniq <- sapply(geneID2GO,function(var)unique(var))
geneList_up <- as.factor(as.integer(rownames(cts) %in%
                                             DE_SvsC_up$GeneNames))
names(geneList_up) <- rownames(cts)
geneList_down <- as.factor(as.integer(rownames(cts) %in% 
                                        DE_SvsC_down$GeneNames))
names(geneList_down) <- rownames(cts)  
go_type <- c("MF","BP","CC")

##up 构建富集GOdata数据
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

##down 构建GOdata数据
down_go <- list()
for(i in 1:length(go_type)){
  type=go_type[i]
  godata <- new("topGOdata",ontology=type,allGenes=geneList_down,
                description=paste("GOdata_SvsC_down",type,sep="\t"),annot=annFUN.gene2GO,
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
  name=paste0("Up_",go_type[i],"_Enrichment_Map")
  tmp$Annot_comb <- paste(tmp$GO.ID,tmp$Term,sep=" : ")
  tmp$qvalue <- as.numeric(tmp$qvalue)
  tmp$Significant <- as.numeric(tmp$Significant)
  tmp$Annot_comb <- factor(tmp$Annot_comb,levels = rev(tmp$Annot_comb))
  p<- ggplot(tmp,aes(qvalue,Annot_comb))+geom_point(aes(size=Significant,color=qvalue))+
    scale_color_gradient(low="red",high="green")+scale_x_reverse()+
    labs(color="Classic Fisher Qvalue",size="Significant Count",x="Classic Fisher Qvalue",
         y="GO Terms",title=name)+theme(plot.title=element_text(hjust = 0.5))+theme_bw()
  ggsave(paste0("./GO_enrichment_results/",output,"_",name,".pdf"),plot=p,width=25,height=15,units = "cm")
  printGraph(up_go[[i]],up_go_results[[i]],
                 firstSigNodes = 10,useInfo = "all",pdfSW=F,
                 fn.prefix=paste0("./GO_enrichment_results/",output,"_",name,"_DAG"))
  write.table(up_go_results_gentable[[i]],file=paste0(
    "./GO_enrichment_results/",output,"_",name,".xls"),sep="\t",quote=F,row.names = F)
} 

for(i in 1:3){
  tmp=down_go_results_table[[i]]
  name=paste0("Down_",go_type[i],"_","Enrichment_Map")
  tmp$Annot_comb <- paste(tmp$GO.ID,tmp$Term,sep=" : ")
  tmp$qvalue <- as.numeric(tmp$qvalue)
  tmp$Significant <- as.numeric(tmp$Significant)
  tmp$Annot_comb <- factor(tmp$Annot_comb,levels = rev(tmp$Annot_comb))
  p<- ggplot(tmp,aes(qvalue,Annot_comb))+geom_point(aes(size=Significant,color=qvalue))+
    scale_color_gradient(low="red",high="green")+scale_x_reverse()+
    labs(color="Classic Fisher Qvalue",size="Significant Count",x="Classic Fisher Qvalue",
         y="GO Terms",title=name)+theme(plot.title=element_text(hjust = 0.5))+theme_bw()
  ggsave(paste0("./GO_enrichment_results/",output,"_",name,".pdf"),plot=p,width=25,height=15,units = "cm")
  printGraph(down_go[[i]],down_go_results[[i]],
             firstSigNodes = 10,useInfo = "all",pdfSW=F,
             fn.prefix=paste0("./GO_enrichment_results/",output,"_",name,"_DAG"))
  write.table(down_go_results_gentable[[i]],file=paste0(
    "./GO_enrichment_results/",output,"_",name,".xls"),sep="\t",quote=F,row.names = F)
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
  
  #pdf(paste0("./GO_enrichment_results/","SvsC_up_",go_type[n],"_star_network.pdf"))
  #plot(net,layout=layout_as_star,main=paste0(
  #  "SvsC_up_",go_type[n],"_star_network"))
  #dev.off()
  
  pdf(paste0("./GO_enrichment_results/",output,"_Up_",go_type[n],"_network.pdf"))
  plot(net,layout=layout_nicely,main=paste0(
    "Up_",go_type[n],"_network"))
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
  
  #pdf(paste0("./GO_enrichment_results/","SvsC_down_",go_type[n],"_star_network.pdf"))
  #plot(net,layout=layout_as_star,main=paste0(
  #  "SvsC_down_",go_type[n],"_star_network"))
  #dev.off()
  
  pdf(paste0("./GO_enrichment_results/",output,"_Down_",go_type[n],"_network.pdf"))
  plot(net,layout=layout_nicely,main=paste0(
    "Down_",go_type[n],"_network"))
  image.plot(legend.only = TRUE,zlim=range(rescale(-log10(vertices$qvalue))),
             col=rbPal(100)[cut(seq(0,1,by=0.001),breaks=100)],
             horizontal = TRUE,legend.shrink=0.2,legend.width = 1,
             legend.args=list(text=expression(-log[10](Qvalue)),cex=0.7,line=0.1),
             axis.args=list(at=c(0,0.5,1),labels=c(0,0.5,1),line=0.05),
             smallplot = c(0.8,0.9,0.85,0.9))
  dev.off()
}

##up heatmap of GO terms with genes sorted by names
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
    labs(title=paste0("Up_",go_type[m],"_Heatmap"),y="GO Terms",x="DE Genes",fill="Qvalue")
  ggsave(paste0("./GO_enrichment_results/",output,"_Up_",go_type[m],"_by_Name_Heatmap.pdf"),
         plot=p,width = 28,height=18,units = "cm")
}

##up heatmap of GO terms with genes sorted by fold
##modify the input DE object manually
for(m in 1:3){
  tmp <- data.frame()
  genes_up <- vector()
  Data <- data.frame()
  tmp=up_go_results_table[[m]]
  for(i in 1:nrow(tmp)){
    genes_up <- append(genes_up, unlist(strsplit(tmp[i,]$Sig_Genes,",")))
  }
  genes_up <- sort(unique(genes_up))
  genes_up <- genes_up[order(DE_SvsC_up[genes_up,"log2.Fold_change."],
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
    scale_fill_gradient(low="red",high="blue")+
    scale_y_continuous(breaks=seq(1,length(tmp$GO.ID)),labels=tmp$GO.ID,expand = c(0,0))+
    scale_x_continuous(breaks=seq(1,length(genes_up)),labels=genes_up,expand = c(0,0))+
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1,size=6.5),
          plot.title=element_text(hjust = 0.5))+
    labs(title=paste0("Up_",go_type[m],"_Heatmap"),y="GO Terms",x="DE Genes",fill="Qvalue")
  ggsave(paste0("./GO_enrichment_results/",output,"_Up_",go_type[m],"_by_Fold_Heatmap.pdf"),
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
    scale_fill_gradient(low="red",high="blue")+
    scale_y_continuous(breaks=seq(1,length(tmp$GO.ID)),labels=tmp$GO.ID,expand = c(0,0))+
    scale_x_continuous(breaks=seq(1,length(genes_down)),labels=genes_down,expand = c(0,0))+
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1,size=6.5),
          plot.title=element_text(hjust = 0.5))+
    labs(title=paste0("Down_",go_type[m],"_Heatmap"),y="GO Terms",x="DE Genes",fill="Qvalue")
  ggsave(paste0("./GO_enrichment_results/",output,"_Down_",go_type[m],"_by_Name_Heatmap.pdf"),
         plot=p,width = 28,height=18,units = "cm")
}

##down heatmap of GO terms with genes by fold
##modify the input DE object manually
for(m in 1:3){
  tmp <- data.frame()
  genes_down <- vector()
  Data <- data.frame()
  tmp=down_go_results_table[[m]]
  for(i in 1:nrow(tmp)){
    genes_down <- append(genes_down, unlist(strsplit(tmp[i,]$Sig_Genes,",")))
  }
  genes_down <- sort(unique(genes_down))
  genes_down <- genes_down[order(DE_SvsC_down[genes_down,"log2.Fold_change."])]
  
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
    labs(title=paste0("Down_",go_type[m],"_Heatmap"),y="GO Terms",x="DE Genes",fill="Qvalue")
  ggsave(paste0("./GO_enrichment_results/",output,"_Down_",go_type[m],"_by_Fold_Heatmap.pdf"),
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
kegg_SvsC_up <- enrichMKEGG(gene=DE_SvsC_up$GeneNames,organism = "kpm",
                            pvalueCutoff = 1,minGSSize= 1,qvalueCutoff = 1)
kegg_SvsC_down <- enrichMKEGG(gene=DE_SvsC_down$GeneNames,organism = "kpm",
                             pvalueCutoff = 1,minGSSize = 1,qvalueCutoff = 1)

##dotplot pic / 
kegg_SvsC_up <- as.data.frame(kegg_SvsC_up)
kegg_SvsC_down <- as.data.frame(kegg_SvsC_down)
kegg_SvsC_up$Terms <- paste0(kegg_SvsC_up$ID,":",kegg_SvsC_up$Description)
kegg_SvsC_down$Terms <- paste0(kegg_SvsC_down$ID,":",kegg_SvsC_down$Description)
kegg_SvsC_results <- list(kegg_SvsC_up,kegg_SvsC_down)
name_kegg <- c("KEGG_SvsC_Up","KEGG_SvsC_Down")

for(i in 1:2){
  tmp <- data.frame()
  name <- c()
  if(nrow(kegg_SvsC_results[[i]])>30){
    tmp <- kegg_SvsC_results[[i]][1:30,]
    }else{
      tmp <- kegg_SvsC_results[[i]]
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

##heatmap for kegg enrichment sorted by names
for(m in 1:2){
  tmp <- data.frame()
  genes <- vector()
  Data <- data.frame()
  if(nrow(kegg_SvsC_results[[m]])>30){
    tmp <- kegg_SvsC_results[[m]][1:30,]
  }else{
    tmp <- kegg_SvsC_results[[m]]
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
  ggsave(paste0("./GO_enrichment_results/",name_kegg[m],"_by_Name_Heapmap.pdf"),
         plot=p,width = 28,height=18,units = "cm") 
}

##heatmap for kegg enrichment sorted by folds
for(m in 1:2){
  tmp <- data.frame()
  genes <- vector()
  Data <- data.frame()
  if(nrow(kegg_SvsC_results[[m]])>30){
    tmp <- kegg_SvsC_results[[m]][1:30,]
  }else{
    tmp <- kegg_SvsC_results[[m]]
  }
  
  for(i in 1:nrow(tmp)){
    genes <- append(genes, unlist(strsplit(tmp$geneID,"/")))
  }
  
  genes <- sort(unique(genes))
  if(m==1){
    genes <- genes[order(DE_SvsC_up[genes,"log2.Fold_change."],
                            decreasing = TRUE)]
  }else{
    genes <- genes[order(DE_SvsC_up[genes,"log2.Fold_change."])]
  }
  
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
  ggsave(paste0("./GO_enrichment_results/",name_kegg[m],"_by_Fold_Heapmap.pdf"),
         plot=p,width = 28,height=18,units = "cm") 
}

##browseKEGG
##browseKEGG(kk, 'hsa04110')
##library("pathview")
#hsa04110 <- pathview(gene.data  = geneList,
#                     pathway.id = "hsa04110",
#                     species    = "hsa",
#                     limit      = list(gene=max(abs(geneList)), cpd=1))
