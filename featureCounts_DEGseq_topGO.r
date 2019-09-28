##trimmomatic
#trimmomatic PE -summary kp_21_trimmed.log kp_21_R1.fastq.gz kp_21_R2.fastq.gz -baseout kp_21_trimmed.fastq.gz \
#ILLUMINACLIP:Sabgon_adapter.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:35 
##sangon_adapter.fa
##>PrefixPE/1
##AGATCGGAAGAGCACACGTCTGAAC
##>PrefixPE/2
##AGATCGGAAGAGCGTCGTGTAGGGA

##align with Rsubread
library(Rsubread)
#buildindex(basename="HS11286_index",reference="GCF_000240185.1_ASM24018v2_genomic.fna")
align(index="HS11286_index",readfile1 = "./kp_21_trimmed.pe_1P.fastq.gz",
      readfile2 = "./kp_21_trimmed.pe_2P.fastq.gz",minFragLength=35,
      output_file = "kp_21_align.bam",type="rna")
align(index="HS11286_index",readfile1 = "./kp_28_trimmed.pe_1P.fastq.gz",
      readfile2 = "./kp_28_trimmed.pe_2P.fastq.gz",minFragLength=35,
      output_file = "kp_28_align.bam",type="rna")
repair("kp_21_align.bam",inFormat="BAM",outFiles="kp_21_align_repair.bam",
       addDummy=T,compress=T,nthreads=2)
repair("kp_28_align.bam",inFormat="BAM",outFiles="kp_28_align_repair.bam",
       addDummy=T,compress=T,nthreads=2)

##featureCounts
#kp_sangon_21_exp <- featureCounts("kp_21_sangon_align.bam",annot.ext = "HS11286.gtf",
#                                  isGTFAnnotationFile = T,GTF.featureType = "transcript",
#                                  GTF.attrType = "Name",isPairedEnd = T,ignoreDup=T,
#                                  requireBothEndsMapped=T,nthreads=2,countChimericFragments=F)
#kp_sangon_28_exp <- featureCounts("kp_28_sangon_align.bam",annot.ext = "HS11286.gtf",
#                                  isGTFAnnotationFile = T,GTF.featureType = "transcript",
#                                  GTF.attrType = "Name",isPairedEnd = T,ignoreDup=T,
#                                  requireBothEndsMapped=T,nthreads=2,countChimericFragments=F)

###align bowtie2 with shell scripts
#bowtie2-build GCF_000240185.1_ASM24018v2_genomic.fna hs11286_index
#bowtie2 -x hs11286_index -1 kp_21_sangon_1P.fastq.gz -2 kp_21_sangon_2P.fastq.gz -S kp_21_sangon_mapped.sam
#bowtie2 -x hs11286_index -1 kp_28_sangon_1P.fastq.gz -2 kp_28_sangon_2P.fastq.gz -S kp_28_sangon_mapped.sam
#samtools view -bS kp_21_sangon_mapped.sam > kp_21_sangon_mapped.bam
#samtools view -bS kp_28_sangon_mapped.sam > kp_28_sangon_mapped.bam
#samtools sort kp_21_sangon_mapped.bam -o kp_21_sangon_mapped_sorted.bam
#samtools sort kp_28_sangon_mapped.bam -o kp_28_sangon_mapped_sorted.bam

####featureCounts ignoreDup=T
##featureCounts
kp_sangon_21_exp <- featureCounts("kp_21_align_repair.bam",annot.ext = "HS11286.gtf",
                                  isGTFAnnotationFile = T,GTF.featureType = "transcript",
                                  GTF.attrType = "Name",isPairedEnd = T,ignoreDup=T,
                                  requireBothEndsMapped=T,nthreads=2,countChimericFragments=F)
kp_sangon_28_exp <- featureCounts("kp_28_align_repair.bam",annot.ext = "HS11286.gtf",
                                  isGTFAnnotationFile = T,GTF.featureType = "transcript",
                                  GTF.attrType = "Name",isPairedEnd = T,ignoreDup=T,
                                  requireBothEndsMapped=T,nthreads=2,countChimericFragments=F)

###collect counts ignoreDup=T
fc_sangon_21 <- data.frame(Gene_ID=rownames(kp_sangon_21_exp$counts),Pos=paste0(
  kp_sangon_21_exp$annotation$Chr,"[",kp_sangon_21_exp$annotation$Strand,"]",
  kp_sangon_21_exp$annotation$Start,"-",kp_sangon_21_exp$annotation$End),
  Length=kp_sangon_21_exp$annotation$Length,Count=kp_sangon_21_exp$counts[,1])
fc_sangon_28 <- data.frame(Gene_ID=rownames(kp_sangon_28_exp$counts),Pos=paste0(
  kp_sangon_28_exp$annotation$Chr,"[",kp_sangon_28_exp$annotation$Strand,"]",
  kp_sangon_28_exp$annotation$Start,"-",kp_sangon_28_exp$annotation$End),
  Length=kp_sangon_28_exp$annotation$Length,Count=kp_sangon_28_exp$counts[,1])

##prepare edgeR input
cts <- data.frame(kp_21=fc_sangon_21$Count,kp_28=fc_sangon_28[rownames(fc_sangon_21),]$Count,
                  row.names=rownames(fc_sangon_21))

##edgeR
library(edgeR)
group=factor(c("kp_21","kp_28"),levels=c("kp_21","kp_28"))
y <- DGEList(cts,group=group)

##filter count 0/5
#y <- y[!rowSums(y$counts == 0 )>0,]
y <- y[rowSums(y$counts >=5)==2,]
y$samples$lib.size=colSums(y$counts)

##normalization TMM
y <- calcNormFactors(y,method="TMM")

##updated fc_sangon cpm rpkm tpm normalized count
cpm_sangon <- data.frame(cpm(y))
fc_sangon_21$CPM <- cpm_sangon[rownames(fc_sangon_21),]$kp_21
fc_sangon_28$CPM <- cpm_sangon[rownames(fc_sangon_28),]$kp_28
#length
rpkm_sangon <- data.frame(rpkm(y,gene.length =
  c(fc_sangon_21[rownames(y$counts),]$Length,
    fc_sangon_28[rownames(y$counts),]$Length)))
fc_sangon_21$RPKM <- rpkm_sangon[rownames(fc_sangon_21),]$kp_21
fc_sangon_28$RPKM <- rpkm_sangon[rownames(fc_sangon_28),]$kp_28
##这里使用的是fc_sangon dataframe中的值,y中对应值为NA
fc_sangon_21$TPM <- ((fc_sangon_21$Count/fc_sangon_21$Length)*1000000)/
  (sum(fc_sangon_21$Count/fc_sangon_21$Length))
fc_sangon_28$TPM <- ((fc_sangon_28$Count/fc_sangon_28$Length)*1000000)/
  (sum(fc_sangon_28$Count/fc_sangon_28$Length))

fc_sangon_21$Norm.Count <- fc_sangon_21$Count/y$samples$norm.factors[1]
fc_sangon_28$Norm.Count <- fc_sangon_28$Count/y$samples$norm.factors[2]

### cts count normalized by norm.factors and is used for DEGseq package analysis
library(DEGseq)
cts_norm <- data.frame(row.names=rownames(cts),
                       name=rownames(cts),kp_21_cts=cts$kp_21/y$samples$norm.factors[1],
                       kp_28_cts=cts$kp_28/y$samples$norm.factors[2])

##根据y过滤掉低表达genes，即两样本基因均满足count>=5
cts_norm <- cts_norm[rownames(y$counts),]

write.table(cts_norm,file="cts_normalized_by_norm_factor.txt",
            sep="\t",row.names = F)
kp_21_cts <- readGeneExp(file="cts_normalized_by_norm_factor.txt",geneCol=1,valCol = 2)
kp_28_cts <- readGeneExp(file="cts_normalized_by_norm_factor.txt",geneCol=1,valCol = 3)
layout(matrix(c(1,2,3,4,5,6),3,2,byrow=T))
par(mar=c(2,2,2,2))
DEGexp(geneExpMatrix1 = kp_28_cts,geneCol1 = 1,expCol1 = 2,groupLabel1 = "kp_28_cts",
       geneExpMatrix2 = kp_21_cts,geneCol2 = 1,expCol2 = 2,groupLabel2 = "kp_21_cts",
       method="MARS",rawCount = F,thresholdKind=3,qValue=0.05,
       outputDir="./kp_sangon_cts_norm_degseq_results")
###the result is same to Rsubread map and counts results

##topGO analysis 
kp_score <- read.delim("./kp_sangon_cts_norm_degseq_results/output_score.txt",
                      header=T,sep="\t")
DE_28vs21_h <- kp_score[kp_score$log2.Fold_change. > 1 & 
                          kp_score$Signature.q.value.Benjamini.et.al..1995....0.05.=="TRUE",]
DE_28vs21_l <- kp_score[kp_score$log2.Fold_change. < -1 &
                          kp_score$Signature.q.value.Benjamini.et.al..1995....0.05.=="TRUE",]

##topGO enrichment
library(topGO)
geneID2GO <- readMappings("HS11286_sangon_go.txt")
##28 vs 21 high
geneList_28vs21_h <- as.factor(as.integer(rownames(cts) %in%
                                             DE_28vs21_h$GeneNames))
names(geneList_28vs21_h) <- rownames(cts)
go_type <- c("MF","BP","CC")
kp_28vs21_h_go <- list()
for(i in 1:length(go_type)){
  type=go_type[i]
  godata <- new("topGOdata",ontology=type,allGenes=geneList_28vs21_h,
                description=paste("GOdata_28vs21_h",type,sep="\t"),annot=annFUN.gene2GO,
                gene2GO=geneID2GO,nodeSize=1)
  ##renew the GOdata
  .geneList_28vs21_h <- as.factor(as.integer(genes(godata) %in% sigGenes(godata)))
  names(.geneList_28vs21_h) <- genes(godata)
  godata <- new("topGOdata", ontology=type,allGenes=.geneList_28vs21_h,
                description=paste("GOdata_28vs21_h",type,sep="\t"),annot=annFUN.gene2GO,
                gene2GO=geneID2GO,nodeSize=1)
  kp_28vs21_h_go[[i]] <- godata
}

#statistic <- c("classic","weight01","elim")
kp_28vs21_h_go_results <- list()
kp_28vs21_h_go_results_table <- list()
kp_28vs21_h_go_results_gentable <- list()
for(i in 1:length(kp_28vs21_h_go)){
  godata <- kp_28vs21_h_go[[i]]
  for(j in 1:length(statistic)){
    result <- runTest(godata,algorithm = "classic",statistic = "fisher")
  }
  kp_28vs21_h_go_results_gentable[[i]] <- GenTable(godata,classic=result,orderBy="classic",ranksOf="classic",
                                                   topNodes=length(usedGO(godata)),numChar=50)
  kp_28vs21_h_go_results_gentable[[i]]$qvalue <- p.adjust(sort(score(result)),"BH")
  kp_28vs21_h_go_results_gentable[[i]] <- kp_28vs21_h_go_results_gentable[[i]][order(kp_28vs21_h_go_results_gentable[[i]]$qvalue),]
  kp_28vs21_h_go_results_table[[i]] <- kp_28vs21_h_go_results_gentable[[i]][1:30,]
  kp_28vs21_h_go_results_gentable[[i]]$Term <- Definition(kp_28vs21_h_go_results_gentable[[i]]$GO.ID)
  
  kp_28vs21_h_go_results_gentable[[i]]$Sig_Genes <- sapply(
    sapply(genesInTerm(godata,kp_28vs21_h_go_results_gentable[[i]]$GO.ID),function(x){
      sigGenes(godata)[sigGenes(godata) %in% x]
    }),function(y){paste(y,collapse = ",")})
  
  kp_28vs21_h_go_results_gentable[[i]]$All_Genes <- sapply(
    genesInTerm(godata,kp_28vs21_h_go_results_gentable[[i]]$GO.ID),
    function(x){paste(x,collapse = ",")})
  
  kp_28vs21_h_go_results[[i]] <- result
}

##28 vs 21 low
geneList_28vs21_l <- as.factor(as.integer(rownames(cts) %in% 
                                            DE_28vs21_l$GeneNames))
names(geneList_28vs21_l) <- rownames(cts)  
kp_28vs21_l_go <- list()
for(i in 1:length(go_type)){
  type=go_type[i]
  godata <- new("topGOdata",ontology=type,allGenes=geneList_28vs21_l,
                description=paste("GOdata_28vs21_l",type,sep="\t"),annot=annFUN.gene2GO,
                gene2GO=geneID2GO,nodeSize=1)
  ##renew the genelist
  .geneList_28vs21_l <- as.factor(as.integer(genes(godata) %in% sigGenes(godata)))
  names(.geneList_28vs21_l) <- genes(godata)
  godata <- new("topGOdata",ontology=type,allGenes=.geneList_28vs21_l,
                description=paste("GOdata_28vs21_l",type,sep="\t"),annot=annFUN.gene2GO,
                gene2GO=geneID2GO,nodeSize=1)
  kp_28vs21_l_go[[i]] <- godata
}

kp_28vs21_l_go_results <- list()
kp_28vs21_l_go_results_table <- list()
kp_28vs21_l_go_results_gentable <- list()                
for(i in 1:length(kp_28vs21_l_go)){
  godata <- kp_28vs21_l_go[[i]]
  for(j in 1:length(statistic)){
    result <- runTest(godata,algorithm = "classic",statistic = "fisher")
    }
  kp_28vs21_l_go_results_gentable[[i]] <- GenTable(godata,classic=result,orderBy="classic",ranksOf="classic",
                                                   topNodes=length(usedGO(godata)),numChar=50)
  kp_28vs21_l_go_results_gentable[[i]]$qvalue <- p.adjust(sort(score(result)),"BH")
  kp_28vs21_l_go_results_gentable[[i]] <- kp_28vs21_l_go_results_gentable[[i]][order(kp_28vs21_h_go_results_gentable[[i]]$qvalue),]
  kp_28vs21_l_go_results_table[[i]] <- kp_28vs21_l_go_results_gentable[[i]][1:30,]
  kp_28vs21_l_go_results_gentable[[i]]$Term <- Definition(kp_28vs21_l_go_results_gentable[[i]]$GO.ID)
  
  kp_28vs21_l_go_results_gentable[[i]]$Sig_Genes <- sapply(
    sapply(genesInTerm(godata,kp_28vs21_l_go_results_gentable[[i]]$GO.ID),function(x){
      sigGenes(godata)[sigGenes(godata) %in% x]
      }),function(y){paste(y,collapse = ",")})
    
  kp_28vs21_l_go_results_gentable[[i]]$All_Genes <- sapply(
    genesInTerm(godata,kp_28vs21_l_go_results_gentable[[i]]$GO.ID),
    function(x){paste(x,collapse = ",")})

  kp_28vs21_l_go_results[[i]] <- result
}

##绘制GO 富集散点图,ggsave根据后缀判断类型(device)
library(ggplot2)
dir.create("./GO_enrichment_results")
for(i in 1:3){
  tmp=kp_28vs21_h_go_results_table[[i]]
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
  printGraph(kp_28vs21_h_go[[i]],kp_28vs21_h_go_results[[i]],
                 firstSigNodes = 10,useInfo = "all",
                 fn.prefix=paste0("./GO_enrichment_results/",name,"DAG"))
  write.table(kp_28vs21_h_go_results_gentable[[i]],file=paste0(
    "./GO_enrichment_results/",name,".xls"),sep="\t",quote=F,row.names = F)
} 

for(i in 1:3){
  tmp=kp_28vs21_l_go_results_table[[i]]
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
  printGraph(kp_28vs21_l_go[[i]],kp_28vs21_l_go_results[[i]],
             firstSigNodes = 10,useInfo = "all",
             fn.prefix=paste0("./GO_enrichment_results/",name,"DAG"))
  write.table(kp_28vs21_l_go_results_gentable[[i]],file=paste0(
    "./GO_enrichment_results/",name,".xls"),sep="\t",quote=F,row.names = F)
} 

##clusterProfiler kegg analysis
library(clusterProfiler)
##kegg_database <- search_kegg_organism("Klebsiella pneumoniae",by="scientific_name",
##                                      ignore.case = T)
##kpm
kegg_28vs21_m_h <- enrichMKEGG(gene=DE_28vs21_h$GeneNames,organism = "kpm",
                            pvalueCutoff = 0.05)
kegg_28vs21_m_l <- enrichMKEGG(gene=DE_28vs21_l$GeneNames,organism = "kpm",
                             pvalueCutoff = 0.05)
kegg_28vs21_h <- enrichKEGG(gene=DE_28vs21_h$GeneNames,organism = "kpm",
                               pvalueCutoff = 0.05)
kegg_28vs21_l <- enrichKEGG(gene=DE_28vs21_l$GeneNames,organism = "kpm",
                               pvalueCutoff = 0.05)

kp <- loadDb("/Users/carlos/.AnnotationHub/74476/org.Klebsiella_pneumoniae_subsp._pneumoniae_HS11286.eg.sqlite")
DE_down_dataframe <- select(kp,keys=as.vector(DE_down),columns=c("ENTREZID","SYMBOL"),keytype="SYMBOL")
ego_down_mf <- enrichGO(DE_down_dataframe$ENTREZID,OrgDb = kp,keyType = "ENTREZID",
                        ont="MF",readable = T)






  





