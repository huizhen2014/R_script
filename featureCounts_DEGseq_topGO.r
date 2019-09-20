##align with Rsubread
#align(index="HS11286_index",readfile1 = "./kp_21_sangon_1P.fastq.gz",readfile2 = "./kp_21_sango    n_2P.fastq.gz",
#      output_file = "kp_21_sangon_align.bam",type="rna")
#align(index="HS11286_index",readfile1 = "./kp_28_sangon_1P.fastq.gz",readfile2 = "./kp_28_sango    n_2P.fastq.gz",
#      output_file = "kp_28_sangon_align.bam",type="rna")
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
##bowtie2-build GCF_000240185.1_ASM24018v2_genomic.fna hs11286_index
##bowtie2 -x hs11286_index -1 kp_21_sangon_1P.fastq.gz -2 kp_21_sangon_2P.fastq.gz -S kp_21_sangon_mapped.sam
##bowtie2 -x hs11286_index -1 kp_28_sangon_1P.fastq.gz -2 kp_28_sangon_2P.fastq.gz -S kp_28_sangon_mapped.sam
##samtools view -bS kp_21_sangon_mapped.sam > kp_21_sangon_mapped.bam
##samtools view -bS kp_28_sangon_mapped.sam > kp_28_sangon_mapped.bam

####featureCounts ignoreDup=T
##featureCounts
kp_sangon_21_exp <- featureCounts("kp_21_sangon_mapped.bam",annot.ext = "HS11286.gtf",
                                  isGTFAnnotationFile = T,GTF.featureType = "transcript",
                                  GTF.attrType = "Name",isPairedEnd = T,ignoreDup=T,
                                  requireBothEndsMapped=T,nthreads=2,countChimericFragments=F)
kp_sangon_28_exp <- featureCounts("kp_28_sangon_mapped.bam",annot.ext = "HS11286.gtf",
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
normMat <- data.frame(kp_21=fc_sangon_21$Length,kp_28=fc_sangon_28[rownames(fc_sangon_21),]$Length,
                      row.names=rownames(fc_sangon_21))

##edgeR
library(edgeR)
group=factor(c("kp_21","kp_28"),levels=c("kp_21","kp_28"))
y <- DGEList(cts,group=group)

##filter count 0
y <- y[!rowSums(y$counts == 0 )>0,]
y$samples$lib.size=colSums(y$counts)
fc_sangon_21 <- fc_sangon_21[fc_sangon_21$Count !=0,]
fc_sangon_28 <- fc_sangon_28[fc_sangon_28$Count !=0,]
                     
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
fc_sangon_28$PRKM <- rpkm_sangon[rownames(fc_sangon_28),]$kp_28
fc_sangon_21$TPM <- ((fc_sangon_21$Count/fc_sangon_21$Length)*1000000)/
  (sum(fc_sangon_21$Count/fc_sangon_21$Length))
fc_sangon_28$TPM <- ((fc_sangon_28$Count/fc_sangon_28$Length)*1000000)/
  (sum(fc_sangon_28$Count/fc_sangon_28$Length))

fc_sangon_21$Norm.Count <- fc_sangon_21$Count/y$samples$norm.factors[1]
fc_sangon_28$Norm.Count <- fc_sangon_28$Count/y$samples$norm.factors[2]

### cts count normalized by norm.factors and is used for DEGseq package analysis
library(DEGseq)
cts_norm <- data.frame(name=rownames(cts),kp_21_cts=cts$kp_21/y$samples$norm.factors[1],
                       kp_28_cts=cts$kp_28/y$samples$norm.factors[2])
cts_norm <- cts_norm[cts_norm$kp_21_cts!=0,]
cts_norm <- cts_norm[cts_norm$kp_28_cts!=0,]
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

##tmpe test 
##DE_28vs21_h <- read.table(file="/Data_analysis/RNA_analysis/K.p_analysis/sangon_analysis_results/result/7_DEA/DEG/up_1",
##                          sep="\t",header=F,row.names = 1)
##DE_28vs21_l <- read.table(file="/Data_analysis/RNA_analysis/K.p_analysis/sangon_analysis_results/result/7_DEA/DEG/down_1",
##                          sep="\t",header=F,row.names = 1)
###

library(topGO)
geneID2GO <- readMappings("KPHS_Genbank_Go.txt")
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
                gene2GO=geneID2GO,nodeSize=5)
  ##renew the GOdata
  .geneList_28vs21_h <- as.factor(as.integer(genes(godata) %in% sigGenes(godata)))
  names(.geneList_28vs21_h) <- genes(godata)
  godata <- new("topGOdata", ontology=type,allGenes=.geneList_28vs21_h,
                description=paste("GOdata_28vs21_h",type,sep="\t"),annot=annFUN.gene2GO,
                gene2GO=geneID2GO,nodeSize=5)
  kp_28vs21_h_go[[i]] <- godata
}

statistic <- c("classic","weight01","elim")
kp_28vs21_h_go_results <- list()
kp_28vs21_h_go_results_gentable <- list()
for(i in 1:length(kp_28vs21_h_go)){
  godata <- kp_28vs21_h_go[[i]]
  tmp=list()
  for(j in 1:length(statistic)){
    s <- statistic[j]
    result <- runTest(godata,algorithm = s,statistic = "fisher")
    tmp[[j]] <- result
  }
  kp_28vs21_h_go_results_gentable[[i]] <- GenTable(godata,classic=tmp[[1]],weight01=tmp[[2]],
                                                   elim=tmp[[3]],orderBy="classic",ranksOf="classic",
                                                   topNodes=20)
  kp_28vs21_h_go_results[[i]] <- tmp
}

##28 vs 21 low
geneList_28vs21_l <- as.factor(as.integer(rownames(cts) %in% 
                                            DE_28vs21_l$GeneNames))
names(geneList_28vs21_l) <- rownames(cts)  
kp_28vs21_l_go <- list()
for(i in 1:length(go_type)){
  type=go_type[i]
  godata <- new("topGOdata",ontology=type,allGenes=geneList_28vs21_l,
                description=paste("GOdata_28vs21_h",type,sep="\t"),annot=annFUN.gene2GO,
                gene2GO=geneID2GO,nodeSize=5)
  ##renew the genelist
  .geneList_28vs21_l <- as.factor(as.integer(genes(godata) %in% sigGenes(godata)))
  names(.geneList_28vs21_l) <- genes(godata)
  godata <- new("topGOdata",ontology=type,allGenes=.geneList_28vs21_l,
                description=paste("GOdata_28vs21_h",type,sep="\t"),annot=annFUN.gene2GO,
                gene2GO=geneID2GO,nodeSize=5)
  kp_28vs21_l_go[[i]] <- godata
}
kp_28vs21_l_go_results <- list()
kp_28vs21_l_go_results_gentable <- list()                
for(i in 1:length(kp_28vs21_l_go)){
  godata <- kp_28vs21_l_go[[i]]
  tmp=list()
  for(j in 1:length(statistic)){
    s <- statistic[j]
    result <- runTest(godata,algorithm = s,statistic = "fisher")
    tmp[[j]] <- result
  }
  kp_28vs21_l_go_results_gentable[[i]] <- GenTable(godata,classic=tmp[[1]],weight01=tmp[[2]],
                                                   elim=tmp[[3]],orderBy="classic",ranksOf="classic",
                                                   topNodes=20)
  kp_28vs21_l_go_results[[i]] <- tmp
}

