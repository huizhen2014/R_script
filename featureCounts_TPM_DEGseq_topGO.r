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
#y <- y[rowSums(y$counts >=5)==2,]
#y$samples$lib.size=colSums(y$counts)

##normalization TMM
y <- calcNormFactors(y,method="TMM")

##updated fc_sangon cpm rpkm tpm normalized count, 这些值计算均为原始count值，
##不是均一化后counts值
cpm_sangon <- data.frame(cpm(y))
fc_sangon_21$CPM <- cpm_sangon[rownames(fc_sangon_21),]$kp_21
fc_sangon_28$CPM <- cpm_sangon[rownames(fc_sangon_28),]$kp_28
#length
rpkm_sangon <- data.frame(rpkm(y,gene.length =
  c(fc_sangon_21[rownames(y$counts),]$Length,
    fc_sangon_28[rownames(y$counts),]$Length)))
fc_sangon_21$RPKM <- rpkm_sangon[rownames(fc_sangon_21),]$kp_21
fc_sangon_28$RPKM <- rpkm_sangon[rownames(fc_sangon_28),]$kp_28

fc_sangon_21$TPM <- ((fc_sangon_21$Count/fc_sangon_21$Length)*1000000)/
  (sum(fc_sangon_21$Count/fc_sangon_21$Length))
fc_sangon_28$TPM <- ((fc_sangon_28$Count/fc_sangon_28$Length)*1000000)/
  (sum(fc_sangon_28$Count/fc_sangon_28$Length))

fc_sangon_21$Norm.Count <- fc_sangon_21$Count/y$samples$norm.factors[1]
fc_sangon_28$Norm.Count <- fc_sangon_28$Count/y$samples$norm.factors[2]

cts_norm <- data.frame(row.names=rownames(cts),
                       name=rownames(cts),kp_21_cts=cts$kp_21/y$samples$norm.factors[1],
                       kp_28_cts=cts$kp_28/y$samples$norm.factors[2])
cts_norm$TPM_21 <- fc_sangon_21[rownames(cts_norm),]$TPM
cts_norm$TPM_28 <- fc_sangon_28[rownames(cts_norm),]$TPM
#cts_norm_filtered <- cts_norm[(cts_norm$TPM_21 >=5 |cts_norm$TPM_28 >=5),]

### cts count normalized by norm.factors and is used for DEGseq package analysis
library(DEGseq)
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

DE_28vs21_h <- kp_score[cts_norm[(cts_norm$TPM_21 >=5 & cts_norm$TPM_28 >=5),]$name,][kp_score$log2.Fold_change. > 1 & 
                          kp_score$Signature.q.value.Benjamini.et.al..1995....0.05.=="TRUE",]
DE_28vs21_l <- kp_score[cts_norm[(cts_norm$TPM_21 >=5 & cts_norm$TPM_28 >=5),]$name,][kp_score$log2.Fold_change. < -1 &
                          kp_score$Signature.q.value.Benjamini.et.al..1995....0.05.=="TRUE",]

##topGO enrichment
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
                gene2GO=geneID2GO,nodeSize=2)
  ##renew the GOdata
  .geneList_28vs21_h <- as.factor(as.integer(genes(godata) %in% sigGenes(godata)))
  names(.geneList_28vs21_h) <- genes(godata)
  godata <- new("topGOdata", ontology=type,allGenes=.geneList_28vs21_h,
                description=paste("GOdata_28vs21_h",type,sep="\t"),annot=annFUN.gene2GO,
                gene2GO=geneID2GO,nodeSize=2)
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
                                                   topNodes=30)
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
                gene2GO=geneID2GO,nodeSize=2)
  ##renew the genelist
  .geneList_28vs21_l <- as.factor(as.integer(genes(godata) %in% sigGenes(godata)))
  names(.geneList_28vs21_l) <- genes(godata)
  godata <- new("topGOdata",ontology=type,allGenes=.geneList_28vs21_l,
                description=paste("GOdata_28vs21_h",type,sep="\t"),annot=annFUN.gene2GO,
                gene2GO=geneID2GO,nodeSize=2)
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
                                                   topNodes=30)
  kp_28vs21_l_go_results[[i]] <- tmp
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

