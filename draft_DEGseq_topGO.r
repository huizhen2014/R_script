##trimmomatic
trimmomatic PE -summary kp_21_trimmed.log kp_21_R1.fastq.gz kp_21_R2.fastq.gz -baseout kp_21_trimmed.fastq.gz \
ILLUMINACLIP:Sangon_adapter.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:35 
##sangon_adapter.fa
##>PrefixPE/1
##AGATCGGAAGAGCACACGTCTGAAC
##>PrefixPE/2
##AGATCGGAAGAGCGTCGTGTAGGGA


###align bowtie2 with shell scripts /featureCounts from Rsubread package: ignoreDup=T 
library(Rsubread)
bowtie2-build GCF_000240185.1_ASM24018v2_genomic.fna hs11286_index
bowtie2 -x hs11286_index -1 kp_21_sangon_1P.fastq.gz -2 kp_21_sangon_2P.fastq.gz -S kp_21_sangon_mapped.sam
bowtie2 -x hs11286_index -1 kp_28_sangon_1P.fastq.gz -2 kp_28_sangon_2P.fastq.gz -S kp_28_sangon_mapped.sam
samtools view -bS kp_21_sangon_mapped.sam > kp_21_sangon_mapped.bam
samtools view -bS kp_28_sangon_mapped.sam > kp_28_sangon_mapped.bam

kp_sangon_21_exp <- featureCounts("kp_21_sangon_mapped.bam",annot.ext = "HS11286.gtf",
                                  isGTFAnnotationFile = T,GTF.featureType = "transcript",
                                  GTF.attrType = "Name",isPairedEnd = T,ignoreDup=T,
                                  requireBothEndsMapped=T,nthreads=2,countChimericFragments=F)
kp_sangon_28_exp <- featureCounts("kp_28_sangon_mapped.bam",annot.ext = "HS11286.gtf",
                                  isGTFAnnotationFile = T,GTF.featureType = "transcript",
                                  GTF.attrType = "Name",isPairedEnd = T,ignoreDup=T,
                                  requireBothEndsMapped=T,nthreads=2,countChimericFragments=F)

##edgeR prepare edgeR input
cts <- data.frame(kp_21=kp_sangon_21_exp$Count,kp_28=kp_sangon_28_exp[rownames(kp_sangon_21_exp),]$Count,
                  row.names=rownames(fc_sangon_21))
library(edgeR)
group=factor(c("kp_21","kp_28"),levels=c("kp_21","kp_28"))
y <- DGEList(cts,group=group)
##filter count 0
y <- y[!rowSums(y$counts == 0 )>0,]
y$samples$lib.size=colSums(y$counts)
##normalization TMM
y <- calcNormFactors(y,method="TMM")
cts$Norm.Count_21 <- cts$kp_21/y$samples$norm.factors[1]
cts$Norm.Count_28 <- cts$kp_28/y$samples$norm.factors[2]

### cts count normalized by norm.factors and is used for DEGseq package analysis
library(DEGseq)
write.table(cts,file="cts_normalized_by_norm_factor.txt",
            sep="\t",row.names = F)
kp_21_cts <- readGeneExp(file="cts_normalized_by_norm_factor.txt",geneCol=1,valCol = 3)
kp_28_cts <- readGeneExp(file="cts_normalized_by_norm_factor.txt",geneCol=1,valCol = 4)
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

###topGO enrichment analysis
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
