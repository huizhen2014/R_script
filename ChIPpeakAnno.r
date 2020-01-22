##Bacteria ChIPpeakAnno 
library(ChIPpeakAnno)
library(GenomicFeatures)
library(IRanges)
library(xlsx)
library(BSgenome.H10407.NCBI.01)

anno_file="h10407_annotation_extraction.txt"
peaks_file="H_NS_peaks.xls"
output="hns_macs2"

anno_extraction <- read.delim(anno_file,sep="\t",
                              header=F)
anno_extraction$V6 <- ifelse(grepl("complement",anno_extraction$V2),"-",
                             "+")
tmp <- gsub("[a-z><=()]*","",anno_extraction$V2)
tmp[grep(",",tmp,fixed=T)] <- 
  sapply(strsplit(tmp[grep(",",tmp,fixed=T)],",",fixed=T),
       function(x)x[[2]])

tmp_start <- sapply(strsplit(tmp,"..",fixed=T),function(x)x[[1]])
tmp_end <- sapply(strsplit(tmp,"..",fixed=T),function(x)x[[2]])
anno_extraction$V7 <- as.numeric(tmp_start)
anno_extraction$V8 <- as.numeric(tmp_end)
row.names(anno_extraction) <- anno_extraction$V3

##construct annotation GRanges file
##locus_tag changed to feature to adapt the annotatePeakInBatch function
anno_extraction_gr <- GRanges(seqnames = anno_extraction$V1,
                              ranges=IRanges(start = as.numeric(tmp_start),
                                             end=as.numeric(tmp_end)),
                              strand=anno_extraction$V6,
                              data.frame(feature=anno_extraction$V3,
                                               protein_id=anno_extraction$V4,
                                               product=anno_extraction$V5))
names(anno_extraction_gr) <- anno_extraction$V3

peaks <- toGRanges(peaks_file,format="MACS2")
peaks_anno <- annotatePeakInBatch(peaks,
                                      AnnotationData = anno_extraction_gr,
                                      maxgap = 1000L)
peaks_anno_strand <- peaks_anno
strand(peaks_anno_strand) <- peaks_anno_strand$feature_strand

peaks_anno_strand_seq <- getAllPeakSequence(peaks_anno_strand,
                                                genome = BSgenome.H10407.NCBI.01,
                                                upstream = 50L,downstream = 50L)
peaks_anno_strand_dataframe <- data.frame(peaks_anno_strand)
peaks_anno_strand_dataframe$protein_id <- anno_extraction[
  peaks_anno_strand_dataframe$feature,"V4"]
peaks_anno_strand_dataframe$product <- anno_extraction[
  peaks_anno_strand_dataframe$feature,"V5"]

write2FASTA(peaks_anno_strand_seq,paste0(output,"_anno_strand_seq.fa"))
write.table(peaks_anno_strand_dataframe,paste0(output,"_anno_strand.txt"),
            col.names=T,row.names=FALSE,sep="\t",quote=FALSE)
write.xlsx(peaks_anno_strand_dataframe,paste0(output,"_anno_strand_dataframe.xlsx"),
           row.names = FALSE,col.names = TRUE)





