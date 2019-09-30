##for gwas significant points plot
library(ggplot2)
hyper_file  <- read.csv("hypermortality_25_09_2019_1530.results.csv")
selected_gwas_genes <- hyper_file[is.na(as.integer(sapply(hyper_file$Gene,function(x){grep("_",x)}))),][,c("Gene","Benjamini_H_p")]
selected_gwas_genes.bak <- selected_gwas_genes
selected_gwas_genes$Gene <- factor(selected_gwas_genes$Gene,levels=rev(selected_gwas_genes$Gene))

#selected_gwas_genes$Italic_gene <- expression(paste(selected_gwas_genes$Gene,
#                                                    italic('italics')))
p <- ggplot(selected_gwas_genes,aes(-log10(Benjamini_H_p),Gene))+geom_point(aes(color=-log10(Benjamini_H_p)))+
  scale_color_gradient(low="green",high="red") + theme_bw()+
  theme(axis.text.y=element_text(size=6,hjust=0.5,face="italic"),
        plot.margin=unit(c(0.5,5,0.5,5),"cm"),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15))+
  labs(y="Gene",x=expression(-log[10](Qvalue)),color=expression(-log[10](Qvalue)))
p 
ggsave(p,file="hypermortality_ampR.pdf",width=9.5,height=7)
