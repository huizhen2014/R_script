##for gwas significant points plot

tmp8 <- read.csv("hypermortality_25_09_2019_1530.results.csv")
selected_gwas_genes <- tmp8[is.na(as.integer(sapply(tmp8$Gene,function(x){grep("_",x)}))),][,c("Gene","Benjamini_H_p")]
selected_gwas_genes.bak <- selected_gwas_genes
selected_gwas_genes$Gene <- factor(selected_gwas_genes$Gene,levels=as.vector(selected_gwas_genes$Gene))

ggplot(selected_gwas_genes,aes(Gene,-log10(Benjamini_H_p)))+geom_point(aes(color=-log10(Benjamini_H_p)))+
  scale_color_gradient(low="green",high="red") + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5,size=4),plot.title=element_text(hjust=0.5))+
  labs(x="Gene",y="-log10(qvalue)",color="-log10(qvalue)")+ggtitle("Plot")+theme_bw()

