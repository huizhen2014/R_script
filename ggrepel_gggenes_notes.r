##ggrepel
##https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html
##https://github.com/wilkox/ggfittext
##https://zhuanlan.zhihu.com/p/83092046

library(ggplot2)
library(ggfittext)
library(gggenes)
library(ggrepel)
library(RColorBrewer)
library(export)

kl <- read.delim("KL1_47_64.summary",header=F,stringsAsFactors = F)
kl_genes <- data.frame(molecule=kl$V4,gene=ifelse(kl$V6 != "",kl$V6,"hypothetical"),
                       start=as.numeric(sapply(strsplit(kl$V5,".",fixed=T),function(x)x[1])),
                       end=as.numeric(sapply(strsplit(kl$V5,".",fixed=T),function(x)x[3])),
                       direction="1")

#k1 <- data.frame(molecule=kl[1:20,]$V4,gene=ifelse(kl[1:20,]$V6 != "",kl[1:20,]$V6,"hypothetical"),
#                       start=as.numeric(sapply(strsplit(kl[1:20,]$V5,".",fixed=T),function(x)x[1])),
#                       end=as.numeric(sapply(strsplit(kl[1:20,]$V5,".",fixed=T),function(x)x[3])),
#                       direction="1")

#k47 <- data.frame(molecule=kl[21:39,]$V4,gene=ifelse(kl[21:39,]$V6 != "",kl[21:39,]$V6,"hypothetical"),
#                 start=as.numeric(sapply(strsplit(kl[21:39,]$V5,".",fixed=T),function(x)x[1])),
#                 end=as.numeric(sapply(strsplit(kl[21:39,]$V5,".",fixed=T),function(x)x[3])),
#                 direction="1")

#k64 <- data.frame(molecule=kl[40:63,]$V4,gene=ifelse(kl[40:63,]$V6 != "",kl[40:63,]$V6,"hypothetical"),
#                  start=as.numeric(sapply(strsplit(kl[40:63,]$V5,".",fixed=T),function(x)x[1])),
#                  end=as.numeric(sapply(strsplit(kl[40:63,]$V5,".",fixed=T),function(x)x[3])),
#                  direction="1")

dummies <- make_alignment_dummies(
  kl_genes,
  aes(xmin=start,xmax=end,y=molecule,id=gene),
  on="galF"
)

colorCount <- length(unique(kl_genes$gene))
getPalette <- colorRampPalette(brewer.pal(12,"Set3"))
ggplot(kl_genes,aes(xmin=start,xmax=end,y=molecule,fill=gene))+
  geom_gene_arrow()+facet_wrap(~molecule,scales="free_y",ncol=1)+
  scale_fill_manual(values=getPalette(colorCount))+ geom_blank(data=dummies)+
  geom_text_repel(aes(x=start+(end-start)/2,label=gene),size=3,direction = "x",
                  nudge_y=0.15,nudge_x=0.3,fontface="italic",segment.size = 0,angle=30)+
                    theme_genes()+theme(legend.position = "none")
#ggsave("kl_genes.pdf",width = 10,height=5)
graph2ppt(file="kl_genes.pptx")















######
kp_ampR <- read.delim("kp_ampR.genes",header=F)
colnames(kp_ampR) <- colnames(example_genes)

#kp_ampR$start_new <- ifelse(kp_ampR$strand == "forward",kp_ampR$start,kp_ampR$end)
#kp_ampR$end_new <- ifelse(kp_ampR$strand == "forward",kp_ampR$end,kp_ampR$start)

##输入文件存在问题，无法区分direction
kp_ampR <- within(kp_ampR,{
  direction[strand == "reverse"] <- "-1"
  direction[strand == "forward"] <- "1"
})
kp_ampR$direction <- as.numeric(kp_ampR$direction)

##构建调色板函数getPalette
getPalette <- colorRampPalette(brewer.pal(9,"Set1"))
colourCount <- length(unique(kp_ampR$gene))

##geom_gene_label min.size指定最小的图形内的字体大小
ggplot(kp_ampR,aes(xmin=start,xmax=end,y=molecule,
                   fill=gene,label=gene,forward=direction))+
  geom_gene_arrow(arrowhead_height = unit(4,"mm"),arrowhead_width = unit(2,"mm"))+
  facet_wrap(~molecule,scales="free",ncol=1) + 
  scale_fill_manual(values=getPalette(colourCount))+
  geom_text_repel(aes(x=start+(end-start)/2),size=3,direction = "x",
                  nudge_y=0.3,fontface="italic")+
  labs(title="AmpR_Proximity_Gene_Distributuion",x="Gene Location",y="Sample")+
  theme(axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        plot.title = element_text(hjust = 0.5)
        )

