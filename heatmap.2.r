##
##gplots heatmap.2
library(pheatmap)
library(gplots)
##绘制差异基因热图, log2FoldChange > 1/< -1  & padj < 0.05
hm_data <- res[(res$log2FoldChange > 1 | res$log2FoldChange< -1) & 
                      ! is.na(res$padj) & res$padj < 0.05,]
##排序
hm_data_order <- order(hm_data$log2FoldChange,decreasing = FALSE)
hm_data <- assay(dds)[rownames(hm_data),]
hm_data <- hm_data[hm_data_order,]

x <- as.matrix(hm_data)
colnames(x) <- c("cs_11","cs_11","control","control")

##对每行进行scale
x_tran <- t(scale(t(x)))

##绘图
pdf("ab_cs_11_vs_ab_c_heatmap.pdf",width=7,height=8)
heatmap.2(x_tran,col="greenred",dendrogram="both",scale="row",tracecol = "transparent",
          margins=c(7,7),lhei=c(1,4.5),lwid=c(1.5,4.5),key.title="Color Key",
          
          key.xlab=NA,key.ylab = NA,
         keysize=1 )
dev.off()
