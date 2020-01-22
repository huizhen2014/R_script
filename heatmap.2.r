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

hm <- as.matrix(hm_data)
colnames(hm) <- c("c507_0440_1","c507_0440_2","c507_0440_3","control","control")

##对每行进行scale
hm_tran <- t(scale(t(x)))

##绘图
pdf("ab_c507_0440_vs_ab_c_heatmap.pdf",width=7,height=8)
heatmap.2(hm_tran,col="greenred",dendrogram="both",tracecol = "transparent",
          margins=c(7,7),lhei=c(1,4.5),lwid=c(1.5,4.5),key.title="Color Key",
          key.xlab=NA,key.ylab = NA,cexCol=1,srtCol=30,
          #cexCol = 1,cexRow = 0.25, 默认根据实际情况调整，实际可自调xy轴字体大小以显示全部信息
         keysize=1 )
dev.off()

