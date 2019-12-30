##
##
library(ggplot2)
data <- data.frame()
#go_type <- c("MF","BP","CC")
for(i in 1:3){
  down_go_results_table[[i]]$Type <- go_type[i]
  down_go_results_table[[i]]$Name <- reorder(
    paste0(down_go_results_table[[i]]$GO.ID," : ",down_go_results_table[[i]]$Term),
    -down_go_results_table[[i]]$qvalue) ##排序，并映射到ggplot的坐标顺序中
  data <- rbind(data, down_go_results_table[[i]][down_go_results_table[[i]]$qvalue < 0.05,])
}

##绘制geom_col() + coord_flip()
ggplot(data,aes(Name,Significant))+geom_col(aes(fill=Type))+ coord_flip()+
  theme(axis.text.x=element_text(size=10,face="bold"))
ggsave("ac_11_go_down_col.pdf",width=9.5,height=7)
dev.off()

