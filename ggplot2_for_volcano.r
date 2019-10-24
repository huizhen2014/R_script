library(ggplot2)
library(ggrepel)
##get the table of edgeR $count
table_43vsLAC_4 <- lrt_0.04_43vsLAC_4$table
##distribute the color according to the logFC and PValue from table_43vsLAC_4
New_table_43vsLAC_4 <- within(table_43vsLAC_4,{
  color <- NA
  color[logFC < -1 & PValue < 0.01 ] <- "green"
  color[logFC > 1 & PValue < 0.01 ] <- "red"
  color[logFC > -1 & logFC < 1 ] <- "black"
  color[PValue > 0.01 ] <- "grey"
})
##添加名称列，可对应显示名称信息,对应将需要显示的名称添加到新列Name，其他为NA即可
##geom_text(aes(label=Name), size=4)
colours <- c(red="red",green="green",black="black",grey="grey")
ggplot(New_table_43vsLAC_4,aes(logFC,-log10(PValue)))+geom_point(aes(color=color))+
  scale_color_manual(values=colours)+labs(title="43 vs LAC_4",x="logFC",y="-log10(PValue")+
  geom_hline(yintercept = range(-log10(0.01),-log10(0.001)),color=c("blue","red"),linetype="dotted")+
  geom_vline(xintercept=range(-1,1),color="yellow",linetype="solid")
##save the pic
ggsave("43_vs_LAC_4_valcono.png",with=5,height=5)

##更新版，同时修改图例显示内容
#lrt_bcv_0.2_table_new_new <- within(lrt_bcv_0.2_table_new,{
# color_new <- NA
# color_new[logFC > 1 & PValue < 0.05 ] <- "logFC>1&PValue<0.05"
# color_new[logFC < -1 & PValue < 0.05 ] <- "logFC<-1&PValue<0.05"
# color_new[logFC > -1 & logFC < 1 ] <- "-1<logFC<1"
# color_new[PValue > 0.05 ]<- "PValue>0.05"
# })
#colours <- c("-1<logFC<1"="black","logFC<-1&PValue<0.05"="green","PValue>0.05"="grey","logFC>1&PValue<0.05"="red")
#ggplot(lrt_bcv_0.2_table_new_new,aes(logFC,-log10(PValue)))+geom_point(aes(color=color_new))+scale_color_manual(values=colours)
