##绘制topGO富集气泡图
##输入data.frame为Gentable返回数据框,x=weight01,y=Term
##allRes <- GenTable(sampleGOdata, calssicFisher=resultFisher, classicKS=resultKS, elimKS=resultKS.elim, orderBy="elimKS", ranksOf="classicFisher", topNodes=10)
##对Term定义factor水平
GOdata_21vs28_l_MF_results$Term <- factor(GOdata_21vs28_l_MF_results$Term,
                                          levels=rev(GOdata_21vs28_l_MF_results$Term))
##对Pvalue定义factor水平
GOdata_21vs28_l_MF_results$weight01 <- factor(GOdata_21vs28_l_MF_results$weight01,levels=
                                               sort(unique(GOdata_21vs28_l_MF_results$weight01),
                                                    decreasing = T))
##恢复factor回原来的数字
##GOdata_21vs28_l_MF_results$weight01 <- as.numeric(as.character(GOdata_21vs28_l_MF_results$weight01))
##将Significant值转换为数字
GOdata_21vs28_l_MF_results$Significant <- as.numeric(GOdata_21vs28_l_MF_results$Significant)
##绘制，此时weight01为因子，所有x坐标骤不连续
ggplot(GOdata_21vs28_l_MF_results,aes(weight01,Term))+
  geom_point(aes(size=Significant,color=-log10(as.numeric(as.character(weight01)))))+
  scale_color_gradient(low="green",high="red")+
  labs(color=expression(-log[10]("Weight Pvalue")),size="Significant Count",
  x="Weight01 Pvalue",y="Pathway",title="GOdata_21vs28_l_MF Pathway Enrichement")+
  theme(plot.title=element_text(hjust=0.5),axis.text.x=element_text(angle=30, vjust=0.5))
##取消weight01因子，x坐标连续
 ggplot(GOdata_21vs28_h_BP_results,aes(rev(as.numeric(as.character(GOdata_21vs28_h_BP_results$weight01))),Term))+
 geom_point(aes(size=Significant,color=-log10(as.numeric(as.character(weight01)))))+
 scale_color_gradient(low="green",high="red")+
 labs(x="Rev Weight01 Fisher Pvalue",y="Pathway",size="Significant Count",color=expression(-log[10]("weight01 pvalue")),title="GOdata_21vs28_h_BP_weight01_fisher_results")+
 theme(plot.title=element_text(hjust=0.5))
##保存结果
ggsave("GOdata_21vs28_h_BP_weight01_fisher_results.pdf",width=28,height=15,units="cm")


##2019年9月26日绘制
library(ggplot2)
tmp=GenTable_results[[1]]
tmp$Annot_comb <- paste(tmp$GO.ID,tmp$Term,sep=" : ")
tmp$classic <- as.numeric(tmp$classic)
tmp$Significant <- as.numeric(tmp$Significant)
tmp$Annot_comb <- factor(tmp$Annot_comb,levels = rev(tmp$Annot_comb))
tmp$value_log <- -log10(tmp$classic)
tmp$value_log_norm <- (tmp$value_log - min(tmp$value_log))/(max(tmp$value_log)-min(tmp$value_log))
##未采用归一至0-1的方式，直接根据p值分配颜色
ggplot(tmp,aes(classic,Annot_comb))+geom_point(aes(size=Significant,color=classic))+
  scale_color_gradient(low="red",high="green")+scale_x_reverse()+
  labs(color="Classic Fisher Pvalue",size="Significant Count",x="Classic Fisher Pvalue",
       y="GO Terms",title="GO Enrichment MF")+theme(plot.title=element_text(hjust=0.5))+theme_bw()



