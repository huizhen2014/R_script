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
ggsave("GOdata_21vs28_h_BP_weight01_fisher_results.pdf")
