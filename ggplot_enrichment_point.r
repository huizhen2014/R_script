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
##绘制
ggplot(GOdata_21vs28_l_MF_results,aes(weight01,Term))+
  geom_point(aes(size=Significant,color=-log10(as.numeric(as.character(weight01)))))+
  scale_color_gradient(low="green",high="red")+
  labs(color=expression("-log[10](Weight Pvalue)"),size="Significant Count",x="Weight01 Pvalue",
       y="Pathway",title="GOdata_21vs28_l_MF Pathway Enrichement")


