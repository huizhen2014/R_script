##Plot the vocalno pic with Different expression genes matrix
library(ggplot2)
library(ggrepel)
library(xlsx)

##Get the Expression table
ab_Exp <- read.xlsx("ab_c507_0440_vs_ab_c_DE_total_3vs2.xlsx",1,header=T,
                    colClasses = "character")
##Modify MS_Exp by adding Product iterms
#ab_Exp$Product <- substr(MS_Exp$Protein.description,1,
#                         regexpr("OS=",MS_Exp$Protein.description)-1)

##Assign the color,trim genes with the padj == NA 
ab_Exp <- ab_Exp[!is.na(ab_Exp$padj),]
ab_Exp <- within(ab_Exp,{
  Expression <- NA ##必须指定, Not Available
  Expression[(log2FoldChange < 1)|(log2FoldChange > -1)] <- "NA"
  Expression[(log2FoldChange >= 1) & (padj < 0.05)] <- "Up"
  Expression[(-1 >= log2FoldChange)& (padj < 0.05)] <- "Down"
}
)

##添加名称列，可对应显示名称信息,对应将需要显示的名称添加到新列Name，其他为NA即可
##geom_text(aes(label=Name), size=4)
##aes尽量在需要用时才指定，ggplot里指定影响面太大
colours <- c("Up"="red","Down"="green","NA"="grey")
#ggplot(MS_Exp_1.2,aes(log(ampR_1.ampR.Ratio,2),-log(ampR_1.ampR.P.value,10)))+
#  geom_point(aes(color=Expression))+
#  scale_color_manual(values=colours,labels=c("High","Low","Unsig"))+
#  theme_bw()+ ylim(0,3.5)+xlim(-6,6)+
#  geom_label_repel(aes(log(ampR_1.ampR.Ratio,2),-log(ampR_1.ampR.P.value,10),
#                       label=ifelse(Expression != "NA",Product,"")),size=3)+
#  labs(x=expression(log[2](Fold)),y=expression(-log[10](Pvalue)))+
#  geom_hline(yintercept =-log10(0.05),color="black",linetype="dotted")+
#  geom_vline(xintercept = range(log(1/1.2,2),log(1.2,2)),color=c("green","red"),linetype="dotted")

ggplot(ab_Exp,aes(log2FoldChange,-log(padj,10)))+
  geom_point(aes(color=Expression))+ xlim(-4,4)+
  scale_color_manual(values=colours,breaks=c("Up","Down","NA"),labels=c("Up","Down","No Significant"))+
  theme_bw()+
  labs(x=expression(Log[2]("Fold Change")),y=expression(-Log[10]("p-Value")),color="Type")+ ## 修改图例标签
  geom_hline(yintercept =-log10(0.05),color="black",linetype="dotted")+
  geom_vline(xintercept = range(-1,log(2,2)),color=c("green","red"),linetype="dotted")


##Save the pic
ggsave("ab_c507_0440_vs_ab_c_vocalno.pdf",width=7,height = 7.5)

