##Plot the vocalno pic with Different expression genes matrix
library(ggplot2)
library(ggrepel)
library(xlsx)

##Get the Expression table
MS_Exp <- read.xlsx("MS_identified_information.xlsx",2,header=T,
                    colClasses = "character")
##Modify MS_Exp by adding Product iterms
MS_Exp$Product <- substr(MS_Exp$Protein.description,1,
                         regexpr("OS=",MS_Exp$Protein.description)-1)

##Assign the color
MS_Exp_1.2 <- within(MS_Exp,{
  Expression <- NA ##必须指定, Not Available
  Expression[(ampR_1.ampR.Ratio <= 1.2)|(ampR_1.ampR.Ratio >= 1/1.2)] <- "NA"
  Expression[(ampR_1.ampR.Ratio > 1.2) & (ampR_1.ampR.P.value < 0.05)] <- "High"
  Expression[(ampR_1.ampR.Ratio < 1/1.2) & (ampR_1.ampR.P.value < 0.05)] <- "Low"
}
)

##添加名称列，可对应显示名称信息,对应将需要显示的名称添加到新列Name，其他为NA即可
##geom_text(aes(label=Name), size=4)
colours <- c("High"="red","Low"="green","NA"="grey")
ggplot(MS_Exp_1.2,aes(log(ampR_1.ampR.Ratio,2),-log(ampR_1.ampR.P.value,10)))+
  geom_point(aes(color=Expression))+
  scale_color_manual(values=colours,labels=c("High","Low","Unsig"))+
  theme_bw()+ ylim(0,3.5)+xlim(-6,6)+
  geom_label_repel(aes(log(ampR_1.ampR.Ratio,2),-log(ampR_1.ampR.P.value,10),
                       label=ifelse(Expression != "NA",Product,"")),size=3)+
  labs(x=expression(log[2](Fold)),y=expression(-log[10](Pvalue)))+
  geom_hline(yintercept =-log10(0.05),color="black",linetype="dotted")+
  geom_vline(xintercept = range(log(1/1.2,2),log(1.2,2)),color=c("green","red"),linetype="dotted")
##Save the pic
ggsave("MS_Exp_1.2.pdf",width=6,height = 7.5)
