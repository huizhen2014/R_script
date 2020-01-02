##plot the MIC Inhibition plot
library(ggplot2)
library(xlsx)
library(reshape2)
library(scales)
## read.xlsx file
mic_e <- read.xlsx("MIC_efficacy.xlsx",1)
colnames(mic_e) <- sub("NA.","Compound",sub("X","",colnames(mic_e)))
mic_e_melt <- melt(mic_e,id.vars="Compound",variable.name="Tigecycline",
value.name="Inhibition")

mic_e_melt$Tigecycline <-
  round(as.numeric(levels(mic_e_melt$Tigecycline)[mic_e_melt$Tigecycline]),digits = 3)
mic_e_melt$Compound <-
  round(mic_e_melt$Compound,digits = 3)

mic_e_melt$Compound <- factor(mic_e_melt$Compound,
  levels=sort(unique(mic_e_melt$Compound)))
mic_e_melt$Tigecycline <- factor(mic_e_melt$Tigecycline,
  levels=sort(unique(mic_e_melt$Tigecycline)))

## ggplot geom_tile()
ggplot(mic_e_melt,aes(x=Tigecycline,y=Compound,fill=Inhibition))+geom_tile()+
  scale_fill_continuous(low="blue",high="white",
                        breaks=c(0,0.5,1),
                        labels=c(0,50,100))+
  theme(axis.text.y=element_text(size=10,face="bold"),
        axis.text.x=element_text(size=10,face="bold"),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15))+
  labs(x=expression("Vancomycin ("*mu*"g/ml)"),
       y=expression("XT-1 ("*mu*"g/ml)"),fill="Inhibition(%)")

#ggsave("HQ_1_vancomycin_FICI.pdf",width=8,height=6.5)

#pdf("CD102_5G_16_coumpunds_FICI_2.pdf",width=11,height=8.5)
#grid.arrange(p16,p17,p18,p20,p22,p27,ncol=2)
#dev.off()

#pdf("CD102_5G_16_coumpunds_FICI_1.pdf",width=11,height=8.5)
#grid.arrange(p1,p2,p5,p6,p7,p8,p10,p15,ncol=2)
#dev.off()
