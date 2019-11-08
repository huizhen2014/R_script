##ggrepel
##https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html
library(ggrepel)
set.seed(42)
dat <- subset(mtcars, wt > 2.75 & wt < 3.45)
dat$car <- rownames(dat)

p <- ggplot(dat, aes(wt, mpg, label=car)) + geom_point(color="red")
#p1 <- p + geom_text(aes(label=car),position=position_jitter(width=2,height=2))
p1 <- p + geom_text(aes(label=car))
p2 <- p + geom_text_repel(aes(label=car))
#p1 <- p + geom_text_repel(aes(label=car),size=4,vjust=1.5,family="Times",fontface="italic")
gridExtra::grid.arrange(p1,p2,ncol=2)

##hide some of the labels
set.seed(42)
dat2 <- subset(mtcars, wt > 3 & wt < 4)
dat2$car <- ""
ix_label <- c(2,3,16)
dat2$car[ix_label] <- rownames(dat2)[ix_label]

ggplot(dat2, aes(wt,mpg,label=car))+
  geom_point(color=ifelse(dat2$car == "", "grey50","red"))+
  geom_text_repel()

##limit labels to a specific area
set.seed(42)
x_limits <- c(3, NA)
##arrow参数：first指定线条末端为箭头(first,end,both)；
##type指定箭头是否为closed或open的三角；
##length指定箭头长度，angle指定箭头的角度大小
ggplot(dat, aes(wt, mpg, label=car, color=factor(cyl)))+
  geom_vline(xintercept = x_limits, linetype=3)+
  geom_point()+
  geom_label_repel(
    arrow=arrow(length=unit(0.03, 'npc'), type="closed",ends="first"),
    force=10,xlim=x_limits)+scale_color_discrete(name="cyl")

##align labels on the top or bottom edge
set.seed(42)
ggplot(mtcars, aes(x=wt, y=1, label=rownames(mtcars)))+
  geom_point(color="red")+
  geom_text_repel(
    nudge_y = 0.05,direction = "x",angle=90,vjust=0.1,segment.size=0.2)+
  xlim(1,6)+ylim(1,0.5)+
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()
  )

##align labels on the left or right edge
set.seed(42)
p <- ggplot(mtcars, aes(y=wt,x=1,label=rownames(mtcars)))+
  geom_point(color="red")+ylim(1,5.5)+theme(
    axis.line.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    plot.title=element_text(hjust=0.5)
  )

p1 <- p + xlim(1,1.375)+geom_text_repel(
  nudge_x = 0.15,
  direction = "y",
  hjust=0,
  segment.size=0.2
)+ggtitle("hjust=0")

p2 <- p + xlim(1,1.375)+geom_text_repel(
  nudge_x = 0.2,
  direction="y",
  hjust=0.5,
  segment.size=0.2
)+ggtitle("hjust=0.5(default)")

p3 <- p + xlim(0.25,1)+geom_text_repel(
    nudge_x=-0.35,
    direction="y",
    hjust=1,
    segment.size=0.2
  )+ggtitle("hjust=1")

gridExtra::grid.arrange(p1,p2,p3,ncol=3)

##align text horizontally with nudge_x and hjust
set.seed(42)
dat <- subset(mtcars,wt>2.75 & wt < 3.45)
dat$car <- rownames(dat)

ggplot(dat, aes(wt, mpg,label=car))+
  geom_text_repel(
    data = subset(dat, wt > 3),
    nudge_x = 3.5 - subset(dat, wt > 3)$wt,
    segment.size=0.2,
    segment.color="grey50",
    direction="y",
    hjust=0
  )+
  geom_text_repel(
    data = subset(dat, wt <3),
    nudge_x = 2.7 - subset(dat, wt < 3)$wt,
    segment.size = 0.2,
    segment.color="grey50",
    direction="y",
    hjust=1
  )+
  scale_x_continuous(
    breaks=c(2.5,2.75,3,3.25,3.5),
    limits=c(2.4,3.8)
  )+geom_point(color="red")

##mathematical expressions
d <- data.frame(x=c(1,2,2,1.75,1.25),
                y=c(1,3,1,2.65,1.23),
                math=c(
                  NA,
                  "integral(f(x), x %->% 0)",
                  NA,
                  "lim(f(x), x %->%0)",
                  NA
                ))
ggplot(d, aes(x,y,label=math))+
  geom_point()+
    geom_label_repel(
      parse=TRUE, ##解析数学表达式
      size=8,
      box.padding=2
    )

##https://github.com/wilkox/ggfittext
##https://zhuanlan.zhihu.com/p/83092046

library(ggplot2)
library(ggfittext)
library(gggenes)
library(ggrepel)
library(RColorBrewer)

kp_ampR <- read.delim("kp_ampR.genes",header=F)
colnames(kp_ampR) <- colnames(example_genes)

#kp_ampR$start_new <- ifelse(kp_ampR$strand == "forward",kp_ampR$start,kp_ampR$end)
#kp_ampR$end_new <- ifelse(kp_ampR$strand == "forward",kp_ampR$end,kp_ampR$start)

##输入文件存在问题，无法区分direction
kp_ampR <- within(kp_ampR,{
  direction[strand == "reverse"] <- "-1"
  direction[strand == "forward"] <- "1"
})
kp_ampR$direction <- as.numeric(kp_ampR$direction)

##构建调色板函数getPalette
getPalette <- colorRampPalette(brewer.pal(9,"Set1"))
colourCount <- length(unique(kp_ampR$gene))

##geom_gene_label min.size指定最小的图形内的字体大小
ggplot(kp_ampR,aes(xmin=start,xmax=end,y=molecule,
                   fill=gene,label=gene,forward=direction))+
  geom_gene_arrow(arrowhead_height = unit(4,"mm"),arrowhead_width = unit(2,"mm"))+
  facet_wrap(~molecule,scales="free",ncol=1) + 
  scale_fill_manual(values=getPalette(colourCount))+
  geom_text_repel(aes(x=start+(end-start)/2),size=3,direction = "x",
                  nudge_y=0.3,fontface="italic")+
  labs(title="AmpR_Proximity_Gene_Distributuion",x="Gene Location",y="Sample")+
  theme(axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        plot.title = element_text(hjust = 0.5)
        )

