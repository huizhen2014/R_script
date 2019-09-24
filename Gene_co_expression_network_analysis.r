##Gene Co-expresion Network Analysis 
##unweighted networks using hard thresholding 
##weighted networks using soft thresholding

##nodes 代表基因，若在恰当的组织样本中对应的基因共表达，则这些nodes会连接起来
##adjacency matrix/相邻矩阵 networks; A=[aij] 为对称矩阵，A值位于[0,1],传统而言，对角线值为0
##对于unweighted networkds,相邻矩阵包含binary信息,相连为1，不相连为0;
##在weighted networks中，相邻矩阵包含weights.

##该tutuorial data来自yeast microarray data set
##预测酵母生存必须基因,若这些基因紧密相连，则倾向于为必须基因，为1，否则为0
##简而言之，基因的显著性为任何定量检测判断该基因的生物显著性程度，network分析的
##一个目的就是关联基因显著性(必要程度)的测量值和内部模型的联系.

##模型构建
##将共表达的一组基因构建成一模型，一般使用average linkage hierachical clustering
##该方法使用topological overlap测量计算dissimilarity
##两个nodes的topological overlap反应了它们连接节点范围的相似性(the topological
##overlap of two nodes reflects their similarity in terms to the commonality
##of the nodes they connect to)
##当通过hierarchical clustering method获得dendrogram时，我们需要选择一高度阈值
##来实现clustering, 该阈值用于cut the branches.
##高度阈值(height cutoff)的选择可通过检测dendrogram中的一些branches对应
##离散度digonal blocks(modules) in the TOM plot. 一般推荐使用cutreeDynamic
##用于brance cuting
##这里采用静态(fixed) cutree method


library(WGCNA)
options(stringsAsFactors=F)
dat0 <- read.csv("YEASTCellCycle4000.csv",header=T,row.names=1)
##gene exprssion data: 行为array(样本),列为基因
datExpr <- t(dat0)
##数据信息，其中essentiality列指明该基因是否对于酵母生存必须
datSummary <- read.csv("YEASTCellCycle4000.csv",header=T,row.names=1)

##限制分析范围为变化最大的基因
var1 <- function(x)var(x,na.rm=T)
vargenes <- apply(datExpr,2,var1)
rankvar <- rank(-vargenes)
restVariance <- rankvar < 4001 ##这里保留所有基因
datExpr <- datExpr[,restVariance]
datSummary <- datSummary[restVariance,]

##使用scale-free topology creterion 来选择cutoff值
##使用线性回归模型fitting index来评估network满足
##scale-free topology的程度

##针对不同的soft thresoholds(powers) 构建scale free
##topology fitting indices, 用于挑选soft threshold
powers1 <- c(seq(1,10,by=1),seq(12,20,by=2))

cex1=1
par(mfrow=c(1,2))
plot(RpowerTable[,1],-sign(RpowerTable[,3])*RpowerTable[,2],
xlab="Soft Threshold(power)", ylab="Scale Free Topology Model Fit, signed R^2",type="n")
text(RpowerTable[,1],-sign(RpowerTable[,3]*RpowerTable[,2],
labels=powers1,cex=cex1,col="red")
#This line correspond to using an R^2 cutoff of h
abline(h=0.85,col="red")
plot(RpowerTable[,1],RpowerTable[,5],xlab="Soft Threshold(power)",
ylab="Mean Connectivity", type="n")
text(RpowerTable[,1],RpowerTable[,5],labels=powers1,cex=cex1,col="red")

##at power=7, the 拟合曲线的R^2值出现拐点，也就是随着power的增加，
##scale free topology不在提高，因此选择beta1=7
##该过程通过最大化sacle free topology model fit(R^2),维持了高的
##平均连接数目，同时参数值导致R^2接近1也可能带来具有非常少连接的networks
##由于生物学上而言，相比于没有hub genes，networks不太可能拥有多个hubs gene
##we multiply R^2 with -1 if the slop of the regression line between 
##log_{10}(p(k))and log_{10}(k)is positive

##使用该power用于power adjacency function
beta1 < -7
Connectivity <- softConnectivity(datExpr,power=beta1)
##构建scale free topology，黑曲线对应scale free topology
##红曲线对应 truncated scale free topology
par(mfrow=c(1,1))
scaleFreePlot(Connectivity,main=paste("soft threshold, power",beta1),truncated=T)

##检出模型
##network 分析最重要的一个步骤就是module detection
##这里使用clustering in combination with the topological overlap matrix






