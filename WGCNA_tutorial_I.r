##I. Network analysis of liver expression data in female mice
##1a Data input and cleaning
library(WGCNA)
options(stringsAsFactors = FALSE)
femData <- read.csv("LiverFemale3600.csv")
##去除辅助信息
datExpr0 <- as.data.frame(t(femData[,-c(1:8)]))
names(atExpr0) <- femData$substanceBXH
rownames(datExpr0) <- names(femData)[-c(1:8)]

##1b checking data for excessive missing values and identification
##of outlier microarray samples
##goodSamplesGenes，函数检测表达数据，去除含有过多缺失值的基因和样本
##行为样本，列为基因
gsg <- goodSamplesGenes(datExpr0,verbose=3)
gsg$allOK
##若gsg$allOK 不为TRUE，选择下面方法过滤
if(!gsg$allOK){
  if(sum(!gsg$goodGenes)>0){
    printFlush(paste("Remove genes:",paste(names(datExpr0)[!gsg$goodGenes],
                                           collapse = ", ")))
  }
  if(sum(!gsg$goodSamples)>0){
    printFlush(paste("Removing samples:",paste(rownames(datExpr0)[!gsg$goodSamples],
                                               collapse = ", ")))
  }
  datExpr0 <- datExpr0[gsg$goodSamples,gsg$goodGenes]
}

##cluster samples to see whether there are any obvious outliers
##dist 计算matrix行之间的距离
sampleTree <- hclust(dist(datExpr0),method="average")
##plot the pic, window dimensions size 12 /9 inches
#dir.create("./Plots")
sizeGrWindow(12,9)
#pdf(file="Plots/sampleClustering.pdf",width=12,height = 9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree,main="Sample clustering to detect outlier",
     sub="",xlab="",cex.lab=1.5,cex.axis=1.5,cex.main=2)
#dev.off()
##图中可见一个离散样本点，F2_221。可手动去除，或采用自动方式,
##choose a height cut that will remove the offending sample
##绘制cut line
abline(h=15,col="red")
##确定位于line下的cluster
##cutreeStatic : A numeric vector giving labels of objects, 
##with 0 meaning unassigned. The largest cluster is conventionally labeled 1, 
##the next largest 2, etc.
clust <- cutreeStatic(sampleTree,cutHeight = 15,minSize = 10)
table(clust)
##保留clust为1的样本
keepSamples <- (clust == 1)
datExpr <- datExpr0[keepSamples,]
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
##此时的datExpr包含表达数据可用于下面网络分析

##1c loading clinical trait data
traitData <- read.csv("ClinicalTraits.csv")
dim(traitData)
names(traitData)
##去除不需要信息
allTraits <- traitData[,-c(31,16)]
allTraits <- allTraits[,c(2,11:36)]
dim(allTraits)
names(allTraits)

femaleSamples <- rownames(datExpr)
traitRows <- match(femaleSamples, allTraits$Mice)
datTraits <- allTraits[traitRows,-1]
rownames(datTraits) <- allTraits[traitRows,1]

##Performs garbage collection until free memory idicators show no change.
collectGarbage()

##现在拥有表达数据datExpr和对应临床数据datTraits
##在进行网络分析前，查看临床数据和样本系统发育树的关系
##re-cluster samples
sampleTree2 <- hclust(dist(datExpr),method="average")
##使用颜色代表特征traits：白色表低表达，红色表高表达，
##灰色表示缺失
traitColors <- numbers2colors(datTraits,signed = FALSE)
##绘图
plotDendroAndColors(sampleTree2,traitColors,
                    groupLabels = names(datTraits),
                    main="Sample dendrogram and trait heatmap")
##最后保存相关表达和特征数据
save(datExpr,datTraits,file="FemaleLiver-01-dataInput.RData")

####################################################################
##2a automatic network construction and module detection
##2a.1 Choosing the soft-threshonding power: analysis of network topology
##根据近似sacle-free topology，选择soft thresholding power
##选择一套预设的soft-thresholding powers
##datExpr 行为样本，列为基因
powers <- c(c(1:10),seq(12,20,by=2))
sft <- pickSoftThreshold(datExpr,powerVector = powers,verbose=5)
##第一列为soft thrshold Power,第二列为fitting index R^2(scale.law.R.2)
##第三列为fitting line slope，第5/6/7列为平均，中值，最大connectivity
##选择标准，R^2和power 出现转弯处，第二列R^@大于0.8，slop位于-1左右

par(mfrow=c(1,2))
cex1 <- 0.9
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (Power)",ylab="Scale Free Topolgy Model Fit, signed R^2",
     type="n",main=paste("Scale independence"))
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.9,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#########################################################################
##2a.2 one-step network construction and module detection
#net <- blockwiseModules(datExpr,power=6,
#                        TOMType = "unsigned",minModuleSize = 30,
#                        reassignThreshold = 0,mergeCutHeight = 0.25,
#                        numericLabels = TRUE,pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE,saveTOMFileBase = "femaleMouseTOM",
#                        verboase=3)
##通过net$colors查看识别的模块
#table(net$colors)
##结果0表示所有模块以外的基因数目,对应模块下表示该模块所拥有的基因数目
##针对模块构建系统发育树，并着色
#mergedColors <- labels2colors(net$colors)
#plotDendroAndColors(net$dendrograms[[1]],mergedColors[net$blockGenes[[1]]],
#                    "Module colors",
#                    dendroLabels = FALSE,hang=0.03,
#                    addGuide = TRUE,guideHang = 0.05)
##使用函数recutBlockwiseTrees重新根据修正后标准计算network并构建系统发育树
##最后保存结果
#moduleLabels <- net$colors
#moduleColors <- labels2colors(net$colors)
#MEs <- net$MEs
#geneTree <- net$dendrograms[[1]]
#save(MEs,moduleLabels,moduleColors,geneTree,
#     file="FemaleLiver-02-networkConstruction-auto.RData")
###############################################################

##2b Step-by-Step network construction and module detection
##2b.1 Choosing the soft-thresholding power: analysis of network topology
sft <- pickSoftThreshold(datExpr, powerVector=powers,verbose=5)

##2b.2 Co-expression similarity and adjacency
##使用soft thresholding power 6
softPower <- 6
##Calculates (correlation or distance) network adjacency from 
##given expression data or from a similarity.
adjacency <- adjacency(datExpr,power=softPower)

##2b.3 Topological Overlap Matrix(TOM)
##通过将adjacency转成Topological Overlap Matrix，减少
##噪音和错误关联的影响,同时计算对应的dissimilarity
##Calculation of the topological overlap matrix, 
##and the corresponding dissimilarity, from a given adjacency matrix.
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1- TOM

##2b.4 Clustering using TOM
##构建基因的hierarchical clustering tree(dendrogram)
geneTree <- hclust(as.dist(dissTOM),method="average")
plot(geneTree,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity",
     labels=FALSE,hang=0.04)
##系统发育树中密集的group表示其紧密相关，高度共表达基因

###The function cutreeStaticColor colors each gene by the branches that 
##result from choosing a particular height cut-off.
##colorh1= cutreeStaticColor(geneTree,cutHeight = 0.97)
##table(colorh1)

##使用cutreeDynamic函数识别基因模块
##power6,相对小对基因模块数量30，中等cluster splitting敏感度deepSplit=2
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree,distM = dissTOM,
                             deepSplit = 2,pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
table(dynamicMods)
##结果0表示所有模块以外的基因数目,对应模块下表示该模块所拥有的基因数目

##针对模块构建系统发育树，并着色
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree,dynamicColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE,hang=0.03,
                    addGuide=TRUE,guideHang = 0.05,
                    main="Gene dendrogram and module colors")

##2b.5 Mergeing of modules whose expression profiles are very similar
##dynamic tree cut可能识别表达形式非常相似的模型，合并这些模型会更合理
##通过计算cluster的eigengens，评估整个模型的共表达形似性并根据其相关性
##聚类
##Calculates module eigengenes (1st principal component) of modules 
##in a given single dataset.
MEList <- moduleEigengenes(datExpr,colors=dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss),method="average")
plot(METree,main="Clustering of module eigengenes",
     xlab="",sub="")

##选择hight cut 0.25，对应相关性为0.75进行合并
MEDissTres <- 0.25
abline(h=MEDissTres,col="red")
##Merges modules in gene expression networks that are too close as measured 
##by the correlation of their eigengenes.
merge <- mergeCloseModules(datExpr,dynamicColors,cutHeight = MEDissTres,
                           verbose = 3)
mergedColors <-merge$colors
mergedMEs <- merge$newMEs

##查看原始的和合并后的模型
plotDendroAndColors(geneTree,cbind(dynamicColors,mergedColors),
                    c('Dynamic Tree Cut','Merged dynamic'),
                    dendroLabels = FALSE, hang=0.03,
                    addGuide = TRUE,guideHang = 0.05)
##保存结果
moduleColors <- mergedColors
colorOrder <- c("grey",standardColors(50))
moduleLabels <- match(moduleColors,colorOrder)
MEs <- mergedMEs
save(MEs,moduleLabels,moduleColors,geneTree,file="
     FemaleLiver-02-networkConstruction-stepBystep.RData")

#####################################################################
##2c Dealing with large data sets: block-wise network construction and module detection
##2c.1 Chooseing the soft-threshonding power: analysis of network topology
##as before

##2c.2 Block-wise network construction and module detection
##先将基因分成小的block进行粗糙的聚类，然后分别对每个block进行
##全network 分析
##power6,相对小对基因模块数量30，中等cluster splitting敏感度deepSplit=2
##mergeCutHeigh表示合并相似模型高度0.25
#bwnet <- blockwiseModules(datExpr,maxBlockSize = 2000,
#                          power=6, TOMType = "unsigned",
#                          minModuleSize = 30,reassignThreshold = 0,
#                          mergeCutHeight = 0.25,
#                          numericLabels = TRUE,saveTOMs = TRUE,
#                          saveTOMFileBase = "femaleMouseTOM-blackwise",
#                          verbose=3)
##和自动方式一样，建议stepbystep进行

##比较使用2a single-block模块方式
#laod(file="FemaleLiver-02-networkConstruction-auto.RData")
#bwLabels <- matchLabels(bwnet$colors,moduleLabels)
#bwModuleColors <- labels2colors(bwLabels)
#plotDendroAndColors(bwnet$dendrograms[[1]],bwModuleColors[bwnet$blocks[[1]]],
#                    "Module colors",main="Gene dendrogram and module colors in block 1",
#                    dendroLabels = FALSE,hang=0.03,addGuide=TRUE,guideHang = 0.05)
#plotDendroAndColors(bwnet$dendrograms[[2]],bwModuleColors[bwnet$blockGenes[[2]]],
#                    "Module colors",main="Gene dendrogram and module colors in block 2",
#                    dendroLabels = FALSE, hang=0.03,
#                    addGuide=TRUE,guideHang = 0.05)
##2c.3 Comparing the single block and block-wise network analysis
#plotDendroAndColors(geneTree,
#                    cbind(moduleColors,bwModuleColors),
#                    c("Single block","2 block"),
#                    main="Single block gene dendrogram and module colors",
#                    dendroLabels = FALSE, hang=0.03,
#                    addGuide = TRUE, guideHang = 0.05)
##sigle-block and block-wise过程结果很相似
##计算eigengenes
#singleBlockMEs <- moduleEigengenes(datExpr,moduleColors)$eigengenes
#blockwiseMEs = moduleEigengenes(datExpr, bwModuleColors)$eigengenes
#single2blockwise = match(names(singleBlockMEs), names(blockwiseMEs))
#signif(diag(cor(blockwiseMEs[,single2blockwise],singleBlockMEs)), 3)
##所有相关性都接近1
#########################################################################

##3. Relating modules to external information and identifying important genes
##3a Quantifying module-trait associations
##简单地关联eigengenes和其外部特征，查看最显著性关联
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

##行为样本，列为基因；color为等基因长度的向量
##cor The matrix of the Pearson correlations of 
##the columns of x with columns of y if y is given, 
##and the correlations of the columns of x if y is not given.
MEs0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs,datTraits,use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)

##根据相关性值着色
textMatrix <- paste(signif(moduleTraitCor,2), "\n(",
                   signif(moduleTraitPvalue,2),")",sep="")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar=c(6,8.5,3,3))
##display the correlation values within a heatmap plot
##Plots a heatmap plot with color legend, row and column annotation, 
##and optional text within th heatmap.
labeledHeatmap(Matrix=moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),##作用于legend /package(fields)
               main = paste("Module-trait relationships")
)
##图中可见多个显著性module-trait关联,这里着重于weight.

##3b Gene relationship to trait and important modules:
##Gene Significanceand Module Membership

##通过定义GS(gene significance)作为基因和特征之间的关联,查看个体基因和感兴趣(weight)特征之间的关系。
##对于每一模块，定义模块成员的MM值，作为模型eigengene和基因表达之间的相关性。
##根据以上两个值来量化所有模块的相关性
weight <- as.data.frame(datTraits$weight_g)
names(weight) <- "weight"
modNames <- substring(names(MEs),3)

geneModuleMembership <- as.data.frame(cor(datExpr,MEs,use="p"))
##Calculates Student asymptotic p-value for given correlations.
##A vector of p-values of the same length as the input correlations.
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),
                                           nSamples))
names(geneModuleMembership) <- paste("MM",modNames,sep="")
names(MMPvalue) <- paste('p.MM',modNames,sep="")

geneTraitSignificance <- as.data.frame(cor(datExpr,weight,use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
names(geneTraitSignificance) <- paste("GS.",names(weight),sep="")
names(GSPvalue) <- paste("p.GS.",names(weight),sep="")

##3c Intramodular analysis: identifying genes with high GS and MM
##根据GS和MM值，识别模块内和weight显著性相关基因和显著性相关模块

module="magenta"  ##棕色模块和weight最相关
column <- match(module,modNames)
moduleGenes <- moduleColors == module

par(mfrow=c(1,1))
##Produce a scatterplot annotated by the correlation, p-value, and regression line.
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab=paste("Module Memebership in",module,"module"),
                   ylab="Gene significance for body weigth",
                   main=paste("Module memebership vs. gene significance\n"),
                   cex.main=1.2,cex.lab=1.2,cex.axis=1.2,col=module
)
##GS MM高度相关，图示高度显著性关联一个特征的基因常是模块内最重要的单元
##位于模块中心未知

##3d Summary output of network analysis results
##合并该统计信息和基因注释信息，输出文件，核查
names(datExpr)
names(datExpr)[moduleColors == "brown"]

annot <- read.csv(file="GeneAnnotation.csv")
probes <- names(datExpr)
probes2annot <- match(probes, annot$substanceBXH)
sum(is.na(probes2annot))

##针对所有探针输出信息：probe ID/gene symbol/locus link ID(entrezid)
##module color/gene significance for weight/module membership/p-values
##构建其实数据框
geneInfo0 <- data.frame(substanceBXH=probes,
                        geneSymbol=annot$gene_symbol[probes2annot],
                        LocuslinkID=annot$LocusLinkID[probes2annot],
                        moduleColor=moduleColors,
                        geneTraitSignificance,GSPvalue
                        )
##根据其针对weight的显著性排序
modOrder <- order(-abs(cor(MEs,weight,use="p")))
##根据排序添加信息
for (mod in 1:ncol(geneModuleMembership)){
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(geneInfo0,geneModuleMembership[,modOrder],
                          MMPvalue[,modOrder])
  names(geneInfo0) <- c(oldNames,paste("MM.",modNames[modOrder],sep=""),
                        paste("p.MM.",modNames[modOrder],sep=""))
}
##先根据module color排序，再根据geneTraitSignificance排序
geneOrder <- order(geneInfo0$moduleColor,-abs(geneInfo0$GS.weight))
geneInfo0 <- geneInfo0[geneOrder,]

write.csv(geneInfo0,file="geneInfo.csv")
##使用fix函数直接在R中查看
##fix(geneInfo0)

##4 Interfacing network analysis with oter data such as functional
##annotation and gene ontology
##查询gene ontologies和显著性modul之间的生物学联系

##4a 直接导出基因列表用于gene ontology和functional enrichment分析
##例如导出brown module 基因
##4b Enrichment analysis directly within R
##可直接使用WGCNA进行GO enrichment analysis
##需安装必须的GO.db和AnnotationDBI,org.Xx.eg.db包
#GOenr <- GOenrichmentAnalysis(moduleGolors,allLLIDs,organism="mouse",nBestP=10)
#tab <- GOenr$bestPTerms[[4]]$enrichment
##略

##5 Network visualization using WGCNA functions
##Visualization of networks within R
##5a Visualizing the gene netwrok
##使用heatmap查看weighted network，每行列对应单个基因
##heatmap可描述相邻或topological的重叠，浅色代表低相邻性(重叠),
##深色代表高相邻性(重叠)。

##重新计算相应值，也可使用之前保存值
dissTOM <- 1 - TOMsimilarityFromExpr(datExpr,power=6)
##使用power转换dissTOM值，使得heatmap更明显显示较强连接
plotTOM <- dissTOM^7
##设置对角线为NA,for a nicer plot
diag(plotTOM) <- NA
#adjacency <- adjacency(datExpr,power=softPower)
#TOM <- TOMsimilarity(adjacency)
#dissTOM <- 1- TOM
#geneTree <- hclust(as.dist(dissTOM),method="average")
TOMplot(plotTOM,geneTree,moduleColors,main="Network heatmap plot, all genes")
##绘制时间和基因数目有关，同时取部分基因绘制得到的dendrogram图和所有
##基因绘制的dendrogram会看起来不同
##绘制400基因图
nSelect <- 400
set.seed(10)
select <- sample(nGenes,size=nSelect)
selectTOM <- dissTOM[select,select]
selectTree <- hclust(as.dist(selectTOM),method="average")
selectColors <- moduleColors[select]
plotDiss <- selectTOM^7
diag(plotDiss) <- NA
TOMplot(plotDiss,selectTree,selectColors,main="Network heatmap plot, selected genes")

##5b Visualizing the netwrok of eigengenes
##使用eigengenes代表profile来研究modules之间的关系
##通过eigengene correlation来量化相似性
##使用plotEigengeneNetworks绘制eigengene network，同时添加相应临床信息
MEs <- moduleEigengenes(datExpr,moduleColors)$eigengenes
weight <- as.data.frame(datTraits$weight_g)
names(weight) <- "weight"
##Reorder given (eigen-)vectors such that similar ones 
##(as measured by correlation) are next to each other.
MET <- orderMEs(cbind(MEs,weight))
par(cex=0.9)
plotEigengeneNetworks(MET,"",marDendro = c(0,4,1,2),marHeatmap = c(3,4,1,2),
                      cex.lab=0.8, xLabelsAngle=90)

##或分开绘制
par(cex=1)
plotEigengeneNetworks(MET,"Eigengene dendrogram",marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
par(cex=1.0)
plotEigengeneNetworks(MET,"Eigengenes adjacency heatmap",marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE,xLabelsAngle=90)

##6 Exporting a gene network to external visulization software
##6a Exporting to VisANT 略

##6b Exporting to Cytoscape
##Cytoscape要求输入edgefile和nodefile，用户可指定weights link以及node颜色
##例如，导出红色和棕色数据供Cytoscape绘制使用
TOM <- TOMsimilarityFromExpr(datExpr,power = 6)
annot <- read.csv(file="GeneAnnotation.csv")
modules <- c("brown",'red')
probes <- names(datExpr)
inModule <- is.finite(match(moduleColors,modules))
modProbes <- probes[inModule]
modGenes <- annot$gene_symbol[match(modProbes,annot$substanceBXH)]
modTOM <- TOM[inModule,inModule]
dimnames(modTom) <- list(modProbes,modProbes)

cyt <- exportNetworkToCytoscape(modTOM,
                                edgeFile = paste("CytoscapeInput-edges-",paste(modules,collapse="-"),
                                                 ".txt",sep=""),
                                nodeFile = paste("CytoscapeInput-nodes-",paste(modules,collapse="_"),
                                                 ".txt",sep=""),
                                weighted=TRUE,
                                threshold=0.02,
                                nodeNames=modProbes,
                                altNodeNames=modGenes,
                                nodeAttr=moduleColors[inModule])














