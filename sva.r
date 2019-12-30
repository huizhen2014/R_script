##
##
##sva
library(sva)
library(bladderbatch)
data(bladderdata)
library(pamr)
library(limma)

##setting up the data from ExpressionSet
pheno <- pData(bladderEset)
edata <- exprs(bladderEset)

mod <- model.matrix(~as.factor(cancer),data=pheno)
mod0 <- model.matrix(~1, data=pheno)

##applying the sva function to estimate batch and other artifacts
n.sv <- num.sv(edata, mod, method="leek")

svobj <- sva(edata, mod, mod0, n.sv=n.sv)

##Adjusting for surrogate variables using the f.pvalue function
pValues <- f.pvalue(edata, mod, mod0)
qvalues <- p.adjust(pValues, method="BH")

modSv <- cbind(mod,svobj$sv)
mod0Sv <- cbind(mod0,svobj$sv)
pValuesSv <- f.pvalue(edata, modSv, mod0Sv)
qValuesSv <- p.adjust(pValuesSv,method="BH")

##Adjusting for surrogate variables using the limma package
fit <- lmFit(edata, modSv)

contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,svobj$n.sv)),
                         "C2"=c(0,-1,1,rep(0,svobj$n.sv)))
fitContrasts <- contrasts.fit(fit, contrast.matrix)

eb <- eBayes(fitContrasts)
topTableF(eb, adjust="BH")

##Applying the ComBat function to adjust for known batches
batch <- pheno$batch

modcombat <- model.matrix(~1, data=pheno)
combat_edata <- ComBat(dat=edata, batch=batch,
                       mod=modcombat,par.prior=T,
                       prior.plots = T)

pValuesComBat <- f.pvalue(combat_edata,mod,mod0)
qValuesComBat <- p.adjust(pValuesComBat,method="BH")

##Removing known batch effects with a linear model
modBatch <- model.matrix(~ as.factor(cancer) + as.factor(batch),
                         data=pheno)
mod0Batch <- model.matrix(~as.factor(batch), data=pheno)
pValuesBatch <- f.pvalue(edata, modBatch, mod0Batch)
qValuesBatch <- p.adjust(pValuesBatch, method="BH")


##Applying the fsva function to remove batch effects for perdiction
set.seed(12345)
trainIndicator <- sample(1:57,size=30,replace=FALSE)
testIndicator <- (1:57)[-trainIndicator]
trainData <- edata[,trainIndicator]
testData <- edata[,testIndicator]
trainPheno <- pheno[trainIndicator,]
testPhene <- pheno[testIndicator,]

library(pamr)
mydata <- list(x=trainData,y=trainPheno$cancer)
mytrain <- pamr.train(mydata)
table(pamr.predict(mytrain, testData,threshold = 2), testPhene$cancer)

trainMod <- model.matrix(~cancer, data=trainPheno)
trainMod0 <- model.matrix(~1, data=trainPheno)
trainSv <- sva(trainData, trainMod, trainMod0)

fsvaobj <- fsva(trainData, trainMod, trainSv, testData)
mydataSv <- list(x=fsvaobj$db, y=trainPheno$cancer)
mytrainSv <- pamr.train(mydataSv)
table(pamr.predict(mytrainSv,fsvaobj$new,threshold = 1),
      testPhene$cancer)

##sva for sequencing (svaseq)
library(zebrafishRNASeq)
data("zfGenes")
filter <- apply(zfGenes,1,function(x)length(x[x>5])>=2)
filtered <- zfGenes[filter,]
genes <- rownames(filtered)[grep("^ENS",rownames(filtered))]
controls <- grepl("^ERCC",rownames(filtered))
group <- as.factor(rep(c("Ctl","Trt"),each=3))
dat0 <- as.matrix(filtered)

##set null and alternative models (ignore batch), as the sva needed
mod1 <- model.matrix(~group)
mod0 <- cbind(mod1[,1])
svseq <- svaseq(dat0, mod1, mod0,n.sv=1)$sv
plot(svseq, pch=19, col="blue")

#suprevised svaseq function
sup_svseq <- svaseq(dat0, mod1, mod0, controls = controls,n.sv=1)$sv
plot(sup_svseq,svseq,pch=19,col="blue")














