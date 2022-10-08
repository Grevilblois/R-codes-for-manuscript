
library("glmnet")
library("survival")
coxSigFile="tcga.uniSigExp.txt"     
geoFile="geo.expTime.txt"        
setwd("E:\\BaiduNetdiskDownload\\ssGSEABreastcancer\\31.model")  

rt=read.table(coxSigFile, header=T, sep="\t", check.names=F, row.names=1)
geo=read.table(geoFile, header=T, sep="\t", check.names=F, row.names=1)
sameGene=intersect(colnames(rt)[3:ncol(rt)], colnames(geo)[3:ncol(geo)])
rt=rt[,c("futime","fustat",sameGene)]
rt$futime[rt$futime<=0]=0.003

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime, rt$fustat))
fit=glmnet(x, y, family="cox", maxit=1000)

pdf("lasso.lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit=cv.glmnet(x, y, family="cox", maxit=1000)
pdf("lasso.cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef=coef(fit, s=cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
write.table(geneCoef, file="lasso.geneCoef.txt", sep="\t", quote=F, row.names=F)

trainFinalGeneExp=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="trainRisk.txt",sep="\t",quote=F,row.names=F)

rt=read.table(geoFile, header=T, sep="\t", check.names=F, row.names=1)
rt$futime=rt$futime/365
testFinalGeneExp=rt[,lassoGene]
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="testRisk.txt",sep="\t",quote=F,row.names=F)
