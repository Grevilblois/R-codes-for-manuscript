
library(limma)
library(sva)
tcgaExpFile="symbol.txt"        
geoExpFile="geoMatrix.txt"    
geneFile="interGenes.txt"       
setwd("E:\\BaiduNetdiskDownload\\ssGSEABreastcancer\\25.intersect")    

rt=read.table(tcgaExpFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
tcga=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
tcga=avereps(tcga)
tcga=log2(tcga+1)

rt=read.table(geoExpFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
geo=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
geo=avereps(geo)

qx=as.numeric(quantile(geo, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
    geo[geo<0]=0
    geo=log2(geo+1)}
geo=normalizeBetweenArrays(geo)

sameGene=intersect(row.names(tcga),row.names(geo))
tcgaOut=tcga[sameGene,]
geoOut=geo[sameGene,]

all=cbind(tcgaOut,geoOut)
batchType=c(rep(1,ncol(tcgaOut)),rep(2,ncol(geoOut)))
outTab=ComBat(all, batchType, par.prior=TRUE)
tcgaOut=outTab[,colnames(tcgaOut)]
tcgaOut[tcgaOut<0]=0
geoOut=outTab[,colnames(geoOut)]
geoOut[geoOut<0]=0

gene=read.table(geneFile, header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(tcgaOut))
tcgaShareExp=tcgaOut[sameGene,]
geoShareExp=geoOut[sameGene,]

tcgaShareExp=rbind(ID=colnames(tcgaShareExp),tcgaShareExp)
write.table(tcgaShareExp,file="tcga.share.txt",sep="\t",quote=F,col.names=F)
geoShareExp=rbind(ID=colnames(geoShareExp),geoShareExp)
write.table(geoShareExp,file="geo.share.txt",sep="\t",quote=F,col.names=F)
