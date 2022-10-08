
library(limma)
library(ggplot2)
expFile="symbol.txt"        
cluFile="cluster.txt"       
logFCfilter=0.585         
fdrFilter=0.05           
setwd("E:\\BaiduNetdiskDownload\\ssGSEABreastcancer\\23.clusterDiff")   

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]

Type=read.table(cluFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(Type)=c("Subtype")
sameSample=intersect(colnames(data), row.names(Type))
data=data[,sameSample,drop=F]
Type=Type[sameSample,,drop=F]

immLow=Type[Type$Subtype=="Immunity_L",,drop=F]
immHigh=Type[Type$Subtype=="Immunity_H",,drop=F]
dataLow=data[,row.names(immLow)]
colnames(dataLow)=paste0(colnames(dataLow),"_L")
dataHigh=data[,row.names(immHigh)]
colnames(dataHigh)=paste0(colnames(dataHigh),"_H")
data=cbind(dataLow,dataHigh)
data=data[rowMeans(data)>1,]
conNum=ncol(dataLow)
treatNum=ncol(dataHigh)
Type=c(rep(1,conNum), rep(2,treatNum))

outTab=data.frame()
for(i in row.names(data)){
	rt=data.frame(expression=data[i,], Type=Type)
	wilcoxTest=wilcox.test(expression ~ Type, data=rt)
	pvalue=wilcoxTest$p.value
	conGeneMeans=mean(data[i,1:conNum])
	treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
	logFC=log2(treatGeneMeans)-log2(conGeneMeans)
	conMed=median(data[i,1:conNum])
	treatMed=median(data[i,(conNum+1):ncol(data)])
	diffMed=treatMed-conMed
	if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
		outTab=rbind(outTab,cbind(gene=i,lowMean=conGeneMeans,highMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
	}
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
outTab=cbind(outTab, fdr=fdr)

outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff, file="diff.txt", sep="\t", row.names=F, quote=F)

diffExp=rbind(ID=colnames(data[as.vector(outDiff[,1]),]), data[as.vector(outDiff[,1]),])
write.table(diffExp, file="diffGeneExp.txt", sep="\t", col.names=F, quote=F)

outTab$fdr=as.numeric(outTab$fdr)
outTab$logFC=as.numeric(outTab$logFC)
Significant=ifelse((outTab$fdr<fdrFilter & abs(outTab$logFC)>logFCfilter), ifelse(outTab$logFC>logFCfilter,"Up","Down"), "Not")

p = ggplot(outTab, aes(logFC, -log10(fdr)))+
    geom_point(aes(col=Significant))+
    scale_color_manual(values=c("green", "black", "red"))+
    labs(title = " ")+
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()

pdf("vol.pdf", width=6.2, height=5.5)
print(p)
dev.off()
