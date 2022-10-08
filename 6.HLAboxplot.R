
#ÒýÓÃ°ü
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
expFile="symbol.txt"
cluFile="cluster.txt"
geneFile="gene.txt"
setwd("E:\\BaiduNetdiskDownload\\ssGSEABreastcancer\\17.HLAboxplot")

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data),as.vector(gene[,1]))
data=t(data[sameGene,])
data=log2(data+1)

group=sapply(strsplit(row.names(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[group==0,]

Type=read.table(cluFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(Type)=c("Subtype")
Type=Type[order(Type[,"Subtype"],decreasing=T),,drop=F]
Type$Subtype=factor(Type$Subtype, levels=unique(Type$Subtype))

sameSample=intersect(row.names(data), row.names(Type))
rt1=cbind(data[sameSample,,drop=F],Type[sameSample,,drop=F])

data=melt(rt1, id.vars=c("Subtype"))
colnames(data)=c("Subtype", "Gene", "Expression")

boxplot=ggboxplot(data, x="Gene", y="Expression", fill="Subtype",
				  xlab="",
				  ylab="Gene expression",
				  legend.title="Subtype",
				  width=0.8,
				  palette=c("blue", "red") )+
				  rotate_x_text(45)+
	stat_compare_means(aes(group=Subtype),
		method="wilcox.test",
		symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")

pdf(file="boxplot.pdf", width=8, height=5.6)
print(boxplot)
dev.off()
