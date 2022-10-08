
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
gene="BRCA1"                 
expFile="symbol.txt"        
riskFile="trainRisk.txt"    
setwd("E:\\BaiduNetdiskDownload\\ssGSEABreastcancer\\36.riskGene")    

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
	
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
	
data=rbind(data, gene=data[gene,])
exp=t(data[c("gene",gene),])
row.names(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\",  row.names(exp))
exp=avereps(exp)
	
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
	
sameSample=intersect(row.names(exp), row.names(risk))
exp=exp[sameSample,]
exp=log2(exp+1)
risk=risk[sameSample,]
data=cbind(as.data.frame(exp), as.data.frame(risk))

data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$risk))
data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	
boxplot=ggboxplot(data, x="risk", y="gene", color="risk",
			      xlab="",
			      ylab=paste0(gene, " expression"),
			      legend.title="",
			      palette = c("blue", "red"),
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)	

pdf(file=paste0(gene, ".boxplot.pdf"), width=5, height=4.5)
print(boxplot)
dev.off()
	
xlab="riskScore"
ylab=gene
x=as.numeric(data[,xlab])
y=as.numeric(data[,ylab])
df1=as.data.frame(cbind(x,y))
p1=ggplot(df1, aes(x, y)) + 
		  xlab("Risk score") + ylab(paste0(gene, " expression"))+
		  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
		  stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))

pdf(file=paste0(gene, ".cor.pdf"), width=5.2, height=5)
print(p2)
dev.off()
