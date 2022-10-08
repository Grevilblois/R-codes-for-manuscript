
library(limma)
library(reshape2)
library(ggpubr)
cluFile="cluster.txt"              
immFile="CIBERSORT-Results.txt"   
pFilter=0.05        
setwd("E:\\BaiduNetdiskDownload\\ssGSEABreastcancer\\19.immCellDiff")  

Type=read.table(cluFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(Type)=c("Subtype")
Type=Type[order(Type[,"Subtype"],decreasing=T),,drop=F]
Type$Subtype=factor(Type$Subtype, levels=unique(Type$Subtype))

immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

group=sapply(strsplit(row.names(immune),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
immune=immune[group==0,]

sameSample=intersect(row.names(immune), row.names(Type))
rt=cbind(immune[sameSample,,drop=F], Type[sameSample,,drop=F])

data=melt(rt, id.vars=c("Subtype"))
colnames(data)=c("Subtype", "Immune", "Expression")

boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Subtype",
				  xlab="",
				  ylab="Fraction",
				  legend.title="Subtype",
				  width=0.8,
				  palette=c("blue","red"))+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")

pdf(file="immune.diff.pdf", width=7.5, height=6)
print(boxplot)
dev.off()

