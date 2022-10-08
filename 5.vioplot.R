
library(reshape2)
library(ggpubr)
clusterFile="cluster.txt"
estimateFile="TMEscores.txt"
setwd("E:\\BaiduNetdiskDownload\\ssGSEABreastcancer\\16.estimateVioplot")

Type=read.table(clusterFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(Type)=c("Subtype")
Type=Type[order(Type[,"Subtype"],decreasing=T),,drop=F]
Type$Subtype=factor(Type$Subtype, levels=unique(Type$Subtype))

score=read.table(estimateFile, header=T, sep="\t", check.names=F, row.names=1)
score=score[,1:3]
score=score[row.names(Type),,drop=F]

rt=cbind(Type, score)

data=melt(rt, id.vars=c("Subtype"))
colnames(data)=c("Subtype", "scoreType", "Score")

p=ggviolin(data, x="scoreType", y="Score", fill = "Subtype",
	     xlab="",
	     ylab="TME score",
	     legend.title="Subtype",
	     add = "boxplot", add.params = list(color="white"),
	     palette = c("blue","red"), width=1)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Subtype),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

pdf(file="vioplot.pdf", width=6, height=5)
print(p1)
dev.off()
