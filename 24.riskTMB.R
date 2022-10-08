
library(limma)
library(ggpubr)
tmbFile="TMB.txt"            
riskFile="trainRisk.txt"     
setwd("E:\\BaiduNetdiskDownload\\ssGSEABreastcancer\\37.riskTMB")   


tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)
	
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
	
sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(tmb, risk)
data$TMB[data$TMB>quantile(data$TMB,0.99)]=quantile(data$TMB,0.99)
	
data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$risk))
data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	
boxplot=ggboxplot(data, x="risk", y="TMB", color="risk",
			      xlab="",
			      ylab="Tumor mutation burden",
			      legend.title="",
			      palette = c("blue", "red"),
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	
pdf(file="riskTMB.pdf", width=5.5, height=4.5)
print(boxplot)
dev.off()
