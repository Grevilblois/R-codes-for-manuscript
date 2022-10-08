
library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)
immFile="CIBERSORT-Results.txt"     
riskFile="trainRisk.txt"            
setwd("E:\\BaiduNetdiskDownload\\ssGSEABreastcancer\\35.immuneCor")    

immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<0.05,]
data=as.matrix(immune[,1:(ncol(immune)-3)])

group=sapply(strsplit(row.names(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=avereps(data)

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,3:(ncol(risk)-2),drop=F]

outTab=data.frame()
for(immune in colnames(data)){
	for(gene in colnames(risk)){
		x=as.numeric(data[,immune])
		y=as.numeric(risk[,gene])
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pvalue=corT$p.value
		text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
		outTab=rbind(outTab,cbind(Gene=gene, Immune=immune, cor, text, pvalue))
	}
}

outTab$cor=as.numeric(outTab$cor)
pdf(file="cor.pdf", width=7, height=6)
ggplot(outTab, aes(Gene, Immune)) + 
	geom_tile(aes(fill = cor), colour = "grey", size = 1)+
	scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
	geom_text(aes(label=text),col ="black",size = 3) +
	theme_minimal() +   
	theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),   
	      axis.text.y = element_text(size = 10, face = "bold")) +       #yÖá×ÖÌå
	labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   
	scale_x_discrete(position = "bottom")      
dev.off()
