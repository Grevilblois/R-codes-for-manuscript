
library(sparcl)
library(Rtsne)
library(ggplot2)
inputFile="ssgseaOut.txt"
setwd("E:\\BaiduNetdiskDownload\\ssGSEABreastcancer\\13.cluster") 

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]

hc=hclust(dist(t(data)))

k=2
y=cutree(hc, k)

median1=median(as.matrix(data[,y==1]))
median2=median(as.matrix(data[,y==2]))
if(median1>median2){
	y1=ifelse(y==1, "Immunity_H", "Immunity_L")
}else{
	y1=ifelse(y==1, "Immunity_L", "Immunity_H")
}
write.table(y1, file="cluster.txt", sep="\t", quote=F, col.names=F)

y2=ifelse(y1=="Immunity_H", "red", "blue")
pdf(file="hclust.pdf", width=30, height=15)
ColorDendrogram(hc, y=y2,
                labels=names(y), branchlength=0.3, 
                xlab=" ", sub=" ", main = " ")
dev.off()

tsneOut=Rtsne(t(data), dims=2, perplexity=10, verbose=F, max_iter=500,check_duplicates=F)
tsne=data.frame(tSNE1 = tsneOut$Y[,1], tSNE2 = tsneOut$Y[,2], Subtype=y1)	

pdf(file="tSNE.pdf", width=6.5, height=5)
p=ggplot(data = tsne, aes(tSNE1, tSNE2)) + geom_point(aes(color = Subtype)) +
		scale_colour_manual(name="Subtype",  values =c("red", "blue"))+
	    theme_bw()+
	    theme(plot.margin=unit(rep(1.5,4),'lines'))+
	    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()

