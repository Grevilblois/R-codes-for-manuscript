
library(pheatmap)  
ssgseaFile="ssgseaOut.txt" 
clusterFile="cluster.txt"  
estimateFile="TMEscores.txt"
setwd("E:\\BaiduNetdiskDownload\\ssGSEABreastcancer\\15.estimateHM") 

Type=read.table(clusterFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(Type)=c("Subtype")
Type=Type[order(Type[,"Subtype"],decreasing=T),,drop=F]
Type$Subtype=factor(Type$Subtype, levels=unique(Type$Subtype))

rt=read.table(ssgseaFile, header=T, sep="\t", check.names=F, row.names=1)
rt=rt[,row.names(Type)]

score=read.table(estimateFile, header=T, sep="\t", check.names=F, row.names=1)
score=score[row.names(Type),,drop=F]

cluster=cbind(Type, score)

ann_colors=list()
clusterCol=c("blue", "red")
names(clusterCol)=levels(factor(Type$Subtype))
ann_colors[["Subtype"]]=clusterCol

pdf("estimateHM.pdf", width=9, height=5)
pheatmap(rt, annotation=cluster,
		annotation_colors = ann_colors,
        color = colorRampPalette(c(rep("blue",3), "white", rep("red",3)))(50),
        cluster_cols =F,
        scale="row",
        show_colnames=F,
        fontsize=8,
        fontsize_row=8,
        fontsize_col=3)
dev.off()
