
library(venn)
library(pheatmap)
diffFile="diffGeneExp.txt"      
immGeneFile="gene.txt"          
setwd("E:\\BaiduNetdiskDownload\\ssGSEABreastcancer\\24.immGeneDiff")    
geneList=list()

data=read.table(diffFile, header=T, sep="\t", check.names=F, row.names=1)
geneNames=gsub("^ | $", "", row.names(data))      
uniqGene=unique(geneNames)                        
geneList[["DEGs"]]=uniqGene

rt=read.table(immGeneFile, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)      
uniqGene=unique(geneNames)              
geneList[["Immune"]]=uniqGene

mycol=c("#029149","#E0367A","#5D90BA","#431A3D","#FFD121","#D8D155","#223D6C","#D20A13","#088247","#11AA4D","#7A142C","#5D90BA","#64495D","#7CC767")
pdf(file="venn.pdf", width=5, height=5)
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F,ilabels=F)
dev.off()

interGenes=Reduce(intersect,geneList)
write.table(file="interGenes.txt",interGenes,sep="\t",quote=F,col.names=F,row.names=F)

group=sapply(strsplit(colnames(data),"\\_"), "[", 2)
Subtype=paste0("Immunity_", group)
names(Subtype)=colnames(data)
Type=as.data.frame(Subtype)
Type$Subtype=factor(Type$Subtype, levels=unique(Type$Subtype))

ann_colors=list()
clusterCol=c("blue", "red")
names(clusterCol)=levels(factor(Type$Subtype))
ann_colors[["Subtype"]]=clusterCol

data=log2(data+0.01)
pdf(file="heatmap.pdf", width=8, height=6.5)
pheatmap(data, 
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(50),
         cluster_cols =F,
         show_colnames = F,
         show_rownames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=8,
         fontsize_col=8)
dev.off()

data=data[interGenes,]
pdf(file="immune.heatmap.pdf", width=8, height=6.5)
pheatmap(data, 
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(50),
         cluster_cols =F,
         show_colnames = F,
         show_rownames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=8,
         fontsize_col=8)
dev.off()
