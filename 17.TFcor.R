
library(limma)
library(ggalluvial)
library(ggplot2)
library(dplyr)
corFilter=0.4               
fdrFilter=0.001              
expFile="diffGeneExp.txt"    
tfFile="TF.txt"              
coxFile="tcga.uniCox.txt"     
setwd("E:\\BaiduNetdiskDownload\\ssGSEABreastcancer\\29.TFcor")   

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.1,]

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
data=data[,group==0]

TF=read.table(tfFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data), as.vector(TF[,1]))
TFexp=data[sameGene,,drop=F]

PIG=read.table(coxFile, header=T, sep="\t", check.names=F)
sameGene=intersect(row.names(data), as.vector(PIG[,1]))
PIGexp=data[sameGene,,drop=F]

outTab=data.frame()
for(i in row.names(PIGexp)){
	if(sd(PIGexp[i,])>0.01){
		for(j in row.names(TFexp)){
			if(i==j){next}
			x=as.numeric(PIGexp[i,])
			y=as.numeric(TFexp[j,])
			corT=cor.test(x, y)
			cor=corT$estimate
			pvalue=corT$p.value
			outTab=rbind(outTab,cbind(PIGs=i, TFs=j, cor, pvalue))
		}
	}
}

pvalue=outTab[,"pvalue"]
fdr=p.adjust(as.numeric(as.vector(pvalue)), method="fdr")
outTab=cbind(outTab, FDR=fdr)
outTab=outTab[(as.numeric(as.vector(outTab$FDR))<fdrFilter & abs(as.numeric(as.vector(outTab$cor)))>corFilter ),]
write.table(file="cor.result.txt",outTab,sep="\t",quote=F,row.names=F)

node=unique(c(outTab[,1],outTab[,2]))
write.table(file="cor.node.txt", node, sep="\t", quote=F, row.names=F, col.names=F)

rt=outTab[,c(1,2)]
corLodes=to_lodes_form(rt, axes=1:ncol(rt), id="Cohort")

pdf(file="ggalluvial.pdf", width=4, height=6)
mycol <- rep(c("#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  	 scale_x_discrete(expand = c(0, 0)) +  
  	 
  	 geom_flow(width = 2/10,aes.flow = "forward") + 
	 geom_stratum(alpha = .9,width = 2/10) +
	 scale_fill_manual(values = mycol) +

	 geom_text(stat = "stratum", size = 2,color="black") +
	 xlab("") + ylab("") + theme_bw() + 
	 theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #È¥µô×ø±êÖá
	 theme(panel.grid =element_blank()) + 
	 theme(panel.border = element_blank()) + 
	 ggtitle("") + guides(fill = FALSE)                            
dev.off()
