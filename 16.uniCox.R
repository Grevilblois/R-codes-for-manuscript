
library(survival)                 
coxPfilter=0.05                   
inputFile="tcga.expTime.txt"    
setwd("E:\\BaiduNetdiskDownload\\ssGSEABreastcancer\\28.uniCox")     
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)   
rt$futime=rt$futime/365


outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
	
	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	
	if(coxP<coxPfilter){
	    sigGenes=c(sigGenes,i)
		outTab=rbind(outTab,
			         cbind(id=i,
			         HR=coxSummary$conf.int[,"exp(coef)"],
			         HR.95L=coxSummary$conf.int[,"lower .95"],
			         HR.95H=coxSummary$conf.int[,"upper .95"],
			         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
			        )
	}
}


write.table(outTab,file="tcga.uniCox.txt",sep="\t",row.names=F,quote=F)


uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="tcga.uniSigExp.txt",sep="\t",row.names=F,quote=F)


bioForest=function(coxFile=null,forestFile=null,forestCol=null){

	rt <- read.table(coxFile,header=T,sep="\t",check.names=F,row.names=1)
	gene <- rownames(rt)
	hr <- sprintf("%.3f",rt$"HR")
	hrLow  <- sprintf("%.3f",rt$"HR.95L")
	hrHigh <- sprintf("%.3f",rt$"HR.95H")
	Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
	pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

	pdf(file=forestFile, width = 6.3,height = 5)
	n <- nrow(rt)
	nRow <- n+1
	ylim <- c(1,nRow)
	layout(matrix(c(1,2),nc=2),width=c(3,2.5))

	xlim = c(0,3)
	par(mar=c(4,2.5,2,1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
	text.cex=0.8
	text(0,n:1,gene,adj=0,cex=text.cex)
	text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
	text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

	par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
	xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
	arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
	abline(v=1,col="black",lty=2,lwd=2)
	boxcolor = ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])
	points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
	axis(1)
	dev.off()
}

bioForest(coxFile="tcga.uniCox.txt",forestFile="forest.pdf",forestCol=c("red","green"))
