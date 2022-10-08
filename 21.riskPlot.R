
library(reshape2)
library(ggpubr)
library(ggExtra)
library(pheatmap)
setwd("E:\\BaiduNetdiskDownload\\ssGSEABreastcancer\\34.riskPlot")    

bioRiskPlot=function(inputFile=null, riskScoreFile=null, survStatFile=null, survCorFile=null){
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)  
	rt$riskScore[rt$riskScore>quantile(rt$riskScore,0.99)]=quantile(rt$riskScore,0.99)
	rt$risk=factor(rt$risk, levels=c("low", "high"))
	rt=rt[order(rt$riskScore),]     
		
	riskClass=rt[,"risk"]
	lowLength=length(riskClass[riskClass=="low"])
	highLength=length(riskClass[riskClass=="high"])
	lowMax=max(rt$riskScore[riskClass=="low"])
	line=rt[,"riskScore"]
	pdf(file=riskScoreFile, width=6, height=5)
	plot(line, type="p", pch=20,
		 xlab="Patients (increasing risk socre)", ylab="Risk score",
		 col=c(rep("blue",lowLength),rep("red",highLength)) )
	abline(h=lowMax,v=lowLength,lty=2)
	legend("topleft", c("High risk", "Low Risk"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
	dev.off()
		
	color=as.vector(rt$fustat)
	color[color==1]="red"
	color[color==0]="blue"
	pdf(file=survStatFile, width=6, height=5)
	plot(rt$futime, pch=19,
		 xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
		 col=color)
	legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
	abline(v=lowLength,lty=2)
	dev.off()
			
	x=as.numeric(rt[,"riskScore"])
	y=as.numeric(rt[,"futime"])
	df1=as.data.frame(cbind(x,y))
	p1=ggplot(df1, aes(x, y)) + 
			  xlab("Risk score") + 
			  ylab("OS (years)") +
			  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
			  stat_cor(method = 'spearman', aes(x =x, y =y))
	p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))

	pdf(file=survCorFile, width=5.2, height=5)
	print(p2)
	dev.off()
}


bioRiskPlot(inputFile="trainRisk.txt",
            riskScoreFile="train.riskScore.pdf",
            survStatFile="train.survStat.pdf",
            survCorFile="train.survCor.pdf")

bioRiskPlot(inputFile="testRisk.txt",
            riskScoreFile="test.riskScore.pdf",
            survStatFile="test.survStat.pdf",
            survCorFile="test.survCor.pdf")
