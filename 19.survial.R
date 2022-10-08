
library(survival)
library(survminer)
setwd("E:\\BaiduNetdiskDownload\\ssGSEABreastcancer\\32.survival")  

bioSurvival=function(inputFile=null,outFile=null){
	
	rt=read.table(inputFile,header=T,sep="\t",check.names=F)

	diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
		

	surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=T,
		           pval=pValue,
		           pval.size=6,
		           legend.title="Risk",
		           legend.labs=c("High risk", "Low risk"),
		           xlab="Time(years)",
		           break.time.by = 1,
		           palette=c("red", "blue"),
		           risk.table=TRUE,
		           risk.table.title="",
		           risk.table.height=.25)
	pdf(file=outFile,onefile = FALSE,width = 6.5,height =5.5)
	print(surPlot)
	dev.off()
}
bioSurvival(inputFile="trainRisk.txt", outFile="trainSurv.pdf")
bioSurvival(inputFile="testRisk.txt", outFile="testSurv.pdf")
