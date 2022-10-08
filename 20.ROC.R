
library(survival)
library(survminer)
library(timeROC)
library(rms)
setwd("E:\\BaiduNetdiskDownload\\ssGSEABreastcancer\\33.ROC")     

bioROC=function(inputFile=null, rocFile=null, calFile=null){

	rt=read.table(inputFile, header=T, sep="\t", check.names=F)

	ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
	               marker=rt$riskScore, cause=1,
	               weighting='aalen',
	               times=c(1,3,5), ROC=TRUE)
	pdf(file=rocFile,width=5,height=5)
	plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
	plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
	plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
	legend('bottomright',
	        c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
	          paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
	          paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
	        col=c("green","blue","red"),lwd=2,bty = 'n')
	dev.off()
	
	pdf(file=calFile, width=5, height=5)

	f <- cph(Surv(futime, fustat) ~ riskScore, x=T, y=T, surv=T, data=rt, time.inc=1)
	cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
	plot(cal, xlim=c(0,1), ylim=c(0,1),
		 xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="green", sub=F)

	f <- cph(Surv(futime, fustat) ~ riskScore, x=T, y=T, surv=T, data=rt, time.inc=3)
	cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
	plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="blue", sub=F, add=T)

	f <- cph(Surv(futime, fustat) ~ riskScore, x=T, y=T, surv=T, data=rt, time.inc=5)
	cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
	plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="red", sub=F, add=T)
	legend('bottomright', c('1-year', '3-year', '5-year'),
	        col=c("green","blue","red"), lwd=1.5, bty = 'n')
	dev.off()
}

bioROC(inputFile="trainRisk.txt", rocFile="train.ROC.pdf", calFile="train.cal.pdf")
bioROC(inputFile="testRisk.txt", rocFile="test.ROC.pdf", calFile="test.cal.pdf")
