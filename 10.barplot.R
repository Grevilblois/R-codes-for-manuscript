
library(reshape)
library(ggplot2)
fdrFilter=0.01       
setwd("E:\\BaiduNetdiskDownload\\ssGSEABreastcancer\\22.barplot")    

files=dir()
files=grep(".tsv$", files, value=T)

data=data.frame()
for(i in files){
	rt=read.table(i, header=T, sep="\t", check.names=F)
	data=rbind(data, rt)
}
data=rename(data, c("FDR q-val"= "FDR")) 
data=data[data$FDR<fdrFilter,]

labels=data[order(data$NES), "NAME"]
data$NAME=factor(data$NAME, levels=labels)

p=ggplot(data=data) + geom_bar(aes(x=NAME, y=NES, fill=FDR), stat='identity')+
    coord_flip() + 
    scale_fill_gradient(low="red", high = "blue")+ 
    xlab("") + ylab("NES") + 
    theme(axis.text.x=element_text(color="black", size=10),axis.text.y=element_text(color="black", size=10)) + 
    theme_bw()

pdf(file="barplot.pdf", width=9, height=6)
print(p)
dev.off()


p = ggplot(data,aes(NES, NAME)) + geom_point(aes(size=SIZE, color=FDR))
p1 = p + 
     scale_colour_gradient(low="red",high="blue") + 
     labs(color="FDR",size="SIZE",x="NES",y="")+
     theme(axis.text.x=element_text(color="black", size=10),axis.text.y=element_text(color="black", size=10)) + 
     theme_bw()

pdf(file="bubble.pdf", width=9, height=6)
print(p1)
dev.off()
