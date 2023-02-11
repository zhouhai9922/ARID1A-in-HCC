######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

#install.packages("survival")
#install.packages("survminer")
#install.packages("timeROC")


#???ð?
library(survival)
library(survminer)
library(timeROC)

predictTime=5       #Ԥ??ʱ??
riskFile="nomoRisk.txt"      #????ͼ?ķ????ļ?
cliFile="clinical.txt"       #?ٴ??????ļ?
setwd("C:\\Users\\ROC-NOM")     #?޸Ĺ???Ŀ¼

#??ȡ?????????ļ?
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#??ȡ?ٴ??????ļ?
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#?ϲ?????
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli1=cli[samSample,,drop=F]
data=cbind(risk1, cli1)
rt=cbind(risk1[,c("futime","fustat","riskScore","Nomogram")], cli1)

######????ROC????######
aucText=c()
bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)
pdf(file="cliROC.pdf", width=6, height=6)
#???Ʒ??յ÷ֵ?ROC????
i=3
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)

#???ٴ????ݽ???ѭ?????????ٴ????ݵ?ROC????
for(i in 4:ncol(rt)){
	ROC_rt=timeROC(T=rt$futime,
				   delta=rt$fustat,
				   marker=rt[,i], cause=1,
				   weighting='aalen',
				   times=c(predictTime),ROC=TRUE)
	plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2, add=TRUE)
	aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}
#????ͼ?????õ?ROC?????µ?????
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()


######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

