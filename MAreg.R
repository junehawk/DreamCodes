# TODO: Add comment
# 
# Author: wchwang
###############################################################################


library(multicore)
data = read.table("/home/wchwang/Dream/data/DREAM7_1/DREAM7_DrugSensitivity1_RPPA_v1.txt", header=TRUE, stringsAsFactors=FALSE);
data_info = subset(data, select=c(HCC1419,AU565,BT549,BT483,MCF7,CAMA1,HCC1395,HCC1806,MCF10F,MCF12A,UACC812,MDAMB157,HCC38,HCC70,MDAMB361,BT20,HCC1569,	HCC202,	ZR7530,	HCC1428,T47D,HCC1937,HCC1954,BT474,HCC1143,MDAMB231,MDAMB453,MDAMB415))
ds = read.table("/home/wchwang/Dream/data/DREAM7_1/DREAM7_1_Drug_Response_Training.txt", header=TRUE, row.names=1);
dss = subset(ds, select=c(HCC1419,AU565,BT549,BT483,MCF7,CAMA1,HCC1395,HCC1806,MCF10F,MCF12A,UACC812,MDAMB157,HCC38,HCC70,MDAMB361,BT20,HCC1569,	HCC202,	ZR7530,	HCC1428,T47D,HCC1937,HCC1954,BT474,HCC1143,MDAMB231,MDAMB453,MDAMB415))
binded = rbind(data_info, dss)
mydata = as.data.frame(t(binded))

FDR = dim(data_info)[1]
num_data= dim(data_info)[1]
num_drug = dim(dss)[1]
#num_drug = 1
drug_start = num_data

lmpvalue <- function(x, i, data = mydata) {
#	cat(i, "\t", x[[1]], "\n")
	fit = lm(mydata[,i+drug_start]~mydata[,x[[1]]+1])
	if(dim(summary(fit)$coefficient)[1] >=2) {
		if(!is.na(summary(fit)$coefficient[2,4]))
			summary(fit)$coefficient[2,4]
		else
			1
	}
	else
		1
}


for (i in 1:num_drug) {
	cat("Drug ", i, "\n")
	x = as.list(c(1:num_data))
	
	result_pvalue = unlist(mclapply(x,lmpvalue,i=i))
	result_gene = data[,1]
	corrected = p.adjust(result_pvalue, method="fdr", n=num_data)
	
	result = cbind(result_gene, result_pvalue, corrected)
	
	toWrite = subset(result, result_pvalue<=0.05)
	if(length(toWrite)>0) {
		write_addr = paste("/home/wchwang/Dream/data/DREAM7_1/RPPA/drug",as.character(i),"_RPPA_reg.txt",sep="")
		write.table(toWrite[order(subset(toWrite, select=c(corrected))),],file=write_addr,sep="\t",row.names = F)
	}
}
