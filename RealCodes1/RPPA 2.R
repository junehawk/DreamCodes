library(multicore)
rppa = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_RPPA.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1);
ds = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_Drug_Response_Training.txt", header=TRUE, row.names=1);
rppa_cellines = colnames(rppa)
drug_cellines = colnames(ds)
rp = subset(rppa, select=intersect(drug_cellines, rppa_cellines))
dss = subset(ds, select=intersect(drug_cellines, rppa_cellines))
binded = rbind(rp, dss)
mydata = as.data.frame(t(binded))

FDR = dim(rp)[1]
num_rppa = dim(rp)[1]
num_drug = dim(dss)[1]
#num_drug = 25
drug_start = num_rppa

lmpvalue <- function(x, i, data = mydata) {
	#cat(i, "\t", x[[1]], "\n")
	fit = lm(mydata[,i+drug_start]~mydata[,x[[1]]])
	if(dim(summary(fit)$coefficient)[1] >=2) {
		if(!is.na(summary(fit)$coefficient[2,4]))
			summary(fit)$coefficient[2,4]
		else
			1
	}
	else
		1
}
rm(ds)
rm(dss)
rm(rp)
rm(binded)
x = as.list(c(1:num_rppa))
for (i in 1:num_drug) {
	result_pvalue = unlist(mclapply(x,lmpvalue,i=i))
	corrected = p.adjust(result_pvalue, method="fdr", n=num_rppa)
	
	result = cbind(rppa, corrected)
	
	result = subset(result, corrected <=0.05)

	cat("Drug ", i, " : ", dim(result)[1])

	if(dim(result)[1]>0) {
		write_addr = paste("Dream7/input1/Drug",as.character(i),"_rppa.txt",sep="")
		if(dim(result)[1]>1)
			result = result[order(subset(result, select=c(corrected))),]
		write.table(result[,-dim(result)[2]],file=write_addr,sep="\t",row.names = T)
	}
	rm(result_pvalue)
	rm(corrected)
	rm(result)
}
rm(mydata)
rm(rppa)
