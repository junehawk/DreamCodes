library(multicore)
meth = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_Methylation.txt", header=TRUE, row.names=1, stringsAsFactors=FALSE);
ds = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_Drug_Response_Training.txt", header=TRUE, row.names=1);
meth_cellines = colnames(meth)
drug_cellines = colnames(ds)
mt = subset(meth, select=intersect(drug_cellines, meth_cellines))
dss = subset(ds, select=intersect(drug_cellines, meth_cellines))
binded = rbind(mt, dss)
mydata = as.data.frame(t(binded))

FDR = dim(mt)[1]
num_mt = dim(mt)[1]
num_drug = dim(dss)[1]
#num_drug = 1
drug_start = num_mt

lmpvalue <- function(x, i, data = mydata) {
	cat("x[[1]]\t", x[[1]], "x\t", x, "\n")
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

result_gene = meth[,1]

rm(meth)
rm(ds)
rm(mt)
rm(dss)
rm(binded)

x = as.list(c(1:num_mt))
for (i in 1:num_drug) {
	#cat("Drug ", i, "\n")
	

	result_pvalue = unlist(mclapply(x,lmpvalue,i=i))	
	corrected = p.adjust(result_pvalue, method="fdr", n=num_mt)
	
	result = cbind(result_gene, result_pvalue, corrected)
	
	toWrite = subset(result, result_pvalue<=0.05)

	candidate = subset(result, corrected<=0.05)
	cat("Drug ", i, " : ", dim(candidate)[1], "\n")
	rm(candidate)

	if(length(toWrite)>0) {
		write_addr = paste("Dream7/test/Drug",as.character(i),"_meth_reg.txt",sep="")
		write.table(toWrite[order(subset(toWrite, select=c(corrected))),],file=write_addr,sep="\t",row.names = F)
	}
	rm(result_pvalue)
	rm(corrected)
	rm(result)
	rm(toWrite)
}
rm(mydata)