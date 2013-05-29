library(multicore)
mrna = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_GeneExpression.txt", header=TRUE, row.names=1, stringsAsFactors=FALSE);
ds = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_Drug_Response_Training.txt", header=TRUE, row.names=1);
rna_cellines = colnames(mrna)
drug_cellines = colnames(ds)
rna = subset(mrna, select=intersect(drug_cellines, rna_cellines))
dss = subset(ds, select=intersect(drug_cellines, rna_cellines))
binded = rbind(rna, dss)
mydata = as.data.frame(t(binded))

FDR = dim(rna)[1]
num_rna = dim(rna)[1]
num_drug = dim(dss)[1]
#num_drug = 1
drug_start = num_rna

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
rm(binded)
rm(rna)
x = as.list(c(1:num_rna))
for (i in 1:num_drug) {
	cat("Drug ", i, "\n")
	

	result_pvalue = unlist(mclapply(x,lmpvalue,i=i))
	corrected = p.adjust(result_pvalue, method="fdr", n=num_rna)
	
	result = cbind(mrna, corrected)
	
	result = subset(result, corrected <=0.05)

	cat("Drug ", i, " : ", dim(result)[1], "\n")

	if(dim(result)[1]>0) {
		write_addr = paste("Dream7/input1/Drug",as.character(i),"_expression.txt",sep="")
		write.table(result[,-dim(result)[2]],file=write_addr,sep="\t",row.names = T)
	}
	rm(result_pvalue)
	rm(corrected)
	rm(result)
}
rm(mydata)
rm(mrna)