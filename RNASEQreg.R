library(multicore)
rnaseq = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_RNAseq_quantification.txt", header=TRUE, row.names=2, stringsAsFactors=FALSE);
ds = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_Drug_Response_Training.txt", header=TRUE, row.names=1);
rna_cellines = colnames(rnaseq)
drug_cellines = colnames(ds)
rna = subset(rnaseq, select=intersect(drug_cellines, rna_cellines))
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

result_gene = rnaseq[,1]
rm(rnaseq)
rm(ds)
rm(rna)
rm(dss)
rm(binded)
x = as.list(c(1:num_rna))
for (i in 1:num_drug) {
	cat("Drug ", i, "\n")
	

	result_pvalue = unlist(mclapply(x,lmpvalue,i=i))
	corrected = p.adjust(result_pvalue, method="fdr", n=num_rna)
	
	result = cbind(result_gene, result_pvalue, corrected)
	
	toWrite = subset(result, result_pvalue<=0.05)

	candidate = subset(result, corrected<=0.05)
	cat("Drug ", i, " : ", dim(candidate)[1])
	rm(candidate)

	if(length(toWrite)>0) {
		write_addr = paste("Dream7/test/Drug",as.character(i),"_rna_reg_2.txt",sep="")
		write.table(toWrite[order(subset(toWrite, select=c(corrected))),],file=write_addr,sep="\t",row.names = F)
	}
	rm(result_pvalue)
	rm(corrected)
	rm(result)
	rm(toWrite)
}
rm(mydata)