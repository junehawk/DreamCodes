library(elasticnet)
library(multicore)

rnaseq = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_RNAseq_quantification.txt", header=TRUE, row.names=2, stringsAsFactors=FALSE);
ds = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_Drug_Response_Training.txt", header=TRUE, row.names=1);
rna_cellines = colnames(rnaseq)
drug_cellines = colnames(ds)
rna = subset(rnaseq, select=intersect(drug_cellines, rna_cellines))
dss = subset(ds, select=intersect(drug_cellines, rna_cellines))

gene_names = subset(rnaseq, select=c(HGNC_ID))
num_drug = dim(dss)[1]
#num_drug = 25
gene_start = 27
rm(rnaseq)
rm(ds)

lassofit <- function(i, ddata = dss, rdata = rna) {
	y = t(ddata[i[[1]],])
	y2 = subset(y, !is.na(y))
	x = t(subset(rdata, select=intersect(rownames(y2), colnames(rna))))
	object = enet(x, y2, lambda=0, normalize=FALSE)
	#if(class(result) == “try-error”)
	#	next;
	rm(y)
	rm(y2)
	rm(x)
	object
}


index = as.list(c(gene_start:num_drug))
result = mclapply(index, lassofit)

for(i in 1:length(result)) {

	Cp = as.vector(result[[i]]$Cp)
	Cp = Cp[1:(length(Cp)-1)]
	num_term = length(result[[i]]$actions)
	genes = c()
	for(j in 1:(num_term-1)) {
		genes = append(gene_names[names(result[[i]]$actions[[j]]),], genes, after=0)
	}
	toWrite = cbind(genes, Cp)
	cat("Drug ", i, " : ", num_term-1, "\n")
	if(length(toWrite)>0) {
		write_addr = paste("Dream7/test/Drug",as.character(i+gene_start-1),"_rna_lasso.txt",sep="")
		write.table(toWrite[order(subset(toWrite, select=c(Cp))),],file=write_addr,sep="\t",row.names = F)
	}
}