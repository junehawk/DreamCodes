library(elasticnet)

rnaseq = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_RNAseq_quantification.txt", header=TRUE, row.names=2, stringsAsFactors=FALSE);
ds = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_Drug_Response_Training.txt", header=TRUE, row.names=1);
rna_cellines = colnames(rnaseq)
drug_cellines = colnames(ds)
rna = subset(rnaseq, select=intersect(drug_cellines, rna_cellines))
dss = subset(ds, select=intersect(drug_cellines, rna_cellines))

FDR = dim(rna)[1]
num_rna = dim(rna)[1]
num_drug = dim(dss)[1]
#num_drug = 1
drug_start = num_rna

rm(rnaseq)
rm(ds)


for (i in 27:num_drug) {
	cat("Drug ", i, "\n")

	y = t(dss[i,])
	y2 = subset(y, !is.na(y))
	x = t(subset(rna, select=intersect(rownames(y2), rna_cellines)))
	object <- enet(x, y2, lambda=0, normalize=FALSE)
	print(object, "\n")
	rm(y)
	rm(y2)
	rm(x)
}
