rnaseq = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_RNAseq_quantification.txt", header=TRUE, row.names=2);
rna = subset(rnaseq, select=c(HCC1954,	AU565,	HCC1937,	CAMA1,	UACC812,	HCC1569,	MCF12A,	HCC38,	SUM229PE,	ZR751,	BT483,	T47D,	ZR7530,	BT549,	MDAMB231,	MDAMB453,	MCF10F,	HCC1428,	HCC1419,	MDAMB361,	HCC202,	MCF7,	MDAMB175VII,	HCC1395,	HCC1143,	HCC70,	BT474,	HCC1806,	HS578T))
ds = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_Drug_Response_Training.txt", header=TRUE, row.names=1);
dss = subset(ds, select=c(HCC1954,	AU565,	HCC1937,	CAMA1,	UACC812,	HCC1569,	MCF12A,	HCC38,	SUM229PE,	ZR751,	BT483,	T47D,	ZR7530,	BT549,	MDAMB231,	MDAMB453,	MCF10F,	HCC1428,	HCC1419,	MDAMB361,	HCC202,	MCF7,	MDAMB175VII,	HCC1395,	HCC1143,	HCC70,	BT474,	HCC1806,	HS578T))

binded = rbind(rna, dss)

mydata = as.data.frame(t(binded))


FDR = dim(rna)[1]
num_rna = dim(rna)[1]
#num_drug = dim(dss)[1]
num_drug = 1
drug_start = num_rna


for (i in 1:num_drug) {
	cat("Drug ", i, "\n")
	result_pvalue=c()
	result_gene = c()
	for (j in 1:num_rna) {
		fit = lm(mydata[,i+drug_start]~mydata[,j])
		if(dim(summary(fit)$coefficient)[1] >=2 )
			pvalue = summary(fit)$coefficient[2,4]
		else
			pvalue = 1
			
			
		if(j%%1000==0)
			cat("Gene ", j, "\t")
		
#		if(pvalue <= 0.05 / FDR){
		if(!is.na(pvalue)){
			#print(summary(fit))
			result_pvalue = append(pvalue,result_pvalue,after=0)
			result_gene = append(toString(rnaseq[j,1]),result_gene,after=0)
		}
	}
#	bestfit_index = which(result_pvalue == min(result_pvalue))
#	bestfit = lm(mydata[,i+drug_start]~mydata[,result_gene[bestfit_index]])
#	print(summary(bestfit))
#	print(result_gene[bestfit_index])
	corrected = p.adjust(result_pvalue, method="fdr", n=num_rna)
	#cat(length(result_pvalue), length(result_gene), length(corrected))
	result = cbind(result_gene, result_pvalue, corrected)
	toWrite = subset(result, result_pvalue<=0.05)
	if(length(toWrite)>0) {
		write_addr = paste("Dream7/test/Drug",as.character(i),"_rna_reg.txt",sep="")
		write.table(toWrite[order(subset(toWrite, select=c(corrected))),],file=write_addr,sep="\t",row.names = F)
	}
}