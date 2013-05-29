library(multicore)
consistency = read.table("/Users/junehawk/Desktop/pathologist/con_all.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t");
consistency = consistency[rowMeans(is.na(consistency))==0,]
ds = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_Drug_Response_Training.txt", header=TRUE, row.names=1);
cons_cellines = colnames(consistency)
drug_cellines = colnames(ds)
cons = subset(consistency, select=intersect(drug_cellines, cons_cellines))
dss = subset(ds, select=intersect(drug_cellines, cons_cellines))
binded = rbind(cons, dss)
mydata = as.data.frame(t(binded))

FDR = dim(cons)[1]
num_path = dim(cons)[1]
num_drug = dim(dss)[1]
#num_drug = 1
drug_start = num_path

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

PathwayName = consistency[,1]

rm(consistency)
rm(ds)
rm(cons)
rm(dss)
rm(binded)

x = as.list(c(1:num_path))
for (i in 1:num_drug) {
	#cat("Drug ", i, "\n")
	

	pvalue = unlist(mclapply(x,lmpvalue,i=i))	
	corrected = p.adjust(pvalue, method="fdr", n=num_path)
	
	result = cbind(PathwayName, pvalue, corrected)
	
	toWrite = subset(result, pvalue<=0.05)

	candidate = subset(result, corrected<=0.1)
	cat("Drug ", i, " : ", dim(candidate)[1], "\n")
	rm(candidate)

	if(length(toWrite)>0) {
		write_addr = paste("Dream7/test/Drug",as.character(i),"_pathCon_reg.txt",sep="")
		write.table(toWrite[order(subset(toWrite, select=c(corrected))),],file=write_addr,sep="\t",row.names = F)
	}
	rm(pvalue)
	rm(corrected)
	rm(result)
	rm(toWrite)
}
rm(mydata)