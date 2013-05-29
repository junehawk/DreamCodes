library(elasticnet)
library(multicore)

act = read.table("/Users/junehawk/Desktop/pathologist/act_all.txt", header=TRUE, sep="\t");
ds = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_Drug_Response_Training.txt", header=TRUE, row.names=1);
act_cellines = colnames(act)
drug_cellines = colnames(ds)
activity = subset(act, select=intersect(drug_cellines, act_cellines))
dss = subset(ds, select=intersect(drug_cellines, act_cellines))

pathway_names = subset(act, select=c(Pathway.Name))
num_drug = dim(dss)[1]
gene_start = 27
rm(act)
rm(ds)

lassofit <- function(i, ddata = dss, adata = activity) {
	cat(i,"\n")
	y = t(ddata[i,])
	y2 = subset(y, !is.na(y))
	x = t(subset(adata, select=intersect(rownames(y2), colnames(adata))))
	x2 = x[,colMeans(is.na(x)) == 0]
	object = enet(x2, y2, lambda=0, normalize=FALSE)
	#if(class(result) == “try-error”)
	#	next;
	rm(y)
	rm(y2)
	rm(x)
	object
}


index = as.list(c(1:num_drug))
result = mclapply(index, lassofit)

for(i in 1:length(result)) {

	Cp = as.vector(result[[i]]$Cp)
	Cp = Cp[1:(length(Cp)-1)]
	num_term = length(result[[i]]$actions)
	genes = c()
	for(j in 1:(num_term-1)) {
		genes = append(pathway_names[names(result[[i]]$actions[[j]]),], genes, after=0)
	}
	toWrite = cbind(genes, Cp)
	cat("Drug ", i, " : ", num_term-1, "\n")
	if(length(toWrite)>0) {
		write_addr = paste("Dream7/test/Drug",as.character(i),"_activity_lasso.txt",sep="")
		write.table(toWrite[order(subset(toWrite, select=c(Cp))),],file=write_addr,sep="\t",row.names = F)
	}
}