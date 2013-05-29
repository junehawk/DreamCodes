library(RSofia)
consistency = read.table("/Users/junehawk/Desktop/pathologist/con_all.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t");
consistency = consistency[rowMeans(is.na(consistency))==0,]
ds = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_Drug_Response_Training.txt", header=TRUE, row.names=1);
cons_cellines = colnames(consistency)
drug_cellines = colnames(ds)
cons = subset(consistency, select=intersect(drug_cellines, cons_cellines))
dss = subset(ds, select=intersect(drug_cellines, cons_cellines))
num_drug = dim(dss)[1]
PathwayName = consistency[,1]
rm(consistency)
rm(ds)

for (i in 1:num_drug) {
	binded = rbind(cons, dss[i,])
	mydata = as.data.frame(t(binded))
	mydata = mydata[rowMeans(is.na(mydata))==0,]
	rm(binded)
	theTarget = colnames(mydata)[dim(mydata)[2]]
	theFormula = as.formula(paste(theTarget," ~ . "))
	model = sofia(theFormula,data=mydata,loop_type="combined-ranking")
	p = predict(model, newdata=mydata, prediction_type="linear")
	result = cor(mydata[,dim(mydata)[2],],p,method="spearman")
	#cat("Drug ", i, ": ", result, "\n")
	cat(result, "\n")
	rm(mydata)
}