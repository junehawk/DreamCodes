library(RSofia)
library(rminer)
activity = read.table("/Users/junehawk/Desktop/pathologist/act_all.txt", header=TRUE, row.names=1, sep="\t");
activity = activity[rowMeans(is.na(activity))==0,]
ds = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_Drug_Response_Training.txt", header=TRUE, row.names=1);
act_cellines = colnames(activity)
drug_cellines = colnames(ds)
act = subset(activity, select=intersect(drug_cellines, act_cellines))
dss = subset(ds, select=intersect(drug_cellines, act_cellines))
num_drug = dim(dss)[1]

for (i in 1:num_drug) {
	binded = rbind(act, dss[i,])
	mydata = as.data.frame(t(binded))
	mydata = mydata[rowMeans(is.na(mydata))==0,]
	looSpearman = array(dim=dim(mydata)[1])
	looMse = array(dim=dim(mydata)[1])
	for (j in 1:dim(mydata)[1]) {
		testset = mydata[j,]
		trainset = mydata[-j,]

		theTarget = colnames(mydata)[dim(mydata)[2]]
		theFormula = as.formula(paste(theTarget," ~ . "))

		trainY <- trainset[[which(names(trainset)==theTarget)]]
		testY <- testset[[which(names(testset)==theTarget)]]
		#model = sofia(theFormula,data=trainset,loop_type="combined-ranking")
		model = sofia(theFormula, data=trainset, learner_type="passive-aggressive", loop_type="combined-ranking")
		p = predict(model, newdata=mydata, prediction_type="linear")
		spearman = cor(mydata[,dim(mydata)[2]],p,method="spearman")
		mse = mmetric(p,mydata[,dim(mydata)[2]],"MSE")
		looSpearman[j] = spearman
		looMse[j] = mse
	}
	cat("Drug", i, "\t", mean(looSpearman), "\t", mean(looMse), "\n")
}