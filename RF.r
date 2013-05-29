library(randomForest)
library(multicore)
library(caret)
#library(doMC)
# Load Pathway activity data
activity = read.table("/Users/junehawk/Desktop/pathologist/act_all.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t");
# Remove Pathways with NA values
activity = activity[rowMeans(is.na(activity))==0,]
# Load Pathway consistency data
consistency = read.table("/Users/junehawk/Desktop/pathologist/con_all.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t");
# Remove Pathway with NA values
consistency = consistency[rowMeans(is.na(consistency))==0,]
# Load Pathway Mutation dataa
mutation = read.table("/Users/junehawk/Dream7/DrugSensitivity1/mutation_data_matrix_by_gene_set.txt", header=T, stringsAsFactors=FALSE, sep="\t")
# Load Drug response data
ds = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_Drug_Response_Training.txt", header=TRUE, row.names=1);
# Extract intersects only
act_cellines = colnames(activity)
con_cellines = colnames(consistency)
mut_cellines = colnames(mutation)
drug_cellines = colnames(ds)
act = subset(activity, select=intersect(intersect(drug_cellines, act_cellines), mut_cellines))
row.names(act) = paste("PA", row.names(act), sep="")
cons = subset(consistency, select=intersect(intersect(drug_cellines, act_cellines), mut_cellines))
row.names(cons) = paste("PC", row.names(cons), sep="")
muts = subset(mutation, select=intersect(intersect(drug_cellines, act_cellines), mut_cellines))
row.names(muts) = paste("PM", row.names(muts), sep="")
dss = subset(ds, select=intersect(intersect(drug_cellines, act_cellines), mut_cellines))

# Number of Drugs
num_drug = dim(dss)[1]

# Stores pathway names
PAName = activity[,1]
PCName = consistency[,1]
PMName = mutation[,1]
rm(activity)
rm(consistency)
rm(mutation)
rm(ds)
####################################
# RFE parameters
####################################
MyRFEcontrol <- rfeControl(
		functions = rfFuncs,
		method = "boot",
		number = 25,
		rerank = FALSE,
		returnResamp = "final",
		saveDetails = FALSE,
		verbose = F)
####################################
# Training parameters
####################################
MyTrainControl=trainControl(
		method = "boot",
		number=25,
		returnResamp = "all",
		selectionFunction = "oneSE"
		)

####################################
# Setup Multicore
####################################
#source:
#http://www.r-bloggers.com/feature-selection-using-the-caret-package/
if ( require("multicore", quietly = F, warn.conflicts = FALSE) ) {
	MyRFEcontrol$workers <- multicore:::detectCores()
	MyRFEcontrol$computeFunction <- mclapply
	MyRFEcontrol$computeArgs <- list(mc.preschedule = FALSE, mc.set.seed = FALSE)

	MyTrainControl$workers <- multicore:::detectCores()
	MyTrainControl$computeFunction <- mclapply
	MyTrainControl$computeArgs <- list(mc.preschedule = FALSE, mc.set.seed = FALSE)
}
#registerDoMC(cores = 6)

set.seed(109)

cat("Drug\tTrainSpearman\tTrainKendall\ttestRMSE\ttestSpearman\tTestKendall")
for (i in 1:num_drug) {
	binded = rbind(act, cons, muts, dss[i,])
	mydata = as.data.frame(t(binded))
	mydata = mydata[rowMeans(is.na(mydata))==0,]
	mydata = mydata[order(rownames(mydata)),]
	rm(binded)
	mysample = mydata[sample(1:nrow(mydata), 3,
  	 replace=FALSE),]
	mysample = mysample[order(rownames(mysample)),]
	tmp = rownames(mydata) %in% rownames(mysample)
	trainData = mydata[!tmp,]
	####################################
	# Select Features RandomForest
	####################################
	RFE <- rfe(trainData[,-dim(trainData)[2]],
			trainData[,dim(trainData)[2]],
			sizes = c(seq(2,10, by=2), seq(10,300,by=10)),
			rfeControl = MyRFEcontrol,
				method='randomForest',
				tuneGrid = expand.grid(.mtry=c(1384, 692, 346, 173, 86, 43, 21, 10, 5, 2)),
				trControl = MyTrainControl)

	NewVars <- RFE$optVariables
#	RFE
#	plot(RFE)

	theTarget = colnames(trainData)[dim(trainData)[2]]
	FL <- as.formula(paste(theTarget, " ~ ", paste(NewVars, collapse= "+"))) #RFE

	model <- train(FL,data=trainData,method="rf",
		tuneGrid = expand.grid(.mtry=c(1384, 692, 346, 173, 86, 43, 21, 10, 5, 2)),
		trControl=MyTrainControl)

	predicted = predict(model, mysample[,-dim(mysample)[2]])
	predictedSeq = mydata
	predictedSeq[tmp,dim(mydata)[2]] = predicted

	cat("Drug", i, "\t", 
		cor(predict(model), trainData[,dim(trainData)[2]], method="spearman"), "\t", 
		cor(predict(model), trainData[,dim(trainData)[2]], method="kendall"), "\t", 
		RMSE(predicted, mysample[,dim(mysample)[2]]), "\t",
		cor(predictedSeq[,dim(mydata)[2]], mydata[,dim(mydata)[2]], method="spearman"), "\t",
		cor(predictedSeq[,dim(mydata)[2]], mydata[,dim(mydata)[2]], method="kendall"), "\n"
		)

	rm(mydata)
}