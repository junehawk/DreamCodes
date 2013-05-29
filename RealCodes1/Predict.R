library(randomForest)
library(caret)
#library(doMC)
#registerDoMC(cores = 4)

workingDIR = "/Users/junehawk/Dream7/"
MAX = 999
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

#Load Drug response data
ds = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_Drug_Response_Training.txt", header=TRUE, row.names=1, check.names=F);
num_drug = dim(ds)[1]
pred = read.table("Dream7/DrugSensitivity1/DREAM7_DrugSensitivity1_Predictions.csv",sep=",", header=TRUE, row.names=1, check.names=F)
targetCellines = setdiff(rownames(pred), colnames(ds))
result = mat.or.vec(53,31)
rownames(result) = rownames(pred)
colnames(result) = colnames(pred)
dss = t(ds)
result[rownames(dss),] = dss
result = t(result)
rm(dss)
varFilename = "Dream7/out3/selectedVars.txt"
#Prediction with saved models
	for(i in 1:num_drug) {
		prediction = rep(100, length(targetCellines))
		names(prediction) = targetCellines

		if(i==26) {
			prediction = rep(5.478, length(targetCellines))
			names(prediction) = targetCellines
			cat(prediction, "\n")
			result[i,names(prediction)] = prediction
			next
		}

		rmseFilename = paste("Dream7/out3/Drug",i,"_ModelRMSE.txt", sep="")
		modelsFilename = paste("Dream7/out3/Drug",i,"_FinalModels.RData", sep="")
		rmseloaded = read.table(rmseFilename, check.names=F)
		rmses = t(rmseloaded)
		load(modelsFilename)
		orderedRMSE = rmses[order(rmses,decreasing=T)]


		#Load information on the models that should be modeled
		modelInfo = read.table("/Users/junehawk/Dream7/input2/ModelInfo 2.txt", header=T, sep="\t", check.names=F)
		dataToLoad = colSums(modelInfo[modelInfo$Drug==i,])
		#Load data and manipulate for training
		data = list()
		activity = read.table("/Users/junehawk/Desktop/pathologist/act_all.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t", check.names=F);
		# Remove Pathways with NA values
		activity = activity[rowMeans(is.na(activity))==0,]
		rownames(activity) = paste("PA", rownames(activity), sep="")
		# Load Pathway consistency data
		consistency = read.table("/Users/junehawk/Desktop/pathologist/con_all.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t", check.names=F);
		# Remove Pathway with NA values
		consistency = consistency[rowMeans(is.na(consistency))==0,]
		rownames(consistency) = paste("PC", rownames(consistency), sep="")
		# Load Pathway Mutation data
		mutation = read.table("/Users/junehawk/Dream7/DrugSensitivity1/mutation_data_matrix_by_gene_set.txt", header=T, stringsAsFactors=FALSE, sep="\t", check.names=F)
		rownames(mutation) = paste("PM", rownames(mutation), sep="")

		data[[length(data)+1]] = activity
		data[[length(data)+1]] = consistency
		data[[length(data)+1]] = mutation
		dataIndex = c("PA", "PC", "PM")

		for(j in 2:(length(dataToLoad)-3)) {
			if(dataToLoad[j]!=0) {
				load_addr = paste(workingDIR,"input2/Drug",as.character(i),"_",names(dataToLoad)[j],".txt",sep="")
				loadedData = read.table(load_addr, header=T, sep="\t", row.names=1, check.names=F)
				rownames(loadedData) = gsub("-","_",paste(rownames(loadedData), "_", names(dataToLoad)[j], sep=""))
				rownames(loadedData) = gsub("/","_",rownames(loadedData))
				data[[length(data)+1]] = loadedData
				dataIndex = c(dataIndex, names(dataToLoad)[j])
			}
		}
		
		set.seed(109)
		#Train each models
		num_model = dim(modelInfo[modelInfo$Drug==i,])[1]
		targetModels = modelInfo[modelInfo$Drug==i,]
		targetModels = targetModels[,-1]

#		for(j in 1:num_model) {
		j=1
		while(j<=num_model) {
			modelToUse = which(rmses == orderedRMSE[j])
			if(orderedRMSE[j]>=999) {
				j = j+length(modelToUse)
				next
			}
			cat("Drug", i, "Make Prediction using Model ", modelToUse, "\n")
			for(k in 1:length(modelToUse)) {
				model = models[[modelToUse[k]]]
				cellines = list()
				dataToSelect = targetModels[modelToUse[k],]
				dataNames = names(dataToSelect[which(dataToSelect==1)])
				for(l in 1:length(dataNames)) {
					cellines[[length(cellines)+1]] = colnames(data[[which(dataIndex==dataNames[l])]])
				}
				cellines[[length(cellines)+1]] = targetCellines
				intersected = Reduce(intersect, cellines)
				testData = NULL
				for(l in 1:length(dataNames)) {
					testData = rbind(testData, subset(data[[which(dataIndex==dataNames[l])]], select=intersected))
				}
				#eliminate NA features
				testData = testData[rowMeans(is.na(testData))==0,]

				testData = as.data.frame(t(testData))
				testData = testData[rowMeans(is.na(testData))==0,]
				testData = testData[order(rownames(testData)),]

				predicted = predict(model, testData)
				names(predicted) = rownames(testData)
				prediction[names(predicted)] = predicted
				rm(testData)
				gc()
				j = j+1
			}
		}
		cat(prediction, "\n")
		result[i,names(prediction)] = prediction
		rm(models)
		rm(data)
	}

cat(result, "\n")

result = t(result)
result2=result[order(rownames(result)),]
for(i in 1:31) {
	result2[,i] = rank(-rank(result2[,i], na.last="keep", ties.method="first"))
}
#Write results to file
write_addr = paste("Dream7/out3/DrugSensitivity1.csv")
write.table(result2,file=write_addr,sep=",", quote=F)
rm(result)
rm(ds)
rm(pred)