library(randomForest)
library(multicore)
library(caret)
workingDIR = "/Users/junehawk/Dream7/"
resumeIndex = 30
MAX = 99999
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

#Prediction with saved models
if(resumeIndex>1)
	for(i in 1:(resumeIndex-1)) {
		prediction = rep(100, length(targetCellines))
		names(prediction) = targetCellines

		if(i==26) {
			prediction = rep(5.478, length(targetCellines))
			names(prediction) = targetCellines
			cat(prediction, "\n")
			result[i,names(prediction)] = prediction
			next
		}

		rmseFilename = paste("Dream7/out2/Drug",i,"_ModelRMSE.txt", sep="")
		modelsFilename = paste("Dream7/out2/Drug",i,"_FinalModels.RData", sep="")
		rmseloaded = read.table(rmseFilename, check.names=F)
		rmses = t(rmseloaded)
		load(modelsFilename)
		orderedRMSE = rmses[order(rmses,decreasing=T)]


		#Load information on the models that should be modeled
		modelInfo = read.table("/Users/junehawk/Dream7/input1/ModelInfo.txt", header=T, sep="\t", check.names=F)
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
				load_addr = paste(workingDIR,"input1/Drug",as.character(i),"_",names(dataToLoad)[j],".txt",sep="")
				loadedData = read.table(load_addr, header=T, sep="\t", row.names=1, check.names=F)
				rownames(loadedData) = sub("-","_",paste(rownames(loadedData), "_", names(dataToLoad)[j], sep=""))
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
#For each drug
for(i in resumeIndex:num_drug) {
#	i=10
	
	if(i==26) {
		cat("Skipping Drug 26\n")
		next
	}
		
	#Load information on the models that should be modeled
	modelInfo = read.table("/Users/junehawk/Dream7/input1/ModelInfo.txt", header=T, sep="\t", check.names=F)
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
			load_addr = paste(workingDIR,"input1/Drug",as.character(i),"_",names(dataToLoad)[j],".txt",sep="")
			loadedData = read.table(load_addr, header=T, sep="\t", row.names=1, check.names=F)
			rownames(loadedData) = sub("-","_",paste(rownames(loadedData), "_", names(dataToLoad)[j], sep=""))
			data[[length(data)+1]] = loadedData
			dataIndex = c(dataIndex, names(dataToLoad)[j])
		}
	}
	
	set.seed(109)
	#Train each models
	num_model = dim(modelInfo[modelInfo$Drug==i,])[1]
	targetModels = modelInfo[modelInfo$Drug==i,]
	targetModels = targetModels[,-1]
	rmses = seq(num_model)
	models = list()
	for(j in 1:num_model) {
	#j=1
		cellines = list()
		dataToSelect = targetModels[j,]
		dataNames = names(dataToSelect[which(dataToSelect==1)])
		for(k in 1:length(dataNames)) {
			cellines[[length(cellines)+1]] = colnames(data[[which(dataIndex==dataNames[k])]])
		}
		cellines[[length(cellines)+1]] = colnames(ds)
		intersected = Reduce(intersect, cellines)
		trainingData = NULL
		for(k in 1:length(dataNames)) {
			trainingData = rbind(trainingData, subset(data[[which(dataIndex==dataNames[k])]], select=intersected))
			#cat(dim(subset(data[[which(dataIndex==dataNames[k])]], select=intersected)),"\n")
		}
		trainingData = rbind(trainingData, subset(ds[i,],select=intersected))
		trainingData = as.data.frame(t(trainingData))
		trainingData = trainingData[rowMeans(is.na(trainingData))==0,]
		trainingData = trainingData[order(rownames(trainingData)),]

		if(dim(trainingData)[1]<12) {
			rmses[j] = MAX
			MAX = MAX+1
			cat("Drug", i, " Model", j, "Rejected due to small sample number: ", dim(trainingData)[1], "\n")
			next
		}
			
		if(dim(trainingData)[1]>=13) {
			mysample = trainingData[sample(1:nrow(trainingData), 3, replace=FALSE),]			
		} else {
			mysample = trainingData[sample(1:nrow(trainingData), dim(trainingData)[1]-10, replace=FALSE),]
		}
		mysample = mysample[order(rownames(mysample)),]
		tmp = rownames(trainingData) %in% rownames(mysample)
		trainData = trainingData[!tmp,]


		#mtry vector
		mtrySeq = NULL
		mtryNum = floor(dim(trainData[,-dim(trainData)[2]])[2]/3)
		while(mtryNum >= 2) {
			mtrySeq = c(mtryNum, mtrySeq)
			mtryNum = floor(mtryNum/2)
		}
		mtryNum = floor(dim(trainData[,-dim(trainData)[2]])[2]/3)*2
		while(mtryNum <= dim(trainData[,-dim(trainData)[2]])[2]) {
			mtrySeq = c(mtrySeq, mtryNum)
			mtryNum = mtryNum*2
		}
		cat("Drug", i, " Model", j, ": train data ", dim(trainData)[1], "sample data ", dim(mysample)[1], " Feature Selection...\n")
		####################################
		# Select Features RandomForest
		####################################
		RFE <- rfe(trainData[,-dim(trainData)[2]],
				trainData[,dim(trainData)[2]],
				sizes = c(seq(2,10, by=2), seq(10,300,by=10)),
				rfeControl = MyRFEcontrol,
					method='randomForest',
					tuneGrid = expand.grid(.mtry=mtrySeq),
					trControl = MyTrainControl)

		NewVars <- RFE$optVariables
	#	RFE
	#	plot(RFE)

		theTarget = colnames(trainData)[dim(trainData)[2]]
		FL <- as.formula(paste(theTarget, " ~ ", paste(NewVars, collapse= "+"))) #RFE

		cat("Drug", i, " Model", j, ": train data ", dim(trainData)[1], "sample data ", dim(mysample)[1], " Model Training...\n")
		model <- train(FL,data=trainData,method="rf",
			tuneGrid = expand.grid(.mtry=mtrySeq),
			trControl=MyTrainControl)

		predicted = predict(model, mysample[,-dim(mysample)[2]])
		predictedSeq = trainingData
		predictedSeq[tmp,dim(trainingData)[2]] = predicted

		cat("Drug", i, "\t", 
			cor(predict(model), trainData[,dim(trainData)[2]], method="spearman"), "\t", 
			cor(predict(model), trainData[,dim(trainData)[2]], method="kendall"), "\t", 
			RMSE(predicted, mysample[,dim(mysample)[2]]), "\t",
			cor(predictedSeq[,dim(trainingData)[2]], trainingData[,dim(trainingData)[2]], method="spearman"), "\t",
			cor(predictedSeq[,dim(trainingData)[2]], trainingData[,dim(trainingData)[2]], method="kendall"), "\n"
		)
		rmses[j] = RMSE(predicted, mysample[,dim(mysample)[2]])
		models[[length(models)+1]] = model
		rm(trainingData)
		rm(mysample)
		rm(trainData)
		gc()
	}
	#Save models to file
	write_addr = paste("Dream7/out2/Drug",i,"_PreModels.RData", sep="")
	save(models, file=write_addr)
	rm(models)

	#Train Models with all samples
	models = list()
	for(j in 1:num_model) {
		cellines = list()
		dataToSelect = targetModels[j,]
		dataNames = names(dataToSelect[which(dataToSelect==1)])
		for(k in 1:length(dataNames)) {
			cellines[[length(cellines)+1]] = colnames(data[[which(dataIndex==dataNames[k])]])
		}
		cellines[[length(cellines)+1]] = colnames(ds)
		intersected = Reduce(intersect, cellines)
		trainingData = NULL
		for(k in 1:length(dataNames)) {
			trainingData = rbind(trainingData, subset(data[[which(dataIndex==dataNames[k])]], select=intersected))
			#cat(dim(subset(data[[which(dataIndex==dataNames[k])]], select=intersected)),"\n")
		}
		trainingData = rbind(trainingData, subset(ds[i,],select=intersected))
		trainingData = as.data.frame(t(trainingData))
		trainingData = trainingData[rowMeans(is.na(trainingData))==0,]
		trainingData = trainingData[order(rownames(trainingData)),]
		#mtry vector
		mtrySeq = NULL
		mtryNum = floor(dim(trainingData[,-dim(trainingData)[2]])[2]/3)
		while(mtryNum >= 2) {
			mtrySeq = c(mtryNum, mtrySeq)
			mtryNum = floor(mtryNum/2)
		}
		mtryNum = floor(dim(trainingData[,-dim(trainingData)[2]])[2]/3)*2
		while(mtryNum <= dim(trainingData[,-dim(trainingData)[2]])[2]) {
			mtrySeq = c(mtrySeq, mtryNum)
			mtryNum = mtryNum*2
		}
	
		cat("Drug", i, " Model", j, "Final Feature Selection...\n")
		RFE <- rfe(trainingData[,-dim(trainingData)[2]],
				trainingData[,dim(trainingData)[2]],
				sizes = c(seq(2,10, by=2), seq(10,300,by=10)),
				rfeControl = MyRFEcontrol,
					method='randomForest',
					tuneGrid = expand.grid(.mtry=mtrySeq),
					trControl = MyTrainControl)

		NewVars <- RFE$optVariables
		#	RFE
		#	plot(RFE)

		theTarget = colnames(trainingData)[dim(trainingData)[2]]
		FL <- as.formula(paste(theTarget, " ~ ", paste(NewVars, collapse= "+"))) #RFE

		cat("Drug", i, " Model", j, "Final Model Training...\n")
		model <- train(FL,data=trainingData,method="rf",
			tuneGrid = expand.grid(.mtry=mtrySeq),
			trControl=MyTrainControl)
		models[[length(models)+1]] = model
		rm(trainingData)
		gc()
	}


	#Predict sensitivity by using each model
	prediction = rep(100, length(targetCellines))
	names(prediction) = targetCellines

	orderedRMSE = rmses[order(rmses,decreasing=T)]

#	for(j in 1:num_model) {
	j=1
	while(j<=num_model) {
		modelToUse = which(rmses == orderedRMSE[j])
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
	result[i,names(prediction)] = prediction
	#Save RMSE values and Models to file
	write_addr = paste("Dream7/out2/Drug",i,"_ModelRMSE.txt", sep="")
	write.table(rmses, file=write_addr)
	
	write_addr = paste("Dream7/out2/Drug",i,"_FinalModels.RData", sep="")
	save(models, file=write_addr)
	rm(models)
	
	#Free memory
	rm(data)

}

result2=result[order(rownames(result)),]
for(i in 1:31) {
	result2[,i] = rank(-rank(result2[,i], na.last="keep", ties.method="first"))
}
#Write results to file
write_addr = paste("Dream7/out2/DrugSensitivity1.csv")
write.table(result2,file=write_addr,sep=",", quote=F)
rm(result)
rm(ds)
rm(pred)