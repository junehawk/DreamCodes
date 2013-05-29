#Directory
setwd('~/wherever')

#Load Required Packages
library('caret')
library('glmnet')
library('ipred')
library('e1071')
library('caTools')

############################
# Load the Data, choose target, create train and test sets
############################

Data <- read.csv("overfitting.csv", header=TRUE)

#Choose Target
Data$Target <- as.factor(ifelse(Data$Target_Practice==1,'X1','X0'))
Data$Target_Evaluate = NULL
Data$Target_Leaderboard = NULL
Data$Target_Practice = NULL

#Order
xnames <- setdiff(names(Data),c('Target','case_id','train'))
Data <- Data[,c('Target','case_id','train',xnames)]

#Split to train and test
trainset = Data[Data$train == 1,]
testset = Data[Data$train == 0,]

#Remove unwanted columns
trainset$case_id = NULL
trainset$train = NULL

#Define Formula
FL <- as.formula(paste("Target ~ ", paste(xnames, collapse= "+")))

####################################
# RFE parameters
####################################
library(ipred)
library(e1071)

#Custom Functions
glmnetFuncs <- caretFuncs #Default caret functions

glmnetFuncs$fit <- function (x, y, first, last, ...) { #Fits a GLMNET model
    library("glmnet")
    glmnet(as.matrix(x), y, family = "multinomial", alpha = 0, lambda = 0.02)
}

glmnetFuncs$pred <- function (object, x) { #Makes predictions (in a format other cart functions recognize) from a glmnet model
	tmp <- predict(object, newx=as.matrix(x))
	tmp <- data.frame(tmp)
	names(tmp) <- sub('.s0','',names(tmp))
	tmp$pred <- ifelse(tmp[,1]>tmp[,2],names(tmp)[1],names(tmp)[2])
	tmp
}

glmnetFuncs$rank <- function (object, x, y) { #Ranks predictions, and numbers them from best to worst
	vimp <- sort(object$beta[[2]][, 1])
	vimp <- as.data.frame(vimp)
	vimp$var <- row.names(vimp)
	vimp$'Overall' <- seq(nrow(vimp),1)
	vimp
}

glmnetFuncs$summary <- function (data, lev = NULL, model = NULL) { #Computes Sens, Spec. and ROC for a given model
    if (is.character(data$obs)) {
        data$obs <- factor(data$obs, levels = lev)
	}
	if (is.character(data$pred)) {
        data$pred <- factor(data$pred, levels = lev)
	}
    twoClassSummary(data, lev = lev, model = NULL)
}

#fit <- glmnetFuncs$fit(x,y) #TEST that the functions work properly
#pred <- glmnetFuncs$pred(fit,x)
#rank <- glmnetFuncs$rank(fit)

MyRFEcontrol <- rfeControl(
		functions = glmnetFuncs,
		method = "repeatedCV",
		number = 10,
		repeats = 5,
		rerank = FALSE,
		returnResamp = "final",
		saveDetails = FALSE,
		verbose = TRUE)

####################################
# Training parameters
####################################
MyTrainControl=trainControl(
		method = "repeatedCV",
		number=10,
		repeats=5,
		returnResamp = "all",
		classProbs = TRUE,
		summaryFunction=twoClassSummary
		)

####################################
# Setup Multicore
####################################
#source:
#http://www.r-bloggers.com/feature-selection-using-the-caret-package/
if ( require("multicore", quietly = TRUE, warn.conflicts = FALSE) ) {
	MyRFEcontrol$workers <- multicore:::detectCores()
	MyRFEcontrol$computeFunction <- mclapply
	MyRFEcontrol$computeArgs <- list(mc.preschedule = FALSE, mc.set.seed = FALSE)
	
	MyTrainControl$workers <- multicore:::detectCores()
	MyTrainControl$computeFunction <- mclapply
	MyTrainControl$computeArgs <- list(mc.preschedule = FALSE, mc.set.seed = FALSE)
}

####################################
# Select Features
####################################

x <- trainset[,xnames]
y <- trainset$Target

RFE <- rfe(x,y,sizes = seq(130,160,by=1), #Optimal # of features seems to be within ~130-160
	metric='ROC',
	maximize=TRUE,
	rfeControl = MyRFEcontrol)
	
NewVars <- RFE$optVariables
FL <- as.formula(paste("Target ~ ", paste(NewVars, collapse= "+")))
RFE
plot(RFE)

####################################
# Fit a GLMNET Model
####################################

G <- expand.grid(.alpha=0,.lambda=seq(0,0.05,by=0.01))

model <- train(FL,data=trainset,method='glmnet',
	metric = "ROC",
	tuneGrid = G,
	trControl=MyTrainControl)
plot(model, metric='ROC')
predictions <- predict(model, newdata=testset, type  = "prob")
colAUC(predictions, testset$Target)

########################################
#Generate a file for submission
########################################
testID  <- testset$case_id
submit_file = cbind(testID,predictions)
submit_file <- submit_file[,c('testID','X1','X0')]
write.csv(submit_file, file="submission.csv", row.names = FALSE)

