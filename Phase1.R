getwd()
setwd("C:/Users/Jeremy/Desktop/final_dem/")
library(caret)
library(doParallel)
library(ISLR)
library(ADNIMERGE)
library(randomForest)
library(corrplot)
library(DMwR)
library(CORElearn)
library(pROC)
library(knitr)
library(gridExtra)
library(RColorBrewer)
library(AppliedPredictiveModeling)
library(xgboost)
library(dplyr)

total_time_start <- proc.time()
#Start parallel process with half the number of cores
cl <- makeCluster((detectCores() -1))
registerDoParallel(cl)


#============loading adnimerge data==========
df <- read.csv('adnimergedata.csv',header = T)

#============loading adnimerge data==========
df <- adnimerge

#=====Selecting only distinct baseline variables======
df <- df[,c(1,9:14,60, 65:109)]
df$RID
df <- df %>% distinct(RID, .keep_all = TRUE)
df <- subset(df, select = -c(RID,FLDSTRENG.bl,FSVERSION.bl))
summary(df)

#removing NA values for diagnosis
sum(is.na(df$DX))
df <- df[!is.na(df$DX),]
table(df$DX)

df$ABETA.bl[df$ABETA.bl == ">1700"] <- "1701"
df$ABETA.bl[df$ABETA.bl == "<200"] <- "199"
df$PTAU.bl[df$PTAU.bl == ">120"] <- "121"
df$PTAU.bl[df$PTAU.bl == "<8"] <- "7"
df$TAU.bl[df$TAU.bl == ">1300"] <- "1301"
df$TAU.bl[df$TAU.bl == "<80"] <- "79"

#converting data types
df$DX <- as.factor(df$DX)
df$PTGENDER <- as.factor(df$PTGENDER)
df$PTMARRY <- as.factor(df$PTMARRY)
df$ABETA.bl <- as.numeric(df$ABETA.bl)
df$TAU.bl <- as.numeric(df$TAU.bl)
df$PTAU.bl <- as.numeric(df$PTAU.bl)


# Shuffle data and creating train and test datasets
data<-df[sample(nrow(df)),]
split <- createDataPartition(data$DX, p = 0.75)[[1]]
train_data <- data[split,]
test_data<- data[-split,]


#==========================Training====================================

#Information gain as feature selection
# Rank the attributes in the training data using information gain
IGain<-attrEval(DX ~ .,data=train_data,estimator="InfGain")
sort(IGain, decreasing=T)

#Selecting features >= 0.01
numberOfAttributes<- length(IGain[IGain >=0.01])
ig_predictors<-names(sort(IGain, decreasing=T)[1:numberOfAttributes])
ig_predictors

#subsetting training data with only predictors from feature selection and diagnosis 
train_data <- train_data[c('DX',ig_predictors)]


#######################################################################################
####################### Splitting to use MCI vs Dementia Smote ########################
#######################################################################################

mvsd_train_data<- train_data[train_data$DX !='CN',]
table(mvsd_train_data$DX)
cn_train_data <- train_data[train_data$DX == 'CN',]


# Balance MCI and Dementia classes in the training data using smote
mvsd_train_data$DX <- factor(mvsd_train_data$DX)
train_data <- SMOTE(DX ~., data= mvsd_train_data, k=5, perc.over = 90)  
train_data <- rbind(train_data,cn_train_data)
table(train_data$DX)

# Impute the missing values in the training data using random forest
train_data <- rfImpute(DX ~ ., train_data)

# dummification of training data for nominal values
dummies <- dummyVars(DX~., data = train_data)
data_numeric <-predict(dummies, train_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(train_data$DX,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
train_data<-data_numeric 


#Shuffle train_data after smote
train_data<-train_data[sample(nrow(train_data)),]

#######################################################################################################
######################### Testing 4 models (gbm, rf, svm, knn) with Information Gain ##################
#######################################################################################################
## For accuracy, Kappa, the area under the ROC curve, sensitivity and specificity
fiveStats <- function(...) c(multiClassSummary(...),
                             defaultSummary(...))

fitControl <- trainControl(method = 'cv',
                           number = 10,
                           classProbs = TRUE,
                           summaryFunction = fiveStats,
                           verboseIter = TRUE)

#========================================X Gradient Boosting===========================================
#starting a clock to test how long the model takes
strt <- proc.time()

bl_gbmGrid <- expand.grid(nrounds = c(1, 10,100),
                          max_depth = c(1, 4,6),
                          eta = c(.1, .4),
                          gamma = 0,
                          colsample_bytree = .7,
                          min_child_weight = 1,
                          subsample = c(.5,.8,1))

bl_gbmFit <-  caret::train (x=train_data [,-1],
                            y=train_data[,1],
                            method = "xgbTree",
                            tuneGrid = bl_gbmGrid,
                            trControl = fitControl,
                            preProcess= c('zv','nzv'),
                            tuneLength =10,
                            verbose = FALSE,
                            metric='logLoss')

#ending the clock and finding the difference in time
end <- proc.time()
gbm_time <- end - strt
gbm_time

bl_gbmFit

plot(varImp(bl_gbmFit))


#=========================================== Testing =========================================

#subsetting test data with only predictors from feature selection 
test_data <- test_data[c('DX',ig_predictors)]

#Imputing missing values in test data
test_data <- rfImpute(DX ~ ., test_data)

#dummification of test data
dummies <- dummyVars(DX~., data = test_data)
data_numeric <-predict(dummies, test_data)
data_numeric <- as.data.frame(data_numeric)
data_numeric <-data.frame(test_data$DX,data_numeric)
names(data_numeric )[1] <- "Diagnosis"
test_data<-data_numeric 


bl_gbm_results <- confusionMatrix(predict(bl_gbmFit, test_data), test_data[,1])
bl_gbm_results

prob_results <- predict(bl_gbmFit, test_data, type = "prob")

bl_gbm_ROC <- multiclass.roc(test_data[,1], prob_results[,1],levels = levels(test_data[,1]))
bl_gbm_ROC



#===================================================Random Forest==================================================
strt <- proc.time()
rfFit <-  caret::train (x=train_data[,-1],
                        y=train_data[,1],
                        method = "rf",
                        trControl = fitControl,
                        ntrees = 500,
                        preProcess= c('zv', 'nzv', 'center', 'scale'),
                        tuneLength =10,
                        metric='logLoss')
end <- proc.time()
rf_time <- end - strt
rf_time

rfFit

rf_results <- confusionMatrix(predict(rfFit, test_data), test_data[,1])
rf_results

prob_results <- predict(rfFit, test_data, type = 'prob')

rf_ROC <- multiclass.roc(test_data[,1], prob_results[,1],levels = levels(test_data[,1]))
rf_ROC


#==============================Support vector machine radial kernel=========================
strt <- proc.time()
svmrFit <-  caret::train(Diagnosis~.,
                         data = train_data,
                         metric = 'logLoss',
                         method = "svmRadial",
                         tuneLength = 12,
                         preProcess= c('zv','nzv', 'center', 'scale'),
                         trControl = fitControl)
end <- proc.time()
svm_time <- end - strt
svm_time

svmrFit 

svmresults <- confusionMatrix(predict(svmrFit, test_data), test_data[,1])
svmresults

# compute the performance using probabilities.
prob_results <- predict(svmrFit, test_data, type = "prob")

#ROC curve
svmrROC <- multiclass.roc(test_data[,1], prob_results[,1],levels = levels(test_data[,1]))
svmrROC



#================================== K-Nearest Neighbor ===========================================
strt <- proc.time()

knn_dat_train <- train_data
knn_dat_test <- test_data

indx <- sapply(knn_dat_train[,-1], is.factor)
knn_dat_train[indx] <- lapply(knn_dat_train[indx], function(x) as.numeric(x))

indx <- sapply(knn_dat_test[,-1], is.factor)
knn_dat_test[indx] <- lapply(knn_dat_test[indx], function(x) as.numeric(x))


knnFit <-  caret::train (x=knn_dat_train [,-1],
                         y=knn_dat_train[,1],
                         method = "knn",
                         trControl = fitControl,
                         tuneLength =10,
                         preProcess= c('zv', 'nzv', 'center', 'scale'),
                         metric='ROC')

end <- proc.time()
knn_time <- end - strt
knn_time

knnFit

knn_results <- confusionMatrix(predict(knnFit, knn_dat_test), knn_dat_test[,1])
knn_results

prob_results <- predict(knnFit, knn_dat_test, type = "prob")

knn_ROC <- multiclass.roc(knn_dat_test[,1], prob_results[,1],levels = levels(knn_dat_test[,1]))
knn_ROC

total_time_end <- proc.time()
total_time_spent <- total_time_end - total_time_start



####################################################################################################
####################### Testing model prediction ReliefF Information gain ##########################
####################################################################################################

# Using the same splits as previous models for comparative purposes
train_data <- data[split,]
test_data<- data[-split,]


# Feature selection via reliefF method with 1000 permutations and .99 statistical significance
feat_selection_start <- proc.time()
perm <- permuteRelief(x = train_data[,names(train_data)!='DX'],
                      y = train_data$DX,
                      nperm= 1000,
                      estimator = 'ReliefFequalK')

pred.95<-names(sort(abs(perm$standardized[which(abs(perm$standardized)>=1.96)]), decreasing=T))
predictors<- pred.95


feat_selection_end <- proc.time()
feat_selection_total <- feat_selection_end - feat_selection_start


#subsetting training data with only predictors from feature selection and diagnosis 
train_data <- train_data[c('DX',predictors)]


#######################################################################################
####################### Splitting to use MCI vs Dementia Smote ########################
#######################################################################################

mvsd_train_data<- train_data[train_data$DX !='CN',]
table(mvsd_train_data$DX)
cn_train_data <- train_data[train_data$DX == 'CN',]


# Balance MCI and Dementia classes in the training data using smote
mvsd_train_data$DX <- factor(mvsd_train_data$DX)
train_data <- SMOTE(DX ~., data= mvsd_train_data, k=5, perc.over = 90)  
train_data <- rbind(train_data,cn_train_data)
table(train_data$DX)

# Impute the missing values in the training data using random forest
train_data <- rfImpute(DX ~ ., train_data)

# dummification of training data for nominal values
dummies <- dummyVars(DX~., data = train_data)
data_numeric <-predict(dummies, train_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(train_data$DX,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
train_data<-data_numeric 


#Shuffle train_data after smote
train_data<-train_data[sample(nrow(train_data)),]


## For accuracy, Kappa, the area under the ROC curve, sensitivity and specificity
fiveStats <- function(...) c(multiClassSummary(...),
                             defaultSummary(...))

fitControl <- trainControl(method = 'cv',
                           number = 10,
                           classProbs = TRUE,
                           summaryFunction = fiveStats,
                           verboseIter = TRUE)

#========================================X Gradient Boosting===========================================
#starting a clock to test how long the model takes
strt <- proc.time()

bl_gbmGrid <- expand.grid(nrounds = c(1, 10,100),
                          max_depth = c(1, 4,6),
                          eta = c(.1, .4),
                          gamma = 0,
                          colsample_bytree = .7,
                          min_child_weight = 1,
                          subsample = c(.5,.8,1))

relF_bl_gbmFit <-  caret::train (x=train_data [,-1],
                                 y=train_data[,1],
                                 method = "xgbTree",
                                 tuneGrid = bl_gbmGrid,
                                 trControl = fitControl,
                                 preProcess= c('zv','nzv'),
                                 tuneLength =10,
                                 verbose = FALSE,
                                 metric='logLoss')

#ending the clock and finding the difference in time
end <- proc.time()
relF_gbm_time <- end - strt
relF_gbm_time
relF_bl_gbmFit



#=========================================== Testing =========================================

#subsetting test data with only predictors from feature selection 
test_data <- test_data[c('DX',predictors)]

#Imputing missing values in test data
test_data <- rfImpute(DX ~ ., test_data)

#dummification of test data
dummies <- dummyVars(DX~., data = test_data)
data_numeric <-predict(dummies, test_data)
data_numeric <- as.data.frame(data_numeric)
data_numeric <-data.frame(test_data$DX,data_numeric)
names(data_numeric )[1] <- "Diagnosis"
test_data<-data_numeric 


relF_bl_gbm_results <- confusionMatrix(predict(relF_bl_gbmFit, test_data), test_data[,1])
relF_bl_gbm_results

relF_prob_results <- predict(relF_bl_gbmFit, test_data, type = "prob")

relF_bl_gbm_ROC <- multiclass.roc(test_data[,1], relF_prob_results[,1],levels = levels(test_data[,1]))
relF_bl_gbm_ROC



#===================================================Random Forest==================================================
strt <- proc.time()
relF_rfFit <-  caret::train (x=train_data[,-1],
                             y=train_data[,1],
                             method = "rf",
                             trControl = fitControl,
                             ntrees = 500,
                             preProcess= c('zv', 'nzv', 'center', 'scale'),
                             tuneLength =10,
                             metric='logLoss')

end <- proc.time()
relF_rf_time <- end - strt
relF_rf_time

relF_rfFit

relF_rf_results <- confusionMatrix(predict(relF_rfFit, test_data), test_data[,1])
relF_rf_results

relF_prob_results <- predict(relF_rfFit, test_data, type = 'prob')

relF_rf_ROC <- multiclass.roc(test_data[,1], relF_prob_results[,1],levels = levels(test_data[,1]))
relF_rf_ROC


#==============================Support vector machine radial kernel=========================
strt <- proc.time()
relF_svmrFit <-  caret::train(Diagnosis~.,
                              data = train_data,
                              metric = 'logLoss',
                              method = "svmRadial",
                              tuneLength = 12,
                              preProcess= c('zv','nzv', 'center', 'scale'),
                              trControl = fitControl)
end <- proc.time()
relF_svm_time <- end - strt
relF_svm_time

relF_svmrFit 

relF_svmresults <- confusionMatrix(predict(relF_svmrFit, test_data), test_data[,1])
relF_svmresults

# compute the performance using probabilities.
relF_prob_results <- predict(relF_svmrFit, test_data, type = "prob")

#ROC curve
relF_svmrROC <- multiclass.roc(test_data[,1], relF_prob_results[,1],levels = levels(test_data[,1]))
relF_svmrROC



#================================== K-Nearest Neighbor ===========================================
strt <- proc.time()

knn_dat_train <- train_data
knn_dat_test <- test_data

indx <- sapply(knn_dat_train[,-1], is.factor)
knn_dat_train[indx] <- lapply(knn_dat_train[indx], function(x) as.numeric(x))

indx <- sapply(knn_dat_test[,-1], is.factor)
knn_dat_test[indx] <- lapply(knn_dat_test[indx], function(x) as.numeric(x))


relF_knnFit <-  caret::train (x=knn_dat_train[,-1],
                              y=knn_dat_train[,1],
                              method = "knn",
                              trControl = fitControl,
                              tuneLength =10,
                              preProcess= c('zv', 'nzv', 'center', 'scale'),
                              metric='logLoss')

end <- proc.time()
relF_knn_time <- end - strt
relF_knn_time

relF_knnFit

relF_knn_results <- confusionMatrix(predict(relF_knnFit, knn_dat_test), knn_dat_test[,1])
relF_knn_results

relF_prob_results <- predict(relF_knnFit, knn_dat_test, type = "prob")

relF_knn_ROC <- multiclass.roc(knn_dat_test[,1], relF_prob_results[,1],levels = levels(knn_dat_test[,1]))
relF_knn_ROC

total_time_end <- proc.time()
total_time_spent <- total_time_end - total_time_start


#================Recording the Results =======================
IG_results <- NULL
perm_results <- NULL

# Information Gain Results
gbm_accuracy <- bl_gbm_results$overall[1]
gbm_kappa <- bl_gbm_results$overall[2]
gbm_AUC <- bl_gbm_ROC$auc[1]
gbm_avg_sens <- (bl_gbm_results[[4]][1]+bl_gbm_results[[4]][2]+
                        bl_gbm_results[[4]][3])/3
gbm_avg_spec <- (bl_gbm_results[[4]][4]+bl_gbm_results[[4]][5]+
                   bl_gbm_results[[4]][6])/3


rf_accuracy <- rf_results$overall[1]
rf_kappa <- rf_results$overall[2]
rf_AUC <- rf_ROC$auc[1]
rf_avg_sens <- (rf_results[[4]][1]+rf_results[[4]][2]+
                  rf_results[[4]][3])/3
rf_avg_spec <- (rf_results[[4]][4]+rf_results[[4]][5]+
                  rf_results[[4]][6])/3



svm_accuracy <- svmresults$overall[1]
svm_kappa <- svmresults$overall[2]
svm_AUC <- svmrROC$auc[1]
svm_avg_sens <- (svmresults[[4]][1]+svmresults[[4]][2]+
                   svmresults[[4]][3])/3
svm_avg_spec <- (svmresults[[4]][4]+svmresults[[4]][5]+
                   svmresults[[4]][6])/3


knn_accuracy <- knn_results$overall[1]
knn_kappa <- knn_results$overall[2]
knn_AUC <- knn_ROC$auc[1]
knn_avg_sens <- (knn_results[[4]][1]+knn_results[[4]][2]+
                   knn_results[[4]][3])/3
knn_avg_spec <- (knn_results[[4]][4]+knn_results[[4]][5]+
                   knn_results[[4]][6])/3


IG_results$accuracy <- rbind(gbm_accuracy,rf_accuracy,svm_accuracy,knn_accuracy)

IG_results$avg_sens <- rbind(gbm_avg_sens,rf_avg_sens,svm_avg_sens,knn_avg_sens)

IG_results$avg_spec <- rbind(gbm_avg_spec,rf_avg_spec,svm_avg_spec,knn_avg_spec)

IG_results$kappa <- rbind(gbm_kappa,rf_kappa,svm_kappa,knn_kappa)

IG_results$AUC <- rbind(gbm_AUC,rf_AUC,svm_AUC,knn_AUC)

IG_results <- data.frame(IG_results)
rownames(IG_results)<- c('gbm','rf','svm','knn')
IG_results


# reliF results
#================Recording the Results =======================

# Information Gain Results
perm_gbm_accuracy <- relF_bl_gbm_results$overall[1]
perm_gbm_kappa <- relF_bl_gbm_results$overall[2]
perm_gbm_AUC <- relF_bl_gbm_ROC$auc[1]
perm_gbm_avg_sens <- (relF_bl_gbm_results[[4]][1]+relF_bl_gbm_results[[4]][2]+
                        relF_bl_gbm_results[[4]][3])/3
perm_gbm_avg_spec <- (relF_bl_gbm_results[[4]][4]+relF_bl_gbm_results[[4]][5]+
                        relF_bl_gbm_results[[4]][6])/3


perm_rf_accuracy <- relF_rf_results$overall[1]
perm_rf_kappa <- relF_rf_results$overall[2]
perm_rf_AUC <- relF_rf_ROC$auc[1]
perm_rf_avg_sens <- (relF_rf_results[[4]][1]+relF_rf_results[[4]][2]+
                       relF_rf_results[[4]][3])/3
perm_rf_avg_spec <- (relF_rf_results[[4]][4]+relF_rf_results[[4]][5]+
                       relF_rf_results[[4]][6])/3



perm_svm_accuracy <- relF_svmresults$overall[1]
perm_svm_kappa <- relF_svmresults$overall[2]
perm_svm_AUC <- relF_svmrROC$auc[1]
perm_svm_avg_sens <- (relF_svmresults[[4]][1]+relF_svmresults[[4]][2]+
                        relF_svmresults[[4]][3])/3
perm_svm_avg_spec <- (relF_svmresults[[4]][4]+relF_svmresults[[4]][5]+
                        relF_svmresults[[4]][6])/3


perm_knn_accuracy <- relF_knn_results$overall[1]
perm_knn_kappa <- relF_knn_results$overall[2]
perm_knn_AUC <- relF_knn_ROC$auc[1]
perm_knn_avg_sens <- (relF_knn_results[[4]][1]+relF_knn_results[[4]][2]+
                        relF_knn_results[[4]][3])/3
perm_knn_avg_spec <- (relF_knn_results[[4]][4]+relF_knn_results[[4]][5]+
                        relF_knn_results[[4]][6])/3


perm_results$accuracy <- rbind(perm_gbm_accuracy,perm_rf_accuracy,perm_svm_accuracy,perm_knn_accuracy)

perm_results$avg_sens <- rbind(perm_gbm_avg_sens,perm_rf_avg_sens,perm_svm_avg_sens,perm_knn_avg_sens)

perm_results$avg_spec <- rbind(perm_gbm_avg_spec,perm_rf_avg_spec,perm_svm_avg_spec,perm_knn_avg_spec)

perm_results$kappa <- rbind(perm_gbm_kappa,perm_rf_kappa,perm_svm_kappa,perm_knn_kappa)

perm_results$AUC <- rbind(perm_gbm_AUC,perm_rf_AUC,perm_svm_AUC,perm_knn_AUC)

perm_results <- data.frame(perm_results)
rownames(perm_results)<- c('gbm','rf','svm','knn')
perm_results


save.image(file = 'Phase1.RData')