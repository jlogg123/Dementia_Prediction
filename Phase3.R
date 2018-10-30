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
library(plyr)
library(dplyr)
library(data.table)
library(e1071)
library(doBy)

#Start parallel process with half the number of cores
cl <- makeCluster((detectCores() -1))
registerDoParallel(cl)

#total_time_start <- proc.time()

#============loading adnimerge data==========
df <- read.csv('adnimergedata.csv',header = T)

#=====================Variable Selection====================
df <- df[,c(1:62)]
names(df)

df <- subset(df, select = -c(FLDSTRENG,FSVERSION,
                             COLPROT,ORIGPROT,PTID,SITE,VISCODE,EXAMDATE,CDRSB,DX.bl))

#removing NA values for diagnosis
sum(is.na(df$DX))
df <- df[!is.na(df$DX),]
table(df$DX)


#============creating a new column which represents change in diagnosis==============
for (i in 2:nrow(df)){
  if (df$RID[i] == df$RID[i-1]){
    if ((df$DX[i] == 'CN')& (df$DX[i-1]=='CN')) {
      df$Diagnosis[i] <- "CN"
    }
    else if ((df$DX[i] == 'CN')& (df$DX[i-1]=='MCI')) {
      df$Diagnosis[i] <- "CN_MCI"
    }
    else if ((df$DX[i] == 'MCI')& (df$DX[i-1]=='CN')) {
      df$Diagnosis[i] <- "MCI_CN"
    }
    else if ((df$DX[i] == 'MCI')& (df$DX[i-1]=='MCI')) {
      df$Diagnosis[i] <- "MCI"
    }
    else if ((df$DX[i] == 'MCI')& (df$DX[i-1]=='Dementia')) {
      df$Diagnosis[i] <- "MCI_Dementia"
    }
    else if ((df$DX[i] == 'Dementia')& (df$DX[i-1]=='MCI')) {
      df$Diagnosis[i] <- "Dementia_MCI"
    }
    else if ((df$DX[i] == 'Dementia')& (df$DX[i-1]=='Dementia')) {
      df$Diagnosis[i] <- "Dementia"
    }
  }
  else if (df$RID[i] != df$RID[i-1]){
    if(df$DX[i] == 'CN'){
      df$Diagnosis[i] <- 'CN'
    }
    else if (df$DX[i] == 'MCI'){
      df$Diagnosis[i] <- 'MCI'
    }
    else if (df$DX[i] == 'Dementia'){
      df$Diagnosis[i] <- 'Dementia'
    }
  }
}

table(df$Diagnosis)

#removing DX variable
df <-df[,-50]

#assigning the first patient a diagnosis status of CN (because that was the actual diagnois)
df$Diagnosis[1] <- 'CN'

df$ABETA <- as.character(df$ABETA)
df$TAU <- as.character(df$TAU)
df$PTAU <- as.character(df$PTAU)

df$ABETA[df$ABETA == ">1700"] <- "1701"
df$ABETA[df$ABETA == "<200"] <- "199"
df$PTAU[df$PTAU == ">120"] <- "121"
df$PTAU[df$PTAU == "<8"] <- "7"
df$TAU[df$TAU == ">1300"] <- "1301"
df$TAU[df$TAU == "<80"] <- "79"

#converting data types
df$PTGENDER <- as.factor(df$PTGENDER)
df$PTMARRY <- as.factor(df$PTMARRY)
df$ABETA <- as.numeric(df$ABETA)
df$TAU <- as.numeric(df$TAU)
df$PTAU <- as.numeric(df$PTAU)

#rearranging Diagnosis column
col_order <-c( "RID",            "Diagnosis",                   "AGE",       "PTGENDER",             
               "PTEDUCAT",              "PTETHCAT",              "PTRACCAT",             
               "PTMARRY",               "APOE4",                 "FDG",                  
               "PIB",                   "AV45",                  "ABETA",                
               "TAU",                   "PTAU"  ,                "ADAS11",               
               "ADAS13" ,               "ADASQ4"  ,              "MMSE"    ,             
               "RAVLT.immediate"  ,     "RAVLT.learning"     ,   "RAVLT.forgetting"  ,   
               "RAVLT.perc.forgetting", "LDELTOTAL",             "DIGITSCOR"            ,
               "TRABSCOR",              "FAQ",                   "MOCA",                 
               "EcogPtMem",             "EcogPtLang",            "EcogPtVisspat",        
               "EcogPtPlan",            "EcogPtOrgan",           "EcogPtDivatt",         
               "EcogPtTotal",           "EcogSPMem",             "EcogSPLang",           
               "EcogSPVisspat",         "EcogSPPlan",            "EcogSPOrgan",          
               "EcogSPDivatt",          "EcogSPTotal",           "IMAGEUID",             
               "Ventricles",            "Hippocampus",           "WholeBrain",           
               "Entorhinal",            "Fusiform",              "MidTemp",              
               "ICV",                 "mPACCdigit",              "mPACCtrailsB")

df<-df[,col_order]
#==========================Feature engineering==========================

### Baseline dataframe
base <- df

### velocity dataframe
vel <- df
## adding new velocity columns
newcols <- names(vel[10:52])
newcols <- paste('vel',newcols,sep = '_')
vel[,newcols] <- NA

for(i in 10:52){
  for(r in 2:nrow(vel)){
    if(vel$RID[r] == vel$RID[r-1]){
      vel[r,i+43] <-vel[r,i]-vel[r-1,i]
    }
    else if (vel$RID[r] != vel$RID[r-1]){
      vel[r,i+43]<- NA
    }
  }
}

### Acceleration Dataframe
accel <- vel

## Adding new columns
newcols <- names(accel[10:52])
newcols <- paste('accel',newcols,sep = '_')
accel[,newcols] <- NA

for(i in 53:95){
  for(r in 2:nrow(accel)){
    if(accel$RID[r] == accel$RID[r-1]){
      accel[r,i+43] <-accel[r,i]-accel[r-1,i]
    }
    else if (accel$RID[r] != accel$RID[r-1]){
      accel[r,i+43]<- NA
    }
  }
}
accel <- accel[,-98]

### Absolute value of Accel Dataframe  
abs_accel <- accel

## Adding new columns
newcols <- names(abs_accel[10:52])
newcols <- paste('abs_accel',newcols,sep = '_')
abs_accel[,newcols] <- NA

for(i in 96:138){
  for(r in 2:nrow(abs_accel)){
    abs_accel[r,i+43] <-abs(abs_accel[r,i])
  }
}
abs_accel <- abs_accel[,-c(98,141)]

########################### Datasets version 1 using quantiles #################################

#baseline dataframe
demo <- base[,c('RID','Diagnosis','AGE','PTGENDER','PTEDUCAT','PTETHCAT','PTRACCAT','PTMARRY','APOE4')]
demo <- demo %>% distinct(RID,Diagnosis,.keep_all = T)
demo <- demo[,-c(1,2)]
base_v1 <- base[,c(1,2,10:ncol(base))]
base_v1 <- summaryBy(.~RID+Diagnosis, FUN = quantile,na.rm = T,data = base_v1)
base_v1 <- cbind(base_v1,base_v1)
base_v1 <- base_v1[,-1]
base_v1 <- base_v1[!is.na(base_v1$Diagnosis),]
base_v1$Diagnosis<- as.factor(base_v1$Diagnosis)

#velocity dataframe
vel_v1 <- vel[,c(1,2,10:ncol(vel))]
vel_v1 <- summaryBy(.~c(RID,Diagnosis), FUN = quantile,na.rm = T,data = vel_v1)
vel_v1 <- cbind(vel_v1,demo)
vel_v1 <- vel_v1[,-1]
vel_v1 <- vel_v1[!is.na(vel_v1$Diagnosis),]
vel_v1$Diagnosis<- as.factor(vel_v1$Diagnosis)

#accel dataframe
accel_v1 <- accel[,c(1,2,10:ncol(accel))]
accel_v1 <- summaryBy(.~ RID + Diagnosis, FUN = quantile,na.rm = T, data = accel_v1)
accel_v1 <- cbind(accel_v1,demo)
accel_v1 <- accel_v1[,-1]
accel_v1 <- accel_v1[!is.na(accel_v1$Diagnosis),]
accel_v1$Diagnosis<- as.factor(accel_v1$Diagnosis)


#abs_accel dataframe
abs_accel_v1 <- abs_accel[,c(1,2,10:ncol(abs_accel))]
abs_accel_v1 <- summaryBy(.~RID+Diagnosis,FUN =quantile,na.rm = T,data = abs_accel_v1)
abs_accel_v1 <- abs_accel_v1[,1:842]
abs_accel_v1 <- cbind(abs_accel_v1,demo)
abs_accel_v1 <- abs_accel_v1[,-1]
abs_accel_v1 <- abs_accel_v1[!is.na(abs_accel_v1$Diagnosis),]
abs_accel_v1$Diagnosis<- as.factor(abs_accel_v1$Diagnosis)




## Shuffle data and creating train and test datasets
base_data<-base_v1[sample(nrow(base_v1)),]
base_split <- createDataPartition(base_v1$Diagnosis, p = 0.75)[[1]]
base_train_data <- base_data[base_split,]
base_test_data<- base_data[-base_split,]

vel_data<-vel_v1[sample(nrow(vel_v1)),]
vel_split <- createDataPartition(vel_v1$Diagnosis, p = 0.75)[[1]]
vel_train_data <- vel_data[vel_split,]
vel_test_data<- vel_data[-vel_split,]

accel_data<-accel_v1[sample(nrow(accel_v1)),]
accel_split <- createDataPartition(accel_v1$Diagnosis, p = 0.75)[[1]]
accel_train_data <- accel_data[accel_split,]
accel_test_data<- accel_data[-accel_split,]

abs_accel_data<-abs_accel_v1[sample(nrow(abs_accel_v1)),]
abs_accel_split <- createDataPartition(abs_accel_v1$Diagnosis, p = 0.75)[[1]]
abs_accel_train_data <- abs_accel_data[abs_accel_split,]
abs_accel_test_data<- abs_accel_data[-abs_accel_split,]


#==========================Training====================================

#Information gain as feature selection
# Rank the attributes in the training data using information gain
#Selecting features >= 0.01
base_IGain<-attrEval(Diagnosis ~ .,data=base_train_data,estimator="InfGain")
base_numberOfAttributes<- length(base_IGain[base_IGain >=0.01])
base_ig_predictors<-names(sort(base_IGain, decreasing=T)[1:base_numberOfAttributes])
base_train_data <- base_train_data[c('Diagnosis',base_ig_predictors)]


vel_IGain<-attrEval(Diagnosis ~ .,data=vel_train_data,estimator="InfGain")
vel_numberOfAttributes<- length(vel_IGain[vel_IGain >=0.01])
vel_ig_predictors<-names(sort(vel_IGain, decreasing=T)[1:vel_numberOfAttributes])
vel_train_data <- vel_train_data[c('Diagnosis',vel_ig_predictors)]


accel_IGain<-attrEval(Diagnosis ~ .,data = accel_train_data,estimator="InfGain")
accel_numberOfAttributes<- length(accel_IGain[accel_IGain >=0.01])
accel_ig_predictors<-names(sort(accel_IGain, decreasing=T)[1:accel_numberOfAttributes])
accel_train_data <- accel_train_data[c('Diagnosis',accel_ig_predictors)]


abs_accel_IGain<-attrEval(Diagnosis ~ .,data=abs_accel_train_data,estimator="InfGain")
abs_accel_numberOfAttributes<- length(abs_accel_IGain[abs_accel_IGain >=0.01])
abs_accel_ig_predictors<-names(sort(abs_accel_IGain, decreasing=T)[1:abs_accel_numberOfAttributes])
abs_accel_train_data <- abs_accel_train_data[c('Diagnosis',abs_accel_ig_predictors)]


# Impute the missing values in the training data using random forest
base_train_data <- rfImpute(Diagnosis ~ ., base_train_data)
vel_train_data <- rfImpute(Diagnosis ~ ., vel_train_data)
accel_train_data <- rfImpute(Diagnosis ~ ., accel_train_data)
abs_accel_train_data <- rfImpute(Diagnosis ~ ., abs_accel_train_data)

# dummification of training data for nominal values
dummies <- dummyVars(Diagnosis~., data = base_train_data)
data_numeric <-predict(dummies, base_train_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(base_train_data$Diagnosis,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
base_train_data<-data_numeric 

dummies <- dummyVars(Diagnosis~., data = vel_train_data)
data_numeric <-predict(dummies, vel_train_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(vel_train_data$Diagnosis,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
vel_train_data<-data_numeric 

dummies <- dummyVars(Diagnosis~., data = accel_train_data)
data_numeric <-predict(dummies, accel_train_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(accel_train_data$Diagnosis,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
accel_train_data<-data_numeric 

dummies <- dummyVars(Diagnosis~., data = abs_accel_train_data)
data_numeric <-predict(dummies, abs_accel_train_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(abs_accel_train_data$Diagnosis,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
abs_accel_train_data<-data_numeric 


#Shuffle train_data after smote
base_train_data<-base_train_data[sample(nrow(base_train_data)),]
vel_train_data<-vel_train_data[sample(nrow(vel_train_data)),]
accel_train_data<-accel_train_data[sample(nrow(accel_train_data)),]
abs_accel_train_data<-abs_accel_train_data[sample(nrow(abs_accel_train_data)),]


#subsetting test data with only predictors from feature selection 
base_test_data <- base_test_data[c('Diagnosis',base_ig_predictors)]
vel_test_data <- vel_test_data[c('Diagnosis',vel_ig_predictors)]
accel_test_data <- accel_test_data[c('Diagnosis',accel_ig_predictors)]
abs_accel_test_data <- abs_accel_test_data[c('Diagnosis',abs_accel_ig_predictors)]

#Imputing missing values in test data
base_test_data <- rfImpute(Diagnosis ~ ., base_test_data)
vel_test_data <- rfImpute(Diagnosis ~ ., vel_test_data)
accel_test_data <- rfImpute(Diagnosis ~ ., accel_test_data)
abs_accel_test_data <- rfImpute(Diagnosis ~ ., abs_accel_test_data)

#dummification of test data
dummies <- dummyVars(Diagnosis~., data = base_test_data)
data_numeric <-predict(dummies, base_test_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(base_test_data$Diagnosis,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
base_test_data<-data_numeric 

dummies <- dummyVars(Diagnosis~., data = vel_test_data)
data_numeric <-predict(dummies, vel_test_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(vel_test_data$Diagnosis,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
vel_test_data<-data_numeric 

dummies <- dummyVars(Diagnosis~., data = accel_test_data)
data_numeric <-predict(dummies, accel_test_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(accel_test_data$Diagnosis,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
accel_test_data<-data_numeric 

dummies <- dummyVars(Diagnosis~., data = abs_accel_test_data)
data_numeric <-predict(dummies, abs_accel_test_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(abs_accel_test_data$Diagnosis,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
abs_accel_test_data<-data_numeric 


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
#strt <- proc.time()

bl_gbmGrid <- expand.grid(nrounds = c(1, 10,100),
                          max_depth = c(1, 4,6),
                          eta = c(.1, .4),
                          gamma = 0,
                          colsample_bytree = .7,
                          min_child_weight = 1,
                          subsample = c(.5,.8,1))

base_bl_gbmFit <-  caret::train (x=base_train_data [,-1],
                                 y=base_train_data[,1],
                                 method = "xgbTree",
                                 tuneGrid = bl_gbmGrid,
                                 trControl = fitControl,
                                 preProcess= c('zv','nzv'),
                                 tuneLength =10,
                                 verbose = FALSE,
                                 metric='logLoss')

vel_bl_gbmFit <-  caret::train (x=vel_train_data [,-1],
                                y=vel_train_data[,1],
                                method = "xgbTree",
                                tuneGrid = bl_gbmGrid,
                                trControl = fitControl,
                                preProcess= c('zv','nzv'),
                                tuneLength =10,
                                verbose = FALSE,
                                metric='logLoss')

accel_bl_gbmFit <-  caret::train (x=accel_train_data [,-1],
                                  y=accel_train_data[,1],
                                  method = "xgbTree",
                                  tuneGrid = bl_gbmGrid,
                                  trControl = fitControl,
                                  preProcess= c('zv','nzv'),
                                  tuneLength =10,
                                  verbose = FALSE,
                                  metric='logLoss')

abs_accel_bl_gbmFit <-  caret::train (x=abs_accel_train_data [,-1],
                                      y=abs_accel_train_data[,1],
                                      method = "xgbTree",
                                      tuneGrid = bl_gbmGrid,
                                      trControl = fitControl,
                                      preProcess= c('zv','nzv'),
                                      tuneLength =10,
                                      verbose = FALSE,
                                      metric='logLoss')


#Testing
base_bl_gbm_results <- confusionMatrix(predict(base_bl_gbmFit, base_test_data), base_test_data[,1])
base_bl_gbm_results
base_prob_results <- predict(base_bl_gbmFit, base_test_data, type = "prob")
base_bl_gbm_ROC <- multiclass.roc(base_test_data[,1], base_prob_results[,1],levels = levels(base_test_data[,1]))
base_bl_gbm_ROC

vel_bl_gbm_results <- confusionMatrix(predict(vel_bl_gbmFit, vel_test_data), vel_test_data[,1])
vel_bl_gbm_results
vel_prob_results <- predict(vel_bl_gbmFit, vel_test_data, type = "prob")
vel_bl_gbm_ROC <- multiclass.roc(vel_test_data[,1], vel_prob_results[,1],levels = levels(vel_test_data[,1]))
vel_bl_gbm_ROC

accel_bl_gbm_results <- confusionMatrix(predict(accel_bl_gbmFit, accel_test_data), accel_test_data[,1])
accel_bl_gbm_results
accel_prob_results <- predict(accel_bl_gbmFit, accel_test_data, type = "prob")
accel_bl_gbm_ROC <- multiclass.roc(accel_test_data[,1], accel_prob_results[,1],levels = levels(accel_test_data[,1]))
accel_bl_gbm_ROC

abs_accel_bl_gbm_results <- confusionMatrix(predict(abs_accel_bl_gbmFit, abs_accel_test_data), abs_accel_test_data[,1])
abs_accel_bl_gbm_results
abs_accel_prob_results <- predict(abs_accel_bl_gbmFit, abs_accel_test_data, type = "prob")
abs_accel_bl_gbm_ROC <- multiclass.roc(abs_accel_test_data[,1], abs_accel_prob_results[,1],levels = levels(abs_accel_test_data[,1]))
abs_accel_bl_gbm_ROC


#===================================================Random Forest==================================================
#strt <- proc.time()
base_rfFit <-  caret::train (x=base_train_data[,-1],
                             y=base_train_data[,1],
                             method = "rf",
                             trControl = fitControl,
                             ntrees = 500,
                             preProcess= c('zv', 'nzv', 'center', 'scale'),
                             tuneLength =10,
                             metric='logLoss')

vel_rfFit <-  caret::train (x=vel_train_data[,-1],
                            y=vel_train_data[,1],
                            method = "rf",
                            trControl = fitControl,
                            ntrees = 500,
                            preProcess= c('zv', 'nzv', 'center', 'scale'),
                            tuneLength =10,
                            metric='logLoss')

accel_rfFit <-  caret::train (x=accel_train_data[,-1],
                              y=accel_train_data[,1],
                              method = "rf",
                              trControl = fitControl,
                              ntrees = 500,
                              preProcess= c('zv', 'nzv', 'center', 'scale'),
                              tuneLength =10,
                              metric='logLoss')

abs_accel_rfFit <-  caret::train (x=abs_accel_train_data[,-1],
                                  y=abs_accel_train_data[,1],
                                  method = "rf",
                                  trControl = fitControl,
                                  ntrees = 500,
                                  preProcess= c('zv', 'nzv', 'center', 'scale'),
                                  tuneLength =10,
                                  metric='logLoss')


base_rf_results <- confusionMatrix(predict(base_rfFit, base_test_data), base_test_data[,1])
base_rf_results
base_prob_results <- predict(base_rfFit, base_test_data, type = 'prob')
base_rf_ROC <- multiclass.roc(base_test_data[,1], base_prob_results[,1],levels = levels(base_test_data[,1]))
base_rf_ROC

vel_rf_results <- confusionMatrix(predict(vel_rfFit, vel_test_data), vel_test_data[,1])
vel_rf_results
vel_prob_results <- predict(vel_rfFit, vel_test_data, type = 'prob')
vel_rf_ROC <- multiclass.roc(vel_test_data[,1], vel_prob_results[,1],levels = levels(vel_test_data[,1]))
vel_rf_ROC

accel_rf_results <- confusionMatrix(predict(accel_rfFit, accel_test_data), accel_test_data[,1])
accel_rf_results
accel_prob_results <- predict(accel_rfFit, accel_test_data, type = 'prob')
accel_rf_ROC <- multiclass.roc(accel_test_data[,1], accel_prob_results[,1],levels = levels(accel_test_data[,1]))
accel_rf_ROC

abs_accel_rf_results <- confusionMatrix(predict(abs_accel_rfFit, abs_accel_test_data), abs_accel_test_data[,1])
abs_accel_rf_results
abs_accel_prob_results <- predict(abs_accel_rfFit, abs_accel_test_data, type = 'prob')
abs_accel_rf_ROC <- multiclass.roc(abs_accel_test_data[,1], abs_accel_prob_results[,1],levels = levels(abs_accel_test_data[,1]))
abs_accel_rf_ROC


#==============================Support vector machine radial kernel=========================
#strt <- proc.time()
base_svmrFit <-  caret::train(Diagnosis~.,
                              data = base_train_data,
                              metric = 'logLoss',
                              method = "svmRadial",
                              tuneLength = 12,
                              preProcess= c('zv','nzv', 'center', 'scale'),
                              trControl = fitControl)

vel_svmrFit <-  caret::train(Diagnosis~.,
                             data = vel_train_data,
                             metric = 'logLoss',
                             method = "svmRadial",
                             tuneLength = 12,
                             preProcess= c('zv','nzv', 'center', 'scale'),
                             trControl = fitControl)

accel_svmrFit <-  caret::train(Diagnosis~.,
                               data = accel_train_data,
                               metric = 'logLoss',
                               method = "svmRadial",
                               tuneLength = 12,
                               preProcess= c('zv','nzv', 'center', 'scale'),
                               trControl = fitControl)

abs_accel_svmrFit <-  caret::train(Diagnosis~.,
                                   data = abs_accel_train_data,
                                   metric = 'logLoss',
                                   method = "svmRadial",
                                   tuneLength = 12,
                                   preProcess= c('zv','nzv', 'center', 'scale'),
                                   trControl = fitControl)
#end <- proc.time()
#svm_time <- end - strt
#svm_time
#svmrFit 

base_svmresults <- confusionMatrix(predict(base_svmrFit, base_test_data), base_test_data[,1])
base_svmresults
base_prob_results <- predict(base_svmrFit, base_test_data, type = "prob")
base_svmrROC <- multiclass.roc(base_test_data[,1], base_prob_results[,1],levels = levels(base_test_data[,1]))
base_svmrROC

vel_svmresults <- confusionMatrix(predict(vel_svmrFit, vel_test_data), vel_test_data[,1])
vel_svmresults
vel_prob_results <- predict(vel_svmrFit, vel_test_data, type = "prob")
vel_svmrROC <- multiclass.roc(vel_test_data[,1], vel_prob_results[,1],levels = levels(vel_test_data[,1]))
vel_svmrROC

accel_svmresults <- confusionMatrix(predict(accel_svmrFit, accel_test_data), accel_test_data[,1])
accel_svmresults
accel_prob_results <- predict(accel_svmrFit, accel_test_data, type = "prob")
accel_svmrROC <- multiclass.roc(accel_test_data[,1],accel_prob_results[,1],levels = levels(accel_test_data[,1]))
accel_svmrROC

abs_accel_svmresults <- confusionMatrix(predict(abs_accel_svmrFit, abs_accel_test_data), abs_accel_test_data[,1])
abs_accel_svmresults
abs_accel_prob_results <- predict(abs_accel_svmrFit, abs_accel_test_data, type = "prob")
abs_accel_svmrROC <- multiclass.roc(abs_accel_test_data[,1], abs_accel_prob_results[,1],levels = levels(abs_accel_test_data[,1]))
abs_accel_svmrROC

#================================== K-Nearest Neighbor ===========================================
#strt <- proc.time()

base_knn_dat_train <- base_train_data
base_knn_dat_test <- base_test_data
base_indx <- sapply(base_knn_dat_train[,-1], is.factor)
base_knn_dat_train[base_indx] <- lapply(base_knn_dat_train[base_indx], function(x) as.numeric(x))
base_indx <- sapply(base_knn_dat_test[,-1], is.factor)
base_knn_dat_test[base_indx] <- lapply(base_knn_dat_test[base_indx], function(x) as.numeric(x))
base_knnFit <-  caret::train (x=base_knn_dat_train [,-1],
                              y=base_knn_dat_train[,1],
                              method = "knn",
                              trControl = fitControl,
                              tuneLength =10,
                              preProcess= c('zv', 'nzv', 'center', 'scale'),
                              metric='logLoss')


vel_knn_dat_train <- vel_train_data
vel_knn_dat_test <- vel_test_data
vel_indx <- sapply(vel_knn_dat_train[,-1], is.factor)
vel_knn_dat_train[vel_indx] <- lapply(vel_knn_dat_train[vel_indx], function(x) as.numeric(x))
vel_indx <- sapply(vel_knn_dat_test[,-1], is.factor)
vel_knn_dat_test[vel_indx] <- lapply(vel_knn_dat_test[vel_indx], function(x) as.numeric(x))
vel_knnFit <-  caret::train (x=vel_knn_dat_train [,-1],
                             y=vel_knn_dat_train[,1],
                             method = "knn",
                             trControl = fitControl,
                             tuneLength =10,
                             preProcess= c('zv', 'nzv', 'center', 'scale'),
                             metric='logLoss')

accel_knn_dat_train <- accel_train_data
accel_knn_dat_test <- accel_test_data
accel_indx <- sapply(accel_knn_dat_train[,-1], is.factor)
accel_knn_dat_train[accel_indx] <- lapply(accel_knn_dat_train[accel_indx], function(x) as.numeric(x))
accel_indx <- sapply(accel_knn_dat_test[,-1], is.factor)
accel_knn_dat_test[accel_indx] <- lapply(accel_knn_dat_test[accel_indx], function(x) as.numeric(x))
accel_knnFit <-  caret::train (x=accel_knn_dat_train [,-1],
                               y=accel_knn_dat_train[,1],
                               method = "knn",
                               trControl = fitControl,
                               tuneLength =10,
                               preProcess= c('zv', 'nzv', 'center', 'scale'),
                               metric='logLoss')

abs_accel_knn_dat_train <- abs_accel_train_data
abs_accel_knn_dat_test <- abs_accel_test_data
abs_accel_indx <- sapply(abs_accel_knn_dat_train[,-1], is.factor)
abs_accel_knn_dat_train[abs_accel_indx] <- lapply(abs_accel_knn_dat_train[abs_accel_indx], function(x) as.numeric(x))
abs_accel_indx <- sapply(abs_accel_knn_dat_test[,-1], is.factor)
abs_accel_knn_dat_test[abs_accel_indx] <- lapply(abs_accel_knn_dat_test[abs_accel_indx], function(x) as.numeric(x))
abs_accel_knnFit <-  caret::train (x=abs_accel_knn_dat_train [,-1],
                                   y=abs_accel_knn_dat_train[,1],
                                   method = "knn",
                                   trControl = fitControl,
                                   tuneLength =10,
                                   preProcess= c('zv', 'nzv', 'center', 'scale'),
                                   metric='logLoss')
#end <- proc.time()
#knn_time <- end - strt
##knn_time
#knnFit

base_knn_results <- confusionMatrix(predict(base_knnFit, base_knn_dat_test), base_knn_dat_test[,1])
base_knn_results
base_prob_results <- predict(base_knnFit,base_knn_dat_test, type = "prob")
base_knn_ROC <- multiclass.roc(base_knn_dat_test[,1], base_prob_results[,1],levels = levels(base_knn_dat_test[,1]))
base_knn_ROC

vel_knn_results <- confusionMatrix(predict(vel_knnFit, vel_knn_dat_test), vel_knn_dat_test[,1])
vel_knn_results
vel_prob_results <- predict(vel_knnFit, vel_knn_dat_test, type = "prob")
vel_knn_ROC <- multiclass.roc(vel_knn_dat_test[,1], vel_prob_results[,1],levels = levels(vel_knn_dat_test[,1]))
vel_knn_ROC

accel_knn_results <- confusionMatrix(predict(accel_knnFit, accel_knn_dat_test), accel_knn_dat_test[,1])
accel_knn_results
accel_prob_results <- predict(accel_knnFit, accel_knn_dat_test, type = "prob")
accel_knn_ROC <- multiclass.roc(accel_knn_dat_test[,1], accel_prob_results[,1],levels = levels(accel_knn_dat_test[,1]))
accel_knn_ROC

abs_accel_knn_results <- confusionMatrix(predict(abs_accel_knnFit, abs_accel_knn_dat_test), abs_accel_knn_dat_test[,1])
abs_accel_knn_results
abs_accel_prob_results <- predict(abs_accel_knnFit, abs_accel_knn_dat_test, type = "prob")
abs_accel_knn_ROC <- multiclass.roc(abs_accel_knn_dat_test[,1], abs_accel_prob_results[,1],levels = levels(abs_accel_knn_dat_test[,1]))
abs_accel_knn_ROC

#total_time_end <- proc.time()
#total_time_spent <- total_time_end - total_time_start


#base_resamps <- resamples(list(RF = base_rfFit,
#                          GBM = base_bl_gbmFit,
#                          SVM = base_svmrFit,
#                          KNN = base_knnFit))


#Results
#theme1 <- trellis.par.get(xl)
#theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
#theme1$plot.symbol$pch = 16
#theme1$plot.line$col = rgb(1, 0, 0, .7)
#theme1$plot.line$lwd <- 2
#trellis.par.set(theme1)
#bwplot(base_resamps, layout = c(0, 1))
#?bwplot#

####################################################################################################
####################### Testing model prediction ReliefF ###########################################
####################################################################################################

# Using the same splits as previous models for comparative purposes
base_train_data <- base_data[base_split,]
base_test_data<- base_data[-base_split,]
vel_train_data <- vel_data[vel_split,]
vel_test_data<- vel_data[-vel_split,]
accel_train_data <- accel_data[accel_split,]
accel_test_data<- accel_data[-accel_split,]
abs_accel_train_data <- abs_accel_data[abs_accel_split,]
abs_accel_test_data<- abs_accel_data[-abs_accel_split,]



# Feature selection via reliefF method with 1000 permutations and .95 statistical significance

base_perm <- permuteRelief(x = base_train_data[,names(base_train_data)!='Diagnosis'],
                           y = base_train_data$Diagnosis,
                           nperm= 1000,
                           estimator = 'ReliefFequalK')
base_pred.95<-names(sort(abs(base_perm$standardized[which(abs(base_perm$standardized)>=1.96)]), decreasing=T))
base_predictors<- base_pred.95


histogram(~ value|'mPACCtrails', 
          data = base_perm$permutations)


vel_perm <- permuteRelief(x = vel_train_data[,names(vel_train_data)!='Diagnosis'],
                          y = vel_train_data$Diagnosis,
                          nperm= 1000,
                          estimator = 'ReliefFequalK')
vel_pred.95<-names(sort(abs(vel_perm$standardized[which(abs(vel_perm$standardized)>=1.96)]), decreasing=T))
vel_predictors<- vel_pred.95


accel_perm <- permuteRelief(x = accel_train_data[,names(accel_train_data)!=c('Diagnosis')],
                            y = accel_train_data$Diagnosis,
                            nperm= 1000,
                            estimator = 'ReliefFequalK')
accel_pred.95<-names(sort(abs(accel_perm$standardized[which(abs(accel_perm$standardized)>=1.96)]), decreasing=T))
accel_predictors<- accel_pred.95


abs_accel_perm <- permuteRelief(x = abs_accel_train_data[,names(abs_accel_train_data)!=c('Diagnosis')],
                                y = abs_accel_train_data$Diagnosis,
                                nperm= 1000,
                                estimator = 'ReliefFequalK')
abs_accel_pred.95<-names(sort(abs(abs_accel_perm$standardized[which(abs(abs_accel_perm$standardized)>=1.96)]), decreasing=T))
abs_accel_predictors<- abs_accel_pred.95

#subsetting training data with only predictors from feature selection and diagnosis 
base_train_data <- base_train_data[c('Diagnosis',base_predictors)]
vel_train_data <- vel_train_data[c('Diagnosis',vel_predictors)]
accel_train_data <- accel_train_data[c('Diagnosis',accel_predictors)]
abs_accel_train_data <- abs_accel_train_data[c('Diagnosis',abs_accel_predictors)]


# Impute the missing values in the training data using random forest
base_train_data <- rfImpute(Diagnosis ~ ., base_train_data)
vel_train_data <- rfImpute(Diagnosis ~ ., vel_train_data)
accel_train_data <- rfImpute(Diagnosis ~ ., accel_train_data)
abs_accel_train_data <- rfImpute(Diagnosis ~ ., abs_accel_train_data)

# dummification of training data for nominal values
dummies <- dummyVars(Diagnosis~., data = base_train_data)
data_numeric <-predict(dummies, base_train_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(base_train_data$Diagnosis,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
base_train_data<-data_numeric 

dummies <- dummyVars(Diagnosis~., data = vel_train_data)
data_numeric <-predict(dummies, vel_train_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(vel_train_data$Diagnosis,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
vel_train_data<-data_numeric 

dummies <- dummyVars(Diagnosis~., data = accel_train_data)
data_numeric <-predict(dummies, accel_train_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(accel_train_data$Diagnosis,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
accel_train_data<-data_numeric 

dummies <- dummyVars(Diagnosis~., data = abs_accel_train_data)
data_numeric <-predict(dummies, abs_accel_train_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(abs_accel_train_data$Diagnosis,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
abs_accel_train_data<-data_numeric 


#Shuffle train_data after smote
base_train_data<-base_train_data[sample(nrow(base_train_data)),]
vel_train_data<-vel_train_data[sample(nrow(vel_train_data)),]
accel_train_data<-accel_train_data[sample(nrow(accel_train_data)),]
abs_accel_train_data<-abs_accel_train_data[sample(nrow(abs_accel_train_data)),]


#subsetting test data with only predictors from feature selection 
base_test_data <- base_test_data[c('Diagnosis',base_predictors)]
vel_test_data <- vel_test_data[c('Diagnosis',vel_predictors)]
accel_test_data <- accel_test_data[c('Diagnosis',accel_predictors)]
abs_accel_test_data <- abs_accel_test_data[c('Diagnosis',abs_accel_predictors)]

#Imputing missing values in test data
base_test_data <- rfImpute(Diagnosis ~ ., base_test_data)
vel_test_data <- rfImpute(Diagnosis ~ ., vel_test_data)
accel_test_data <- rfImpute(Diagnosis ~ ., accel_test_data)
abs_accel_test_data <- rfImpute(Diagnosis ~ ., abs_accel_test_data)

#dummification of test data
dummies <- dummyVars(Diagnosis~., data = base_test_data)
data_numeric <-predict(dummies, base_test_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(base_test_data$Diagnosis,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
base_test_data<-data_numeric 

dummies <- dummyVars(Diagnosis~., data = vel_test_data)
data_numeric <-predict(dummies, vel_test_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(vel_test_data$Diagnosis,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
vel_test_data<-data_numeric 

dummies <- dummyVars(Diagnosis~., data = accel_test_data)
data_numeric <-predict(dummies, accel_test_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(accel_test_data$Diagnosis,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
accel_test_data<-data_numeric 

dummies <- dummyVars(Diagnosis~., data = abs_accel_test_data)
data_numeric <-predict(dummies, abs_accel_test_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(abs_accel_test_data$Diagnosis,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
abs_accel_test_data<-data_numeric 



#######################################################################################################
######################### Testing 4 models (gbm, rf, svm, knn) with RelieF ############################
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
#strt <- proc.time()

bl_gbmGrid <- expand.grid(nrounds = c(1, 10,100),
                          max_depth = c(1, 4,6),
                          eta = c(.1, .4),
                          gamma = 0,
                          colsample_bytree = .7,
                          min_child_weight = 1,
                          subsample = c(.5,.8,1))

perm_base_bl_gbmFit <-  caret::train (x=base_train_data [,-1],
                                      y=base_train_data[,1],
                                      method = "xgbTree",
                                      tuneGrid = bl_gbmGrid,
                                      trControl = fitControl,
                                      preProcess= c('zv','nzv'),
                                      tuneLength =10,
                                      verbose = FALSE,
                                      metric='logLoss')

perm_vel_bl_gbmFit <-  caret::train (x=vel_train_data [,-1],
                                     y=vel_train_data[,1],
                                     method = "xgbTree",
                                     tuneGrid = bl_gbmGrid,
                                     trControl = fitControl,
                                     preProcess= c('zv','nzv'),
                                     tuneLength =10,
                                     verbose = FALSE,
                                     metric='logLoss')

perm_accel_bl_gbmFit <-  caret::train (x=accel_train_data [,-1],
                                       y=accel_train_data[,1],
                                       method = "xgbTree",
                                       tuneGrid = bl_gbmGrid,
                                       trControl = fitControl,
                                       preProcess= c('zv','nzv'),
                                       tuneLength =10,
                                       verbose = FALSE,
                                       metric='logLoss')

perm_abs_accel_bl_gbmFit <-  caret::train (x=abs_accel_train_data [,-1],
                                           y=abs_accel_train_data[,1],
                                           method = "xgbTree",
                                           tuneGrid = bl_gbmGrid,
                                           trControl = fitControl,
                                           preProcess= c('zv','nzv'),
                                           tuneLength =10,
                                           verbose = FALSE,
                                           metric='logLoss')


#Testing
perm_base_bl_gbm_results <- confusionMatrix(predict(perm_base_bl_gbmFit, base_test_data), base_test_data[,1])
perm_base_bl_gbm_results
perm_base_prob_results <- predict(perm_base_bl_gbmFit, base_test_data, type = "prob")
perm_base_bl_gbm_ROC <- multiclass.roc(base_test_data[,1], perm_base_prob_results[,1],levels = levels(base_test_data[,1]))
perm_base_bl_gbm_ROC

perm_vel_bl_gbm_results <- confusionMatrix(predict(perm_vel_bl_gbmFit, vel_test_data), vel_test_data[,1])
perm_vel_bl_gbm_results
perm_vel_prob_results <- predict(perm_vel_bl_gbmFit, vel_test_data, type = "prob")
perm_vel_bl_gbm_ROC <- multiclass.roc(vel_test_data[,1], perm_vel_prob_results[,1],levels = levels(vel_test_data[,1]))
perm_vel_bl_gbm_ROC

perm_accel_bl_gbm_results <- confusionMatrix(predict(perm_accel_bl_gbmFit, accel_test_data), accel_test_data[,1])
perm_accel_bl_gbm_results
perm_accel_prob_results <- predict(perm_accel_bl_gbmFit, accel_test_data, type = "prob")
perm_accel_bl_gbm_ROC <- multiclass.roc(accel_test_data[,1], perm_accel_prob_results[,1],levels = levels(accel_test_data[,1]))
perm_accel_bl_gbm_ROC

perm_abs_accel_bl_gbm_results <- confusionMatrix(predict(perm_abs_accel_bl_gbmFit, abs_accel_test_data), abs_accel_test_data[,1])
perm_abs_accel_bl_gbm_results
perm_abs_accel_prob_results <- predict(perm_abs_accel_bl_gbmFit, abs_accel_test_data, type = "prob")
perm_abs_accel_bl_gbm_ROC <- multiclass.roc(abs_accel_test_data[,1], perm_abs_accel_prob_results[,1],levels = levels(abs_accel_test_data[,1]))
perm_abs_accel_bl_gbm_ROC


#===================================================Random Forest==================================================
#strt <- proc.time()
perm_base_rfFit <-  caret::train (x=base_train_data[,-1],
                                  y=base_train_data[,1],
                                  method = "rf",
                                  trControl = fitControl,
                                  ntrees = 500,
                                  preProcess= c('zv', 'nzv', 'center', 'scale'),
                                  tuneLength =10,
                                  metric='logLoss')

perm_vel_rfFit <-  caret::train (x=vel_train_data[,-1],
                                 y=vel_train_data[,1],
                                 method = "rf",
                                 trControl = fitControl,
                                 ntrees = 500,
                                 preProcess= c('zv', 'nzv', 'center', 'scale'),
                                 tuneLength =10,
                                 metric='logLoss')

perm_accel_rfFit <-  caret::train (x=accel_train_data[,-1],
                                   y=accel_train_data[,1],
                                   method = "rf",
                                   trControl = fitControl,
                                   ntrees = 500,
                                   preProcess= c('zv', 'nzv', 'center', 'scale'),
                                   tuneLength =10,
                                   metric='logLoss')

perm_abs_accel_rfFit <-  caret::train (x=abs_accel_train_data[,-1],
                                       y=abs_accel_train_data[,1],
                                       method = "rf",
                                       trControl = fitControl,
                                       ntrees = 500,
                                       preProcess= c('zv', 'nzv', 'center', 'scale'),
                                       tuneLength =10,
                                       metric='logLoss')

perm_base_rf_results <- confusionMatrix(predict(perm_base_rfFit, base_test_data), base_test_data[,1])
perm_base_rf_results
perm_base_prob_results <- predict(perm_base_rfFit, base_test_data, type = 'prob')
perm_base_rf_ROC <- multiclass.roc(base_test_data[,1], perm_base_prob_results[,1],levels = levels(base_test_data[,1]))
perm_base_rf_ROC

perm_vel_rf_results <- confusionMatrix(predict(perm_vel_rfFit, vel_test_data), vel_test_data[,1])
perm_vel_rf_results
perm_vel_prob_results <- predict(perm_vel_rfFit, vel_test_data, type = 'prob')
perm_vel_rf_ROC <- multiclass.roc(vel_test_data[,1], perm_vel_prob_results[,1],levels = levels(vel_test_data[,1]))
perm_vel_rf_ROC

perm_accel_rf_results <- confusionMatrix(predict(perm_accel_rfFit, accel_test_data), accel_test_data[,1])
perm_accel_rf_results
perm_accel_prob_results <- predict(perm_accel_rfFit, accel_test_data, type = 'prob')
perm_accel_rf_ROC <- multiclass.roc(accel_test_data[,1], perm_accel_prob_results[,1],levels = levels(accel_test_data[,1]))
perm_accel_rf_ROC

perm_abs_accel_rf_results <- confusionMatrix(predict(perm_abs_accel_rfFit, abs_accel_test_data), abs_accel_test_data[,1])
perm_abs_accel_rf_results
perm_abs_accel_prob_results <- predict(perm_abs_accel_rfFit, abs_accel_test_data, type = 'prob')
perm_abs_accel_rf_ROC <- multiclass.roc(abs_accel_test_data[,1], perm_abs_accel_prob_results[,1],levels = levels(abs_accel_test_data[,1]))
perm_abs_accel_rf_ROC


#==============================Support vector machine radial kernel=========================
#strt <- proc.time()
perm_base_svmrFit <-  caret::train(Diagnosis~.,
                              data = base_train_data,
                              metric = 'logLoss',
                              method = "svmRadial",
                              tuneLength = 12,
                              preProcess= c('zv','nzv', 'center', 'scale'),
                              trControl = fitControl)

perm_vel_svmrFit <-  caret::train(Diagnosis~.,
                                  data = vel_train_data,
                                  metric = 'logLoss',
                                  method = "svmRadial",
                                  tuneLength = 12,
                                  preProcess= c('zv','nzv', 'center', 'scale'),
                                  trControl = fitControl)

perm_accel_svmrFit <-  caret::train(Diagnosis~.,
                                    data = accel_train_data,
                                    metric = 'logLoss',
                                    method = "svmRadial",
                                    tuneLength = 12,
                                    preProcess= c('zv','nzv', 'center', 'scale'),
                                    trControl = fitControl)

perm_abs_accel_svmrFit <-  caret::train(Diagnosis~.,
                                        data = abs_accel_train_data,
                                        metric = 'logLoss',
                                        method = "svmRadial",
                                        tuneLength = 12,
                                        preProcess= c('zv','nzv', 'center', 'scale'),
                                        trControl = fitControl)

perm_base_svmresults <- confusionMatrix(predict(perm_base_svmrFit, base_test_data), base_test_data[,1])
perm_base_svmresults
perm_base_prob_results <- predict(perm_base_svmrFit, base_test_data, type = "prob")
perm_base_svmrROC <- multiclass.roc(base_test_data[,1],perm_base_prob_results[,1],levels = levels(base_test_data[,1]))
perm_base_svmrROC

perm_vel_svmresults <- confusionMatrix(predict(perm_vel_svmrFit, vel_test_data), vel_test_data[,1])
perm_vel_svmresults
perm_vel_prob_results <- predict(perm_vel_svmrFit, vel_test_data, type = "prob")
perm_vel_svmrROC <- multiclass.roc(vel_test_data[,1], perm_vel_prob_results[,1],levels = levels(vel_test_data[,1]))
perm_vel_svmrROC

perm_accel_svmresults <- confusionMatrix(predict(perm_accel_svmrFit, accel_test_data), accel_test_data[,1])
perm_accel_svmresults
perm_accel_prob_results <- predict(perm_accel_svmrFit, accel_test_data, type = "prob")
perm_accel_svmrROC <- multiclass.roc(accel_test_data[,1],perm_accel_prob_results[,1],levels = levels(accel_test_data[,1]))
perm_accel_svmrROC

perm_abs_accel_svmresults <- confusionMatrix(predict(perm_abs_accel_svmrFit, abs_accel_test_data), abs_accel_test_data[,1])
perm_abs_accel_svmresults
perm_abs_accel_prob_results <- predict(perm_abs_accel_svmrFit, abs_accel_test_data, type = "prob")
perm_abs_accel_svmrROC <- multiclass.roc(abs_accel_test_data[,1], perm_abs_accel_prob_results[,1],levels = levels(abs_accel_test_data[,1]))
perm_abs_accel_svmrROC

#================================== K-Nearest Neighbor ===========================================
#strt <- proc.time()

base_knn_dat_train <- base_train_data
base_knn_dat_test <- base_test_data
base_indx <- sapply(base_knn_dat_train[,-1], is.factor)
base_knn_dat_train[base_indx] <- lapply(base_knn_dat_train[base_indx], function(x) as.numeric(x))
base_indx <- sapply(base_knn_dat_test[,-1], is.factor)
base_knn_dat_test[base_indx] <- lapply(base_knn_dat_test[base_indx], function(x) as.numeric(x))
perm_base_knnFit <-  caret::train (x=base_knn_dat_train [,-1],
                                   y=base_knn_dat_train[,1],
                                   method = "knn",
                                   trControl = fitControl,
                                   tuneLength =10,
                                   preProcess= c('zv', 'nzv', 'center', 'scale'),
                                   metric='logLoss')


vel_knn_dat_train <- vel_train_data
vel_knn_dat_test <- vel_test_data
vel_indx <- sapply(vel_knn_dat_train[,-1], is.factor)
vel_knn_dat_train[vel_indx] <- lapply(vel_knn_dat_train[vel_indx], function(x) as.numeric(x))
vel_indx <- sapply(vel_knn_dat_test[,-1], is.factor)
vel_knn_dat_test[vel_indx] <- lapply(vel_knn_dat_test[vel_indx], function(x) as.numeric(x))
perm_vel_knnFit <-  caret::train (x=vel_knn_dat_train [,-1],
                                  y=vel_knn_dat_train[,1],
                                  method = "knn",
                                  trControl = fitControl,
                                  tuneLength =10,
                                  preProcess= c('zv', 'nzv', 'center', 'scale'),
                                  metric='logLoss')

accel_knn_dat_train <- accel_train_data
accel_knn_dat_test <- accel_test_data
accel_indx <- sapply(accel_knn_dat_train[,-1], is.factor)
accel_knn_dat_train[accel_indx] <- lapply(accel_knn_dat_train[accel_indx], function(x) as.numeric(x))
accel_indx <- sapply(accel_knn_dat_test[,-1], is.factor)
accel_knn_dat_test[accel_indx] <- lapply(accel_knn_dat_test[accel_indx], function(x) as.numeric(x))
perm_accel_knnFit <-  caret::train (x=accel_knn_dat_train [,-1],
                                    y=accel_knn_dat_train[,1],
                                    method = "knn",
                                    trControl = fitControl,
                                    tuneLength =10,
                                    preProcess= c('zv', 'nzv', 'center', 'scale'),
                                    metric='logLoss')

abs_accel_knn_dat_train <- abs_accel_train_data
abs_accel_knn_dat_test <- abs_accel_test_data
abs_accel_indx <- sapply(abs_accel_knn_dat_train[,-1], is.factor)
abs_accel_knn_dat_train[abs_accel_indx] <- lapply(abs_accel_knn_dat_train[abs_accel_indx], function(x) as.numeric(x))
abs_accel_indx <- sapply(abs_accel_knn_dat_test[,-1], is.factor)
abs_accel_knn_dat_test[abs_accel_indx] <- lapply(abs_accel_knn_dat_test[abs_accel_indx], function(x) as.numeric(x))
perm_abs_accel_knnFit <-  caret::train (x=abs_accel_knn_dat_train [,-1],
                                        y=abs_accel_knn_dat_train[,1],
                                        method = "knn",
                                        trControl = fitControl,
                                        tuneLength =10,
                                        preProcess= c('zv', 'nzv', 'center', 'scale'),
                                        metric='logLoss')

perm_base_knn_results <- confusionMatrix(predict(perm_base_knnFit, base_knn_dat_test), base_knn_dat_test[,1])
perm_base_knn_results
perm_base_prob_results <- predict(perm_base_knnFit,base_knn_dat_test, type = "prob")
perm_base_knn_ROC <- multiclass.roc(base_knn_dat_test[,1], perm_base_prob_results[,1],levels = levels(base_knn_dat_test[,1]))
perm_base_knn_ROC

perm_vel_knn_results <- confusionMatrix(predict(perm_vel_knnFit, vel_knn_dat_test), vel_knn_dat_test[,1])
perm_vel_knn_results
perm_vel_prob_results <- predict(perm_vel_knnFit, vel_knn_dat_test, type = "prob")
perm_vel_knn_ROC <- multiclass.roc(vel_knn_dat_test[,1], perm_vel_prob_results[,1],levels = levels(vel_knn_dat_test[,1]))
perm_vel_knn_ROC

perm_accel_knn_results <- confusionMatrix(predict(perm_accel_knnFit, accel_knn_dat_test), accel_knn_dat_test[,1])
perm_accel_knn_results
perm_accel_prob_results <- predict(perm_accel_knnFit, accel_knn_dat_test, type = "prob")
perm_accel_knn_ROC <- multiclass.roc(accel_knn_dat_test[,1], perm_accel_prob_results[,1],levels = levels(accel_knn_dat_test[,1]))
perm_accel_knn_ROC

perm_abs_accel_knn_results <- confusionMatrix(predict(perm_abs_accel_knnFit, abs_accel_knn_dat_test), abs_accel_knn_dat_test[,1])
perm_abs_accel_knn_results
perm_abs_accel_prob_results <- predict(perm_abs_accel_knnFit, abs_accel_knn_dat_test, type = "prob")
perm_abs_accel_knn_ROC <- multiclass.roc(abs_accel_knn_dat_test[,1], perm_abs_accel_prob_results[,1],levels = levels(abs_accel_knn_dat_test[,1]))
perm_abs_accel_knn_ROC

#================Recording the Results =======================

IG_results <- NULL
perm_results <- NULL

#======================== Feature selection with Information gain results ===================
# Base Dataset
base_gbm_accuracy <- base_bl_gbm_results$overall[1]
base_gbm_kappa <- base_bl_gbm_results$overall[2]
base_gbm_AUC <- base_bl_gbm_ROC$auc[1]
base_gbm_avg_sens <- (base_bl_gbm_results[[4]][1]+base_bl_gbm_results[[4]][2]+
                      base_bl_gbm_results[[4]][3]+base_bl_gbm_results[[4]][4]+
                      base_bl_gbm_results[[4]][5]+base_bl_gbm_results[[4]][6]+base_bl_gbm_results[[4]][7])/7

base_gbm_avg_spec <- (base_bl_gbm_results[[4]][8]+base_bl_gbm_results[[4]][9]+
                      base_bl_gbm_results[[4]][10]+base_bl_gbm_results[[4]][11]+
                      base_bl_gbm_results[[4]][12]+base_bl_gbm_results[[4]][13]+base_bl_gbm_results[[4]][14])/7



base_rf_accuracy <- base_rf_results$overall[1]
base_rf_kappa <- base_rf_results$overall[2]
base_rf_AUC <- base_rf_ROC$auc[1]
base_rf_avg_sens <- (base_rf_results[[4]][1]+base_rf_results[[4]][2]+
                       base_rf_results[[4]][3]+base_rf_results[[4]][4]+
                       base_rf_results[[4]][5]+base_rf_results[[4]][6]+base_rf_results[[4]][7])/7

base_rf_avg_spec <- (base_rf_results[[4]][8]+base_rf_results[[4]][9]+
                       base_rf_results[[4]][10]+base_rf_results[[4]][11]+
                       base_rf_results[[4]][12]+base_rf_results[[4]][13]+base_rf_results[[4]][14])/7

base_svm_accuracy <- base_svmresults$overall[1]
base_svm_kappa <- base_svmresults$overall[2]
base_svm_AUC <- base_svmrROC$auc[1]
base_svm_avg_sens <- (base_svmresults[[4]][1]+base_svmresults[[4]][2]+
                        base_svmresults[[4]][3]+base_svmresults[[4]][4]+
                        base_svmresults[[4]][5]+base_svmresults[[4]][6]+base_svmresults[[4]][7])/7

base_svm_avg_spec <- (base_svmresults[[4]][8]+base_svmresults[[4]][9]+
                        base_svmresults[[4]][10]+base_svmresults[[4]][11]+
                        base_svmresults[[4]][12]+base_svmresults[[4]][13]+base_svmresults[[4]][14])/7

base_knn_accuracy <- base_knn_results$overall[1]
base_knn_kappa <- base_knn_results$overall[2]
base_knn_AUC <- base_knn_ROC$auc[1]
base_knn_avg_sens <- (base_knn_results[[4]][1]+base_knn_results[[4]][2]+
                        base_knn_results[[4]][3]+base_knn_results[[4]][4]+
                        base_knn_results[[4]][5]+base_knn_results[[4]][6]+base_knn_results[[4]][7])/7

base_knn_avg_spec <- (base_knn_results[[4]][8]+base_knn_results[[4]][9]+
                        base_knn_results[[4]][10]+base_knn_results[[4]][11]+
                        base_knn_results[[4]][12]+base_knn_results[[4]][13]+base_knn_results[[4]][14])/7



#Velocity Dataset
vel_gbm_accuracy <- vel_bl_gbm_results$overall[1]
vel_gbm_kappa <- vel_bl_gbm_results$overall[2]
vel_gbm_AUC <- vel_bl_gbm_ROC$auc[1]
vel_gbm_avg_sens <- (vel_bl_gbm_results[[4]][1]+vel_bl_gbm_results[[4]][2]+
                       vel_bl_gbm_results[[4]][3]+vel_bl_gbm_results[[4]][4]+
                       vel_bl_gbm_results[[4]][5]+vel_bl_gbm_results[[4]][6]+vel_bl_gbm_results[[4]][7])/7

vel_gbm_avg_spec <- (vel_bl_gbm_results[[4]][8]+vel_bl_gbm_results[[4]][9]+
                       vel_bl_gbm_results[[4]][10]+vel_bl_gbm_results[[4]][11]+
                       vel_bl_gbm_results[[4]][12]+vel_bl_gbm_results[[4]][13]+vel_bl_gbm_results[[4]][14])/7


vel_rf_accuracy <- vel_rf_results$overall[1]
vel_rf_kappa <- vel_rf_results$overall[2]
vel_rf_AUC <- vel_rf_ROC$auc[1]
vel_rf_avg_sens <- (vel_rf_results[[4]][1]+vel_rf_results[[4]][2]+
                      vel_rf_results[[4]][3]+vel_rf_results[[4]][4]+
                      vel_rf_results[[4]][5]+vel_rf_results[[4]][6]+vel_rf_results[[4]][7])/7

vel_rf_avg_spec <- (vel_rf_results[[4]][8]+vel_rf_results[[4]][9]+
                      vel_rf_results[[4]][10]+vel_rf_results[[4]][11]+
                      vel_rf_results[[4]][12]+vel_rf_results[[4]][13]+vel_rf_results[[4]][14])/7


vel_svm_accuracy <- vel_svmresults$overall[1]
vel_svm_kappa <- vel_svmresults$overall[2]
vel_svm_AUC <- vel_svmrROC$auc[1]
vel_svm_avg_sens <- (vel_svmresults[[4]][1]+vel_svmresults[[4]][2]+
                       vel_svmresults[[4]][3]+vel_svmresults[[4]][4]+
                       vel_svmresults[[4]][5]+vel_svmresults[[4]][6]+vel_svmresults[[4]][7])/7

vel_svm_avg_spec <- (vel_svmresults[[4]][8]+vel_svmresults[[4]][9]+
                       vel_svmresults[[4]][10]+vel_svmresults[[4]][11]+
                       vel_svmresults[[4]][12]+vel_svmresults[[4]][13]+vel_svmresults[[4]][14])/7


vel_knn_accuracy <- vel_knn_results$overall[1]
vel_knn_kappa <- vel_knn_results$overall[2]
vel_knn_AUC <- vel_knn_ROC$auc[1]
vel_knn_avg_sens <- (vel_knn_results[[4]][1]+vel_knn_results[[4]][2]+
                       vel_knn_results[[4]][3]+vel_knn_results[[4]][4]+
                       vel_knn_results[[4]][5]+vel_knn_results[[4]][6]+vel_knn_results[[4]][7])/7

vel_knn_avg_spec <- (vel_knn_results[[4]][8]+vel_knn_results[[4]][9]+
                       vel_knn_results[[4]][10]+vel_knn_results[[4]][11]+
                       vel_knn_results[[4]][12]+vel_knn_results[[4]][13]+vel_knn_results[[4]][14])/7


#Acceleration Dataset
accel_gbm_accuracy <- accel_bl_gbm_results$overall[1]
accel_gbm_kappa <- accel_bl_gbm_results$overall[2]
accel_gbm_AUC <- accel_bl_gbm_ROC$auc[1]
accel_gbm_avg_sens <- mean(accel_bl_gbm_results[[4]][1:7])
accel_gbm_avg_spec <- mean(accel_bl_gbm_results[[4]][8:14])



accel_rf_accuracy <- accel_rf_results$overall[1]
accel_rf_kappa <- accel_rf_results$overall[2]
accel_rf_AUC <- accel_rf_ROC$auc[1]
accel_rf_avg_sens <- mean(accel_rf_results[[4]][1:7])
accel_rf_avg_spec <- mean(accel_rf_results[[4]][8:14])

accel_svm_accuracy <- accel_svmresults$overall[1]
accel_svm_kappa <- accel_svmresults$overall[2]
accel_svm_AUC <- accel_svmrROC$auc[1]
accel_svm_avg_sens <- mean(accel_svmresults[[4]][1:7])
accel_svm_avg_spec <- mean(accel_svmresults[[4]][8:14])

accel_knn_accuracy <- accel_knn_results$overall[1]
accel_knn_kappa <- accel_knn_results$overall[2]
accel_knn_AUC <- accel_knn_ROC$auc[1]
accel_knn_avg_sens <- mean(accel_knn_results[[4]][1:7])
accel_knn_avg_spec <- mean(accel_knn_results[[4]][8:14])



#Absolute Acceleration Dataset
abs_accel_gbm_accuracy <- abs_accel_bl_gbm_results$overall[1]
abs_accel_gbm_kappa <- abs_accel_bl_gbm_results$overall[2]
abs_accel_gbm_AUC <- abs_accel_bl_gbm_ROC$auc[1]
abs_accel_gbm_avg_sens <- mean(abs_accel_bl_gbm_results[[4]][1:7])
abs_accel_gbm_avg_spec <- mean(abs_accel_bl_gbm_results[[4]][8:14])


abs_accel_rf_accuracy <- abs_accel_rf_results$overall[1]
abs_accel_rf_kappa <- abs_accel_rf_results$overall[2]
abs_accel_rf_AUC <- abs_accel_rf_ROC$auc[1]
abs_accel_rf_avg_sens <- mean(abs_accel_rf_results[[4]][1:7])
abs_accel_rf_avg_spec <- mean(abs_accel_rf_results[[4]][8:14])

abs_accel_svm_accuracy <- abs_accel_svmresults$overall[1]
abs_accel_svm_kappa <- abs_accel_svmresults$overall[2]
abs_accel_svm_AUC <- abs_accel_svmrROC$auc[1]
abs_accel_svm_avg_sens <- mean(abs_accel_svmresults[[4]][1:7])
abs_accel_svm_avg_spec <- mean(abs_accel_svmresults[[4]][8:14])

abs_accel_knn_accuracy <- abs_accel_knn_results$overall[1]
abs_accel_knn_kappa <- abs_accel_knn_results$overall[2]
abs_accel_knn_AUC <- abs_accel_knn_ROC$auc[1]
abs_accel_knn_avg_sens <- mean(abs_accel_knn_results[[4]][1:7])
abs_accel_knn_avg_spec <- mean(abs_accel_knn_results[[4]][8:14])

IG_results$accuracy <- rbind(base_gbm_accuracy,base_rf_accuracy,base_svm_accuracy,base_knn_accuracy,
                             vel_gbm_accuracy,vel_rf_accuracy,vel_svm_accuracy,vel_knn_accuracy,
                             accel_gbm_accuracy, accel_rf_accuracy,accel_svm_accuracy,accel_knn_accuracy,
                             abs_accel_gbm_accuracy,abs_accel_rf_accuracy,abs_accel_svm_accuracy,abs_accel_knn_accuracy)

IG_results$avg_sens <- rbind(base_gbm_avg_sens,base_rf_avg_sens,base_svm_avg_sens,base_knn_avg_sens,
                             vel_gbm_avg_sens, vel_rf_avg_sens, vel_svm_avg_sens, vel_knn_avg_sens,
                             accel_gbm_avg_sens, accel_rf_avg_sens, accel_svm_avg_sens, accel_knn_avg_sens,
                             abs_accel_gbm_avg_sens, abs_accel_rf_avg_sens, abs_accel_svm_avg_sens, abs_accel_knn_avg_sens)

IG_results$avg_spec <- rbind(base_gbm_avg_spec,base_rf_avg_spec,base_svm_avg_spec,base_knn_avg_spec,
                             vel_gbm_avg_spec, vel_rf_avg_spec, vel_svm_avg_spec, vel_knn_avg_spec,
                             accel_gbm_avg_spec, accel_rf_avg_spec, accel_svm_avg_spec, accel_knn_avg_spec,
                             abs_accel_gbm_avg_spec, abs_accel_rf_avg_spec, abs_accel_svm_avg_spec, abs_accel_knn_avg_spec)

IG_results$kappa <- rbind(base_gbm_kappa,base_rf_kappa,base_svm_kappa,base_knn_kappa,
                          vel_gbm_kappa,vel_rf_kappa,vel_svm_kappa,vel_knn_kappa,
                          accel_gbm_kappa,accel_rf_kappa,accel_svm_kappa,accel_knn_kappa,
                          abs_accel_gbm_kappa,abs_accel_rf_kappa,abs_accel_svm_kappa,abs_accel_knn_kappa)

IG_results$AUC <- rbind(base_gbm_AUC,base_rf_AUC,base_svm_AUC,base_knn_AUC,
                        vel_gbm_AUC,vel_rf_AUC,vel_svm_AUC,vel_knn_AUC,
                        accel_gbm_AUC,accel_rf_AUC, accel_svm_AUC,accel_knn_AUC,
                        abs_accel_gbm_AUC,abs_accel_rf_AUC,abs_accel_svm_AUC,abs_accel_knn_AUC)

IG_results <- data.frame(IG_results)
rownames(IG_results)<- c('base_gbm','base_rf','base_svm','base_knn',
                         'vel_gbm','vel_rf','vel_svm','vel_knn',
                         'accel_gbm','accel_rf','accel_svm','accel_knn',
                         'abs_accel_gbm','abs_accel_rf','abs_accel_svm','abs_accel_knn')


#======================== Feature selection with relieF results ===================
perm_base_gbm_accuracy <- perm_base_bl_gbm_results$overall[1]
perm_base_gbm_kappa <- perm_base_bl_gbm_results$overall[2]
perm_base_gbm_AUC <- perm_base_bl_gbm_ROC$auc[1]
perm_base_gbm_avg_sens <- mean(perm_base_bl_gbm_results[[4]][1:7])
perm_base_gbm_avg_spec <- mean(perm_base_bl_gbm_results[[4]][8:14])

perm_base_rf_accuracy <- perm_base_rf_results$overall[1]
perm_base_rf_kappa <- perm_base_rf_results$overall[2]
perm_base_rf_AUC <- perm_base_rf_ROC$auc[1]
perm_base_rf_avg_sens <- mean(perm_base_rf_results[[4]][1:7])
perm_base_rf_avg_spec <- mean(perm_base_rf_results[[4]][8:14])

perm_base_svm_accuracy <- perm_base_svmresults$overall[1]
perm_base_svm_kappa <- perm_base_svmresults$overall[2]
perm_base_svm_AUC <- perm_base_svmrROC$auc[1]
perm_base_svm_avg_sens <- mean(perm_base_svmresults[[4]][1:7])
perm_base_svm_avg_spec <- mean(perm_base_svmresults[[4]][8:14])

perm_base_knn_accuracy <- perm_base_knn_results$overall[1]
perm_base_knn_kappa <- perm_base_knn_results$overall[2]
perm_base_knn_AUC <- perm_base_knn_ROC$auc[1]
perm_base_knn_avg_sens <- mean(perm_base_knn_results[[4]][1:7])
perm_base_knn_avg_spec <- mean(perm_base_knn_results[[4]][8:14])



perm_vel_gbm_accuracy <- perm_vel_bl_gbm_results$overall[1]
perm_vel_gbm_kappa <- perm_vel_bl_gbm_results$overall[2]
perm_vel_gbm_AUC <- perm_vel_bl_gbm_ROC$auc[1]
perm_vel_gbm_avg_sens <- mean(perm_vel_bl_gbm_results[[4]][1:7])
perm_vel_gbm_avg_spec <- mean(perm_vel_bl_gbm_results[[4]][8:14])

perm_vel_rf_accuracy <- perm_vel_rf_results$overall[1]
perm_vel_rf_kappa <- perm_vel_rf_results$overall[2]
perm_vel_rf_AUC <- perm_vel_rf_ROC$auc[1]
perm_vel_rf_avg_sens <- mean(perm_vel_rf_results[[4]][1:7])
perm_vel_rf_avg_spec <- mean(perm_vel_rf_results[[4]][8:14])



perm_vel_svm_accuracy <- perm_vel_svmresults$overall[1]
perm_vel_svm_kappa <- perm_vel_svmresults$overall[2]
perm_vel_svm_AUC <- perm_vel_svmrROC$auc[1]
perm_vel_svm_avg_sens <- mean(perm_vel_svmresults[[4]][1:7])
perm_vel_svm_avg_spec <- mean(perm_vel_svmresults[[4]][8:14])

perm_vel_knn_accuracy <- perm_vel_knn_results$overall[1]
perm_vel_knn_kappa <- perm_vel_knn_results$overall[2]
perm_vel_knn_AUC <- perm_vel_knn_ROC$auc[1]
perm_vel_knn_avg_sens <- mean(perm_vel_knn_results[[4]][1:7])
perm_vel_knn_avg_spec <- mean(perm_vel_knn_results[[4]][8:14])


perm_accel_gbm_accuracy <- perm_accel_bl_gbm_results$overall[1]
perm_accel_gbm_kappa <- perm_accel_bl_gbm_results$overall[2]
perm_accel_gbm_AUC <- perm_accel_bl_gbm_ROC$auc[1]
perm_accel_gbm_avg_sens <- mean(perm_accel_bl_gbm_results[[4]][1:7])
perm_accel_gbm_avg_spec <- mean(perm_accel_bl_gbm_results[[4]][8:14])

perm_accel_rf_accuracy <- perm_accel_rf_results$overall[1]
perm_accel_rf_kappa <- perm_accel_rf_results$overall[2]
perm_accel_rf_AUC <- perm_accel_rf_ROC$auc[1]
perm_accel_rf_avg_sens <- mean(perm_accel_rf_results[[4]][1:7])
perm_accel_rf_avg_spec <- mean(perm_accel_rf_results[[4]][8:14])

perm_accel_svm_accuracy <- perm_accel_svmresults$overall[1]
perm_accel_svm_kappa <- perm_accel_svmresults$overall[2]
perm_accel_svm_AUC <- perm_accel_svmrROC$auc[1]
perm_accel_svm_avg_sens <- mean(perm_accel_svmresults[[4]][1:7])
perm_accel_svm_avg_spec <- mean(perm_accel_svmresults[[4]][8:14])

perm_accel_knn_accuracy <- perm_accel_knn_results$overall[1]
perm_accel_knn_kappa <- perm_accel_knn_results$overall[2]
perm_accel_knn_AUC <- perm_accel_knn_ROC$auc[1]
perm_accel_knn_avg_sens <- mean(perm_accel_knn_results[[4]][1:7])
perm_accel_knn_avg_spec <- mean(perm_accel_knn_results[[4]][8:14])

perm_abs_accel_gbm_accuracy <- perm_abs_accel_bl_gbm_results$overall[1]
perm_abs_accel_gbm_kappa <- perm_abs_accel_bl_gbm_results$overall[2]
perm_abs_accel_gbm_AUC <- perm_abs_accel_bl_gbm_ROC$auc[1]
perm_abs_accel_gbm_avg_sens <- mean(perm_abs_accel_bl_gbm_results[[4]][1:7])
perm_abs_accel_gbm_avg_spec <- mean(perm_abs_accel_bl_gbm_results[[4]][8:14])

perm_abs_accel_rf_accuracy <- perm_abs_accel_rf_results$overall[1]
perm_abs_accel_rf_kappa <- perm_abs_accel_rf_results$overall[2]
perm_abs_accel_rf_AUC <- perm_abs_accel_rf_ROC$auc[1]
perm_abs_accel_rf_avg_sens <- mean(perm_abs_accel_rf_results[[4]][1:7])
perm_abs_accel_rf_avg_spec <- mean(perm_abs_accel_rf_results[[4]][8:14])

perm_abs_accel_svm_accuracy <- perm_abs_accel_svmresults$overall[1]
perm_abs_accel_svm_kappa <- perm_abs_accel_svmresults$overall[2]
perm_abs_accel_svm_AUC <- perm_abs_accel_svmrROC$auc[1]
perm_abs_accel_svm_avg_sens <- mean(perm_abs_accel_svmresults[[4]][1:7])
perm_abs_accel_svm_avg_spec <- mean(perm_abs_accel_svmresults[[4]][8:14])

perm_abs_accel_knn_accuracy <- perm_abs_accel_knn_results$overall[1]
perm_abs_accel_knn_kappa <- perm_abs_accel_knn_results$overall[2]
perm_abs_accel_knn_AUC <- perm_abs_accel_knn_ROC$auc[1]
perm_abs_accel_knn_avg_sens <- mean(perm_abs_accel_knn_results[[4]][1:7])
perm_abs_accel_knn_avg_spec <- mean(perm_abs_accel_knn_results[[4]][8:14])

perm_results$accuracy <- rbind(perm_base_gbm_accuracy,perm_base_rf_accuracy,perm_base_svm_accuracy,perm_base_knn_accuracy,
                               perm_vel_gbm_accuracy,perm_vel_rf_accuracy,perm_vel_svm_accuracy,perm_vel_knn_accuracy,
                               perm_accel_gbm_accuracy, perm_accel_rf_accuracy,perm_accel_svm_accuracy,perm_accel_knn_accuracy,
                               perm_abs_accel_gbm_accuracy,perm_abs_accel_rf_accuracy,perm_abs_accel_svm_accuracy,perm_abs_accel_knn_accuracy)

perm_results$avg_sens <- rbind(perm_base_gbm_avg_sens,perm_base_rf_avg_sens,perm_base_svm_avg_sens,perm_base_knn_avg_sens,
                               perm_vel_gbm_avg_sens, perm_vel_rf_avg_sens, perm_vel_svm_avg_sens, perm_vel_knn_avg_sens,
                               perm_accel_gbm_avg_sens, perm_accel_rf_avg_sens, perm_accel_svm_avg_sens, perm_accel_knn_avg_sens,
                               perm_abs_accel_gbm_avg_sens, perm_abs_accel_rf_avg_sens, perm_abs_accel_svm_avg_sens, perm_abs_accel_knn_avg_sens)

perm_results$avg_spec <- rbind(perm_base_gbm_avg_spec,perm_base_rf_avg_spec,perm_base_svm_avg_spec,perm_base_knn_avg_spec,
                               perm_vel_gbm_avg_spec, perm_vel_rf_avg_spec, perm_vel_svm_avg_spec, perm_vel_knn_avg_spec,
                               perm_accel_gbm_avg_spec, perm_accel_rf_avg_spec, perm_accel_svm_avg_spec, perm_accel_knn_avg_spec,
                               perm_abs_accel_gbm_avg_spec, perm_abs_accel_rf_avg_spec, perm_abs_accel_svm_avg_spec, perm_abs_accel_knn_avg_spec)

perm_results$kappa <- rbind(perm_base_gbm_kappa,perm_base_rf_kappa,perm_base_svm_kappa,perm_base_knn_kappa,
                            perm_vel_gbm_kappa,perm_vel_rf_kappa,perm_vel_svm_kappa,perm_vel_knn_kappa,
                            perm_accel_gbm_kappa,perm_accel_rf_kappa,perm_accel_svm_kappa,perm_accel_knn_kappa,
                            perm_abs_accel_gbm_kappa,perm_abs_accel_rf_kappa,perm_abs_accel_svm_kappa,perm_abs_accel_knn_kappa)

perm_results$AUC <- rbind(perm_base_gbm_AUC,perm_base_rf_AUC,perm_base_svm_AUC,perm_base_knn_AUC,
                          perm_vel_gbm_AUC,perm_vel_rf_AUC,perm_vel_svm_AUC,perm_vel_knn_AUC,
                          perm_accel_gbm_AUC,perm_accel_rf_AUC, perm_accel_svm_AUC,perm_accel_knn_AUC,
                          perm_abs_accel_gbm_AUC,perm_abs_accel_rf_AUC,perm_abs_accel_svm_AUC,perm_abs_accel_knn_AUC)




perm_results <-data.frame(perm_results)
rownames(perm_results)<- c('perm_base_gbm','perm_base_rf','perm_base_svm','perm_base_knn',
                           'perm_vel_gbm','perm_vel_rf','perm_vel_svm','perm_vel_knn',
                           'perm_accel_gbm','perm_accel_rf','perm_accel_svm','perm_accel_knn',
                           'perm_abs_accel_gbm','perm_abs_accel_rf','perm_abs_accel_svm','perm_abs_accel_knn')
perm_results


save.image(file = 'Phase3.RData')
