setwd("C:/Users/Jeremy/Desktop/final_dem/")
getwd()
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
library(gridExtra)
library(AppliedPredictiveModeling)
library(xgboost)
library(dplyr)
library(data.table)
library(plyr)


cl <- makeCluster((detectCores()-1))
registerDoParallel(cl)

#loading relevant information 
load('Phase2.RData')

#conducting monte carlo of 100 iterations

IG_MC_results <-NULL
perm_MC_results <- NULL

## For accuracy, Kappa, the area under the ROC curve, sensitivity and specificity
fiveStats <- function(...) c(multiClassSummary(...),
                             defaultSummary(...))

fitControl <- trainControl(method = 'cv',
                           number = 10,
                           classProbs = TRUE,
                           summaryFunction = fiveStats,
                           verboseIter = TRUE)


bl_gbmGrid <- expand.grid(nrounds = c(1, 10,100),
                          max_depth = c(1, 4,6),
                          eta = c(.1, .4),
                          gamma = 0,
                          colsample_bytree = .7,
                          min_child_weight = 1,
                          subsample = c(.5,.8,1))

#========================RelieF Monte Carlo============================


for (i in 1:100) {
  
  mc_base_data <- rbind(base_train_data, base_test_data)
  base_data<-mc_base_data[sample(nrow(mc_base_data)),]
  base_split <- createDataPartition(mc_base_data$Diagnosis, p = 0.75)[[1]]
  base_train_data <- base_data[base_split,]
  base_test_data<- base_data[-base_split,]
  
  perm_base_bl_gbmFit <-  caret::train (x=base_train_data [,-1],
                                        y=base_train_data[,1],
                                        method = "xgbTree",
                                        tuneGrid = bl_gbmGrid,
                                        trControl = fitControl,
                                        preProcess= c('zv','nzv'),
                                        tuneLength =10,
                                        verbose = FALSE,
                                        metric='logLoss')
  
  perm_base_bl_gbm_results <- confusionMatrix(predict(perm_base_bl_gbmFit, base_test_data), base_test_data[,1])
  perm_base_prob_results <- predict(perm_base_bl_gbmFit, base_test_data, type = "prob")
  perm_base_bl_gbm_ROC <- multiclass.roc(base_test_data[,1], perm_base_prob_results[,1],levels = levels(base_test_data[,1]))
  
  perm_MC_results <- rbind(perm_MC_results, data.frame( perm_base_bl_gbm_results[[3]][1],
                                                        (perm_base_bl_gbm_results[[4]][1]+perm_base_bl_gbm_results[[4]][2]+perm_base_bl_gbm_results[[4]][3])/3,
                                                        (perm_base_bl_gbm_results[[4]][4]+perm_base_bl_gbm_results[[4]][5]+perm_base_bl_gbm_results[[4]][6])/3,
                                                        perm_base_bl_gbm_results[[3]][2], 
                                                        perm_base_bl_gbm_ROC[7]))
  print(i)
}

colnames(perm_MC_results)<- c('accuracy','avg sens','avg spec','kappa','auc')
boxplot(perm_MC_results, main = 'Monte Carlo xGBM relieF Results',
        xlab = 'Metrics', ylab = 'Percentage')

############################################################################
#=================== Information Gain Monte Carlo ==========================
#=== Note that the reassignment of training and test data is required ======
#===========================================================================
############################################################################


df <- read.csv('adnimergedata.csv',header = T)

#=====================Variable Selection====================
df <- df[,c(1:62)]
names(df)

df <- subset(df, select = -c(FLDSTRENG,FSVERSION,
                             COLPROT,ORIGPROT,PTID,SITE,VISCODE,EXAMDATE,CDRSB,DX.bl))

#removing NA values for diagnosis
sum(is.na(df$DX))
df <- df[!is.na(df$DX),]
colnames(df)[50]<- 'Diagnosis'


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

demo <- base[,c('RID','Diagnosis','AGE','PTGENDER','PTEDUCAT','PTETHCAT','PTRACCAT','PTMARRY','APOE4')]
demo <- demo %>% distinct(RID,Diagnosis,.keep_all = T)
demo <- demo[,-c(1,2)]
base_v1 <- base[,c(1,2,10:ncol(base))]
base_v1 <- summaryBy(.~RID+Diagnosis, FUN = quantile,na.rm = T,data = base_v1)
base_v1 <- cbind(base_v1,base_v1)
base_v1 <- base_v1[,-1]
base_v1 <- base_v1[!is.na(base_v1$Diagnosis),]
base_v1$Diagnosis<- as.factor(base_v1$Diagnosis)

base_data<-base_v1[sample(nrow(base_v1)),]
base_split <- createDataPartition(base_v1$Diagnosis, p = 0.75)[[1]]
base_train_data <- base_data[base_split,]
base_test_data<- base_data[-base_split,]

base_IGain<-attrEval(Diagnosis ~ .,data=base_train_data,estimator="InfGain")
base_numberOfAttributes<- length(base_IGain[base_IGain >=0.01])
base_ig_predictors<-names(sort(base_IGain, decreasing=T)[1:base_numberOfAttributes])
base_train_data <- base_train_data[c('Diagnosis',base_ig_predictors)]



base_mvsd_train_data<- base_train_data[base_train_data$Diagnosis !='CN',]
base_cn_train_data <- base_train_data[base_train_data$Diagnosis == 'CN',]
base_mvsd_train_data$Diagnosis <- factor(base_mvsd_train_data$Diagnosis)
base_train_data <- SMOTE(Diagnosis ~., data= base_mvsd_train_data, k=5, perc.over = 90)  
base_train_data <- rbind(base_train_data,base_cn_train_data)

base_train_data <- rfImpute(Diagnosis ~ ., base_train_data)

dummies <- dummyVars(Diagnosis~., data = base_train_data)
data_numeric <-predict(dummies, base_train_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(base_train_data$Diagnosis,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
base_train_data<-data_numeric 
base_train_data<-base_train_data[sample(nrow(base_train_data)),]

base_test_data <- base_test_data[c('Diagnosis',base_ig_predictors)]

base_test_data <- rfImpute(Diagnosis ~ ., base_test_data)
dummies <- dummyVars(Diagnosis~., data = base_test_data)
data_numeric <-predict(dummies, base_test_data)
data_numeric <- as.data.frame(data_numeric )
data_numeric <-data.frame(base_test_data$Diagnosis,data_numeric )
names(data_numeric )[1] <- "Diagnosis"
base_test_data<-data_numeric 



for (i in 1:100) {
  
  mc_base_data <- rbind(base_train_data, base_test_data)
  base_data<-mc_base_data[sample(nrow(mc_base_data)),]
  base_split <- createDataPartition(mc_base_data$Diagnosis, p = 0.75)[[1]]
  base_train_data <- base_data[base_split,]
  base_test_data<- base_data[-base_split,]
  
  base_bl_gbmFit <-  caret::train (x=base_train_data [,-1],
                                   y=base_train_data[,1],
                                   method = "xgbTree",
                                   tuneGrid = bl_gbmGrid,
                                   trControl = fitControl,
                                   preProcess= c('zv','nzv'),
                                   tuneLength =10,
                                   verbose = FALSE,
                                   metric='logLoss')
  
  base_bl_gbm_results <- confusionMatrix(predict(base_bl_gbmFit, base_test_data), base_test_data[,1])
  base_prob_results <- predict(base_bl_gbmFit, base_test_data, type = "prob")
  base_bl_gbm_ROC <- multiclass.roc(base_test_data[,1], base_prob_results[,1],levels = levels(base_test_data[,1]))
  IG_MC_results <- rbind(IG_MC_results, data.frame(base_bl_gbm_results[[3]][1],
                                                   (base_bl_gbm_results[[4]][1]+base_bl_gbm_results[[4]][2]+base_bl_gbm_results[[4]][3])/3,
                                                   (base_bl_gbm_results[[4]][4]+base_bl_gbm_results[[4]][5]+base_bl_gbm_results[[4]][6])/3,
                                                   base_bl_gbm_results[[3]][2], 
                                                   base_bl_gbm_ROC[7]))
  print(i)
}

colnames(IG_MC_results)<- c('accuracy','avg sens','avg spec','kappa','auc')
boxplot(IG_MC_results, main = 'Monte Carlo xGBM IG Results',
        xlab = 'Metrics', ylab = 'Percentage')

save.image(file = 'Phase2_MC.RData')
