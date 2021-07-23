rm(list = ls())
library(e1071)
library(dplyr)
library(pROC)
library(ROCR)
library(readr)
str_data=read_csv("st_train_all_fature_meansd.csv")
DATA_1=read_csv("UniSwiss-Tst_allfeature.csv")
Feature_id=read_csv("174_feature.csv")
DATA_2=read_csv("st_train_all_feature.csv")
data_10=DATA_1[,-1]
data_11=apply(data_10, 2, scale)
data_12=data_11[,Feature_id$a2]
data_21=DATA_2[,Feature_id$a2]
data_13=cbind(DATA_1[,1],data_12)
data_22=cbind(DATA_2[,2],data_21)
train_data <- data_22
train.data.x=train_data[,-1]
train.data.y=as.factor(train_data[,1])
test_data <- data_13
test.data.x=test_data[,-1]
test.data.y=as.factor(test_data[,1])
svmfit.opt <- svm(train.data.x,train.data.y,scale = FALSE,decision.values = TRUE,probability = TRUE,kernel="radial",gamma=2^-7,cost=2^4) 
pred_01 <- predict(svmfit.opt, test.data.x,probability = TRUE,decision.values = TRUE)
pred_02<-attr(pred_01,"decision.values")
pred_02 <- 1/(1+exp((-1)*pred_02[,1]))
Pre_pro<- pred_02
Class_test<-as.numeric(levels(test.data.y)[test.data.y])
Class_pro = Class_test
save(svmfit.opt, file='RF-SVM_Predictor.RData')
load('RF-SVM_Predictor.RData')
