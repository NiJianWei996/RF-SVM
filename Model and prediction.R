data_1=read.csv('PDB1075pacc.csv')
data_2=read.csv('PDB1075分布.csv')
data_3=read.csv('PDB1075转换1.csv')
data_4=read.csv('PDB1075转换2.csv')
data_5=read.csv('pdb1075_local.csv')
DATA_1<-cbind(data_1,data_2,data_3,data_4,data_5)
colnames(DATA_1)=paste("X",c(1:459),sep = "")
Index=c(1:1075)
Class=rep(c(1,-1),c(525,550))
Train_data=cbind(Index,Class,DATA_1)
write.csv(Train_data,'train_data_all.csv',row.names = F)
############测试机构成####################
data_21=read.csv('PDB186pacc.csv')
data_22=read.csv('PDB186分布.csv')
data_23=read.csv('PDB186转换1.csv')
data_24=read.csv('PDB186转换2.csv')
data_25=read.csv('PDB186_local.csv')
DATA_21<-cbind(data_21,data_22,data_23,data_24,data_25)
colnames(DATA_21)=paste("X",c(1:459),sep = "")
Index=c(1:186)
Class=rep(c(1,-1),c(93,93))
Test_data=cbind(Index,Class,DATA_21)
write.csv(Test_data,'test_data_all.csv',row.names = F)
##################随机森林标准化################
DATA_tem01=Train_data[,3:461]
PDB_mean=apply(DATA_tem01,2,mean)
PDB_sd=apply(DATA_tem01,2, sd)
Train_forest=matrix(NA,1075,459)
for (i in 1:459) {
  Train_forest[,i]=(DATA_tem01[,i]-PDB_mean[i])/PDB_sd[i]
}
Index=c(1:1075)
Class=rep(c(1,0),c(525,550))
colnames(Train_forest)=paste("X",c(1:459),sep = "")
Train_forest=cbind(Index,Class,Train_forest)
write.csv(Train_forest,'train_randomforest.csv',row.names = F)
##########随机森林#########特征选择########
library(ggplot2)
library("randomForest")
library(pROC)
library(ROCR)
library(openxlsx)
######随机森林重要性#######
rate<-0     #设置模型误判率向量初始值
data_tem02<-Train_forest[,-1]
data_tem02<-as.data.frame(data_tem02)
data_tem02[,'Class']=as.factor(data_tem02[,'Class'])
x.data=data_tem02[,-1]
y.data=data_tem02[,1]
for(i in 1:100){
  set.seed(3019)
  rf_train<-randomForest(x.data,y.data,mtry=i,ntree=1000,proximity=TRUE,importance=TRUE)
  rate[i]<-mean(rf_train$err.rate)   #计算基于OOB数据的模型误判率均值
}

min(rate)
which.min(rate)
x1=c(1:100)
rate_1=as.data.frame(cbind(x1,rate))
colnames(rate_1)=c('x','y')
library(ggplot2)
p1=ggplot(rate_1,aes(x=x,y=y))+geom_point(shape=18,size=2)+scale_x_continuous("mtry")+
  scale_y_continuous('error rate')
p1=p1+theme_bw()+ theme(axis.text=element_text(size=12),axis.title=element_text(size=15))
p1=p1+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p1

set.seed(3019)
rf_train<-randomForest(x.data,y.data,mtry=54,ntree=1000)
data_004=as.data.frame(rf_train$err.rate)
type=rep(c("OOB",'postive','negative'),c(1000,1000,1000))
data_004x=rep(c(1:1000),3)
data_004y=c(as.vector(data_004[,1]),as.vector(data_004[,2]),as.vector(data_004[,3]))
data_005=as.data.frame(cbind(as.factor(type),as.numeric(data_004x),as.numeric(data_004y)))
p1=ggplot(data_005,aes(x=data_004x,y=data_004y,shape=type,colour=type))+geom_line()+geom_point(size=0.5,fill='white')+
  scale_x_continuous("tree numbers",limits = c(1,1000),breaks = seq(0,1000,100))+
  scale_y_continuous('error rate')+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p1=p1+theme(legend.position = c(0.9,0.7)) #绘制模型误差与决策树数量关系图  ###38，400
p1
######计算重要性#######
shuju=as.data.frame(Train_forest)
SL01<-1075
aa01<-round(SL01/5)
bb01<-SL01-4*aa01
set.seed(3019)
folds<-split(shuju[,'Index'], sample(rep(1:5, c(aa01,aa01,aa01,aa01,bb01))))
shuju=shuju[,-1]
Gini_01<-0
importance_RF<-NULL
for (i in 1:5) {
  train_data_12 <- as.data.frame(shuju[-folds[[i]],])
  train_data_12x=train_data_12[,-1]
  train_data_12y=as.factor(train_data_12[,1])
  test_data_12 <- as.data.frame(shuju[folds[[i]],])
  test_data_12x=test_data_12[,-1]
  test_data_12y=as.factor(test_data_12[,1])
  RF<-randomForest(train_data_12x,train_data_12y,mtry=54,ntree=400,importance=T)
  im_01 <-importance(x=RF)
  pre <- predict(RF,newdata=test_data_12x)
  Gini_01<-cbind(Gini_01,im_01[,'MeanDecreaseGini'])
  importance_RF[[i]]<-im_01
}
Gini_02<-Gini_01[,-1]
gini_zong=apply(Gini_02,1,mean)
feature_id<-paste("X",c(1:459),sep = "")
aadd<-cbind(feature_id,Gini_02,gini_zong)
tem_names<-c("feature_id","IM1","IM2","IM3","IM4","IM5","IM_zong")
colnames(aadd)<-tem_names
########保存##########
write.csv(rate,file = "随机森林调参——残差.csv",row.names = F)
write.csv(aadd,file = "特征重要性.csv",row.names = F)
write.xlsx(importance_RF,'随机森林重要性综合汇总.xlsx')
########################SVM标准化#################
###############svm标准化函数###################
Standardization_svm_fun<-function(data){
  shuju=data
  a=nrow(data)
  b=ncol(data)
  data_svm_all=matrix(NA,a,b)
  PDB_svm_min=apply(shuju,2,min)
  PDB_svm_max=apply(shuju,2,max)
  for (i in 1:b) {
    data_svm_all[,i]=-1+(1+1)*(shuju[,i]-PDB_svm_min[i])/(PDB_svm_max[i]-PDB_svm_min[i])
  }
  return(data_svm_all)
}
###############训练集SVM标准化###########
train_svm_all=Standardization_svm_fun(DATA_1)
Index=c(1:1075)
Class=rep(c(1,-1),c(525,550))
colnames(train_svm_all)=paste("X",c(1:459),sep = "")
Train_SVM=cbind(Index,Class,train_svm_all)
write.csv(Train_SVM,'svm训练集_all.csv',row.names = F)
##############测试集SVM#############
testdata_svm_all02=Standardization_svm_fun(DATA_21)
Index=c(1:186)
Class=rep(c(1,-1),c(93,93))
colnames(testdata_svm_all02)=paste("X",c(1:459),sep = "")
Test_SVM02=cbind(Index,Class,testdata_svm_all02)
write.csv(Test_SVM02,'svm测试集_all.csv',row.names = F)
top_50<-read.csv(file = "top_50_features.csv",stringsAsFactors = F)
top_100<-read.csv(file = "top_100_features.csv",stringsAsFactors = F)
top_150<-read.csv(file = "top_150_features.csv",stringsAsFactors = F)
top_200<-read.csv(file = "top_200_features.csv",stringsAsFactors = F)
a=top_50[,1]
train_data_50=Train_SVM[,c('Index','Class',a)]
a=top_100[,1]
train_data_100=Train_SVM[,c('Index','Class',a)]
test_data_100=Test_data[,c('Index','Class',a)]
a=top_150[,1]
train_data_150=Train_SVM[,c('Index','Class',a)]
a=top_200[,1]
train_data_200=Train_SVM[,c('Index','Class',a)]
as=nchar(top_150[,1])
a=as.numeric(substring(top_150[,1],2,as))
a1=a[a<=31]
a2=a[a>31&a<=199]
a3=a[a>199&a<=279]
a4=a[a>279&a<=459]
a_1=paste('X',a1,sep = "")
a_2=paste('X',a2,sep = "")
a_3=paste('X',a3,sep = "")
a_4=paste('X',a4,sep = "")
train_data_PseAAC=Train_SVM[,c('Index','Class',a_1)]
train_data_AAD=Train_SVM[,c('Index','Class',a_2)]
train_data_DDS=Train_SVM[,c('Index','Class',a_3)]
train_data_Local_dpp=Train_SVM[,c('Index','Class',a_4)]
write.csv(train_data_PseAAC,'top150_伪氨基酸组成.csv',row.names = F)
write.csv(train_data_AAD,'top150_序列分布.csv',row.names = F)
write.csv(train_data_DDS,'top150_序列转换.csv',row.names = F)
write.csv(train_data_Local_dpp,'top150_Local_dpp.csv',row.names = F)
write.csv(train_data_50,'train_data_50.csv',row.names = F)
write.csv(train_data_100,'train_data_100.csv',row.names = F)
write.csv(train_data_150,'train_data_150.csv',row.names = F)
write.csv(train_data_200,'train_data_200.csv',row.names = F)

##########################选取最有特征集###################
library(openxlsx)
svm_CV_fun<-function(data,cv.1=c(5,10),g.max,g.min,c.max,c.min){
  #包
  #g.max=10
  #g.min=-10
  #c.max=10
  #c.min=-10
  #cv.1=5
  library(e1071)
  library(dplyr)
  library(pROC)
  library(ROCR)
  #基本信息
  shuju=as.data.frame(data)
  shuju_1=as.data.frame(data[,-1])
  Index=shuju[,"Index"]
  a=nrow(shuju)
  SL01=a
  b=ncol(shuju)
  g_1=g.min
  g_2=g.max
  c_1=c.min
  c_2=c.max
  cv_1=cv.1
  d_1=length(c(g_1:g_2))
  d_2=length(c(c_1:c_2))
  row.num_1=abs(g.min-1)
  col.num_1=abs(c.min-1)
  AUC<-matrix(0,nrow = d_1,ncol = d_2)
  Gamma<-0
  Cost<-0
  for (i in g_1:g_2) {
    Gamma[i+row.num_1]=as.character(2^i)
  }
  for (j in c_1:c_2) {
    Cost[j+col.num_1]=as.character(2^j)
  }
  colnames(AUC)<-Cost
  rownames(AUC)<-Gamma
  AUC_f<-rep(0,cv_1)
  Pre_pro=0
  Class_21=0
  ngrids <- 100#设置格点数为100
  TP <- rep(0,ngrids)
  FP <- rep(0,ngrids)
  FN <- rep(0,ngrids)
  TN <- rep(0,ngrids)
  Sn<-rep(0,ngrids)
  Sp<-rep(0,ngrids)
  F_measure<-rep(0,ngrids)
  Pre<-rep(0,ngrids)
  ACC <- rep(0,ngrids)
  MCC<-rep(0,ngrids)
  p0<-0
  Result=NULL
  
  #交叉
  if(cv.1==5){
    aa01<-round(a/5)
    bb01<-SL01-4*aa01
    set.seed(3019)
    folds<-split(Index, sample(rep(1:5, c(aa01,aa01,aa01,aa01,bb01))))
  }
  if(cv.1==10){
    aa01<-round(a/10)
    bb01<-SL01-9*aa01
    set.seed(3019)
    folds<-split(Index, sample(rep(1:10, c(aa01,aa01,aa01,aa01,aa01,aa01,aa01,aa01,aa01,bb01))))
  }
  #SVM参数择优
  for(g in g_1:g_2){
    for(c in c_1:c_2){
      for (i in 1:cv_1) {
        train_data <- shuju_1[-folds[[i]],]
        train.data.x=train_data[,-1]
        train.data.y=as.factor(train_data[,1])
        test_data <- shuju_1[folds[[i]],]
        test.data.x=test_data[,-1]
        test.data.y=as.factor(test_data[,1])
        svmfit.opt <- svm(train.data.x,train.data.y,scale = FALSE,decision.values = TRUE,probability = TRUE,kernel="radial",gamma=2^g,cost=2^c) 
        pred_01 <- predict(svmfit.opt, test.data.x,probability = TRUE,decision.values = TRUE)
        pred_02<-attr(pred_01,"decision.values")
        pred_02 <- 1/(1+exp((-1)*pred_02[,1]))
        pred_03 <- prediction(pred_02,as.numeric(levels(test.data.y)[test.data.y]))
        perf_01 <- performance(pred_03,"tpr","fpr")
        auc_01 <- performance(pred_03,"auc")
        auc_02 <- unlist(slot(auc_01,"y.values"))
        AUC_f[i]<-auc_02
      }
      AUC[g+row.num_1,c+col.num_1]=mean(AUC_f)
    }
  }
  AUC_best=max(AUC)
  max_nrow=as.numeric(which(AUC==AUC[which.max(AUC)],arr.ind=T)[1,1])
  max_ncol=as.numeric(which(AUC==AUC[which.max(AUC)],arr.ind = T)[1,2])
  g_best=Gamma[max_nrow]
  c_best=Cost[max_ncol]
  Result[[1]]=as.data.frame(AUC)
  Result[[2]]=as.data.frame(cbind(AUC_best,max_nrow,max_ncol,g_best,c_best))
  for (i in 1:cv_1) {
    train_data <- shuju_1[-folds[[i]],]
    train.data.x=train_data[,-1]
    train.data.y=as.factor(train_data[,1])
    test_data <- shuju_1[folds[[i]],]
    test.data.x=test_data[,-1]
    test.data.y=as.factor(test_data[,1])
    svmfit.opt <- svm(train.data.x,train.data.y,scale = FALSE,decision.values = TRUE,probability = TRUE,kernel="radial",gamma=g_best,cost=c_best) 
    pred_01 <- predict(svmfit.opt, test.data.x,probability = TRUE,decision.values = TRUE)
    pred_02<-attr(pred_01,"decision.values")
    pred_02 <- 1/(1+exp((-1)*pred_02[,1]))
    Pre_pro<-c(Pre_pro,pred_02)
    Class_test<-as.numeric(levels(test.data.y)[test.data.y])
    Class_21<-c(Class_21,Class_test)
  }
  Pre_pro<-Pre_pro[-1]
  Class_21<-Class_21[-1]
  Confusion_Matrix<-cbind(Pre_pro,Class_21)
  Confusion_Matrix<-as.data.frame(Confusion_Matrix)
  Result[[3]]=Confusion_Matrix
  #阈值选择
  for(i in 1:ngrids){
    p0[i] <- i/ngrids    #选取阈值p0
    y.pred <- as.numeric(1*(Pre_pro > p0[i]))
    Confusion_Matrix<-cbind(Confusion_Matrix,y.pred)
    TP[i] <- sum(nrow(filter(Confusion_Matrix,Class_21==1&y.pred==1)))
    FP[i] <- sum(nrow(filter(Confusion_Matrix,Class_21==-1&y.pred==1)))
    FN[i] <- sum(nrow(filter(Confusion_Matrix,Class_21==1&y.pred==0)))
    TN[i] <- sum(nrow(filter(Confusion_Matrix,Class_21==-1&y.pred==0)))
    Sn[i]<-TP[i]/(TP[i]+FN[i])
    Sp[i]<-TN[i]/(TN[i]+FP[i])
    F_measure[i]<-(2*TP[i])/(2*TP[i]+FN[i]+FP[i])
    Pre[i]<-Sn[i]/(Sn[i]+0.07*(1-Sp[i]))
    ACC[i] <- (TP[i]+TN[i])/(TP[i]+TN[i]+FP[i]+FN[i])
    MCC[i]<-(TP[i]*TN[i]-FP[i]*FN[i])/sqrt((TP[i]+FP[i])*(TP[i]+FN[i])*(TN[i]+FP[i])*(TN[i]+FN[i]))
    Confusion_Matrix<-Confusion_Matrix[,-3]
  }
  ZB<-cbind(p0,TP,TN,FP,FN,Sn,Sp,F_measure,Pre,ACC,MCC)
  Result[[4]]=as.data.frame(ZB)
  return(Result)
}
result_100=svm_CV_fun(data = train_data_100,cv.1 = 5,g.max = 5,g.min = -5,c.max = 5,c.min = -5)
result_150=svm_CV_fun(data = train_data_150,cv.1 = 5,g.max = 5,g.min = -6,c.max = 5,c.min = -6)
result_50=svm_CV_fun(data = train_data_50,cv.1 = 5,g.max = 5,g.min = -6,c.max = 5,c.min = -6)
result_200=svm_CV_fun(data = train_data_200,cv.1 = 5,g.max = 5,g.min = -6,c.max = 5,c.min = -6)
result_top150_PseAAC=svm_CV_fun(data = train_data_PseAAC,cv.1 = 5,g.max = 10,g.min = -10,c.max = 10,c.min = -10)
result_top150_AAD=svm_CV_fun(data = train_data_AAD,cv.1 = 5,g.max = 10,g.min = -10,c.max = 10,c.min = -10)
result_top150_DDS=svm_CV_fun(data = train_data_DDS,cv.1 = 5,g.max = 10,g.min = -10,c.max = 10,c.min = -10)
result_top150_Local_dpp=svm_CV_fun(data = train_data_Local_dpp,cv.1 = 5,g.max = 10,g.min = -10,c.max = 10,c.min = -10)
write.xlsx(result_top150_PseAAC,'result_top150_PseAA.xlsx')
write.xlsx(result_top150_AAD,'result_top150_AAD.xlsx')
write.xlsx(result_top150_DDS,'result_top150_DDS.xlsx')
write.xlsx(result_top150_Local_dpp,'result_top150_Local_dpp.xlsx')
write.xlsx(result_50,"top50_结果.xlsx")
write.xlsx(result_100,'top100_结果.xlsx')
write.xlsx(result_150,'top150_结果.xlsx')
write.xlsx(result_200,'top200_结果.xlsx')
#最优子集150
library(ggplot2)
DATA_ROC=read.csv("特征子集_ROC曲线.csv")
library(RColorBrewer)
DATA_ROC[,3]=as.factor(DATA_ROC[,3])
a=brewer.pal(4,'Set1')
lines.colour=c(a,'#AAAAAA')
p1=ggplot(DATA_ROC,aes(x=FPR,y=TPR,colour=Feature_subset,shape=Feature_subset))+geom_line()+geom_point(size=1.5,fill='white')
p1=p1+theme(legend.position = c(0.8,0.3))+scale_color_manual(values =c("TOP50"="#E41A1C","TOP100"="#377EB8",'TOP150'='#4DAF4A','TOP200'='#984EA3','Line_of_reference'="#AAAAAA"))
p1=p1+scale_x_continuous("False positive rate")+scale_y_continuous('True positive rate')
p1=p1+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p1


DATA_ROC02=read.csv('ROC_分组特征在训练集上.csv')
a_1=brewer.pal(5,'Set1')
lines.colour_1=c(a_1,'#AAAAAA')
p1=ggplot(DATA_ROC02,aes(x=FPR,y=TPR,colour=Feature_type,shape=Feature_type))+geom_line()+geom_point(size=1.5,fill='white')
p1=p1+theme(legend.position = c(0.8,0.3))+scale_color_manual(values =c("Local-DPP"="#E41A1C","PseAAC"="#377EB8",'AAD'='#4DAF4A','DTF'='#984EA3',"TOP150"="#FF7F00",'Line_of_reference'="#AAAAAA"))
p1=p1+scale_x_continuous("False positive rate")+scale_y_continuous('True positive rate')
p1=p1+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p1=p1+theme(legend.position = c(0.8,0.3))+theme(axis.text.x = element_text(face="bold", size=10),axis.text.y = element_text(face="bold", size=10))+
  theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14))
p1

################Test_SVM####PDB186表现######################
testdata_svm_all02=Standardization_svm_fun(DATA_21)
Index=c(1:186)
Class=rep(c(1,-1),c(93,93))
colnames(testdata_svm_all02)=paste("X",c(1:459),sep = "")
Test_SVM02=cbind(Index,Class,testdata_svm_all02)
top_150<-read.csv(file = "top_150_features.csv",stringsAsFactors = F)
a=top_150[,1]
test_data_150=Test_SVM02[,c('Index','Class',a)]
train_data_150=Train_SVM[,c('Index','Class',a)]
library(e1071)
library(dplyr)
library(pROC)
library(ROCR)
train_data_150_02=as.data.frame(train_data_150)
test_data_150_02=as.data.frame(test_data_150)
train.data.x=train_data_150_02[,-c(1,2)]
train.data.y=as.factor(train_data_150_02[,2])
test.data.x=test_data_150_02[,-c(1,2)]
test.data.y=as.factor(test_data_150_02[,2])
svmfit.opt <- svm(train.data.x,train.data.y,scale = FALSE,decision.values = TRUE,probability = TRUE,kernel="radial",gamma=0.0625,cost=2) 
pred_01 <- predict(svmfit.opt, test.data.x,probability = TRUE,decision.values = TRUE)
pred_02<-attr(pred_01,"decision.values")
pred_02 <- 1/(1+exp((-1)*pred_02[,1]))
pred_03 <- prediction(pred_02,as.numeric(levels(test.data.y)[test.data.y]))
perf_01 <- performance(pred_03,"tpr","fpr")
auc_01 <- performance(pred_03,"auc")
auc_02 <- unlist(slot(auc_01,"y.values"))
Pre_pro<-pred_02
Class_21<-as.numeric(levels(test.data.y)[test.data.y])
Confusion_Matrix<-cbind(Pre_pro,Class_21)
Confusion_Matrix<-as.data.frame(Confusion_Matrix)
p0 <- 0.643 #选取阈值p0
y.pred <- as.numeric(1*(Pre_pro > p0))
Confusion_Matrix<-cbind(Confusion_Matrix,y.pred)
TP <- sum(nrow(filter(Confusion_Matrix,Class_21==1&y.pred==1)))
FP <- sum(nrow(filter(Confusion_Matrix,Class_21==-1&y.pred==1)))
FN <- sum(nrow(filter(Confusion_Matrix,Class_21==1&y.pred==0)))
TN <- sum(nrow(filter(Confusion_Matrix,Class_21==-1&y.pred==0)))
Sn<-TP/(TP+FN)
Sp<-TN/(TN+FP)
F_measure<-(2*TP)/(2*TP+FN+FP)
Pre<-Sn/(Sn+0.07*(1-Sp))
ACC <- (TP+TN)/(TP+TN+FP+FN)
MCC<-(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
Confusion_Matrix<-Confusion_Matrix[,-3]
ZB=cbind(TP,FP,FN,TN,Sn,Sp,F_measure,Pre,ACC,MCC)
write.csv(ZB,'PDB186结果.csv',row.names = F)


