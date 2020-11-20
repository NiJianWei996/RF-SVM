F_L<-read.csv('8种氨基酸性质.csv',row.names = 1)
F_seq<-read.csv(file = "PDB1075_Train_Data.csv",header = T)
F_Class<-rep(c(1,0),c(525,550))#####Class
library(stringr)
library(data.table)
library(pROC)
library(openxlsx)
library(e1071)
library(pROC)
library(ROCR)
##氨基酸频率##
SL<-as.numeric(nrow(F_seq))
ABC<-strsplit(as.character(F_seq$seq.content),"",fixed=TRUE)
LENGTH<-matrix(NA,SL,1)
for (i in 1:SL) {
  LENGTH[i,]=lengths(ABC[i])
}
a<-lengths(ABC[1])
amino<-matrix(NA,SL,20)
x<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
tem_E<-as.numeric(length(ABC))
for (jj in 1:tem_E) {
  tem_EE<-ABC[[jj]]
  for (i in 1:20) {
    amino[jj,i]<-sum(str_count(tem_EE, pattern = x[i]))
  }
}
zan_data<-cbind(amino,LENGTH)
result04<-matrix(NA,SL,20)
for (i in 1:20) {
  result04[,i]<-zan_data[,i]/zan_data[,21]
}
Result_5<-result04
write.csv(Result_5,"PDB1075_氨基酸频率.csv",row.names=F)
###理化性质###
#####标准化理化性质#######
D_LHXZ<-F_L
#####x_id<-c(1:20)
x_n<- as.numeric(nrow(D_LHXZ))
result01<-matrix(NA,x_n,20)
x_mean<-0
x_DATA01<-as.matrix(D_LHXZ)
for (i in 1:x_n) {
  x_mean[i]<-mean(x_DATA01[i,])
}
x_var<-0
for (i in 1:x_n) {
  x_he<-0
  for (j in 1:20) {
    x_he[j]<-(x_DATA01[i,j]-x_mean[i])^2
  }
  x_var[i]<-sqrt(sum(x_he)/20)
}
for (j in 1:7) {
  for (k in 1:20) {
    result01[j,k]<-(D_LHXZ[j,k]-x_mean[j])/x_var[j]
  }
}
Result_6<-result01###标准化后的理化性质
write.csv(Result_6,"标准化后的理化性质.csv",row.names = F)
######求3到11的seita######
D_seq<-F_seq$seq.content
x_seql<-as.numeric(nrow(F_seq))
x_names<-colnames(D_LHXZ)
Result_21<-NULL
#提取一个理化性质
for (i_2 in 1:x_n) {
  xx_l<-Result_6[i_2,]
  data_t1<-matrix(NA,SL,11)
  #提取一条序列
  for (i_1 in 1:x_seql) {
    xx_seq<-D_seq[i_1]
    xx_seq<-strsplit(xx_seq, "", fixed=TRUE)
    xx_seq<-as.vector(xx_seq[[1]])
    xxx_seql<-as.numeric(length(xx_seq))
    for (i in 1:20) {
      xx_seq<-gsub(x_names[i],xx_l[i],xx_seq)
    }
    xx_seq<-gsub("X",0,xx_seq)
    for (i_3 in 1:11) {
      var_21<-0
      for (i_4 in 1:(xxx_seql-i)) {
        var_21[i_4]<-(as.numeric(xx_seq[i_4+i_3])-as.numeric(xx_seq[i_4]))^2
      }
      var_22<-sum(var_21)/(xxx_seql-i)
      data_t1[i_1,i_3]<-var_22/7
    }
  }
  Result_21[[i_2]]<-data_t1
}
write.xlsx(Result_21,'8种理化性质的R.csv')
####求seita###11个拉姆达
Result_22<-matrix(NA,SL,11)
for (i in 1:11) {
  Result_22[,i]=Result_21[[1]][,i]+Result_21[[2]][,i]+Result_21[[3]][,i]+
    Result_21[[4]][,i]+Result_21[[5]][,i]+Result_21[[6]][,i]+Result_21[[7]][,i]
}
write.csv(Result_22,'不加欧米伽的seita.csv',row.names = F)
#########计算欧米伽########
ww<-seq(0.01,0.2,0.01)#欧米伽取值
Class_1<-rep(c(1,-1),c(525,550))
Index_1<-c(1:SL)
Result_23<-ww
var_23<-0
Var_21=0
#部分参数设置#
#支持向量机计算AUC#
AUC<-matrix(0,nrow = 6,ncol = 7)
AUC_f<-0
#十折交叉验证#
set.seed(3019)
SL01<-SL
aa01<-round(SL01/5)
bb01<-SL01-4*aa01
folds<-split(Index_1, sample(rep(1:5, c(aa01,aa01,aa01,aa01,bb01))))
#支持向量机#
###不同lamada###
for (i_1 in 1:5) {
  lamada=2*i_1+1
  data_s1<-Result_22[,c(1:lamada)]#前lamada个seita
  data_s2<-Result_5#氨基酸频率
  result_21<-matrix(NA,SL,20+lamada)
  result_22<-0
  for (i_2 in 1:length(ww)) {
    a<-ww[i_2]
    tem_A=apply(as.matrix(data_s2[,c(1:20)]), 1,sum)#氨基酸频率和
    tem_B<-matrix(a,lamada,1)#欧米伽
    tem_C<-as.matrix(data_s1)%*%tem_B#欧米伽和
    tem_D<-tem_A+tem_C#总和
    for (i_3 in 1:20) {
      result_21[,i_3]<-data_s2[,i_3]/tem_D[,1]
    }
    for (i_4 in 1:lamada) {
      result_21[,20+i_4]<-a*data_s1[,i_4]/tem_D[,1]
    }
    data_s3<-cbind(Index_1,result_21,Class_1)
    shuju<-data_s3[,-1]
    for(g in 3:8){
      for(c in -3:3){
        for (i in 1:5) {
          train_data <- shuju[as.vector(-folds[[i]]),]
          test_data <- shuju[as.vector(folds[[i]]),]
          svmfit.opt <- svm(train_data[,1:(20+lamada)],train_data[,'Class_1'],scale = FALSE,decision.values = TRUE,probability = TRUE,kernel="radial",gamma=2^g,cost=2^c) 
          pred_01 <- predict(svmfit.opt, test_data[,1:(20+lamada)],probability = TRUE,decision.values = TRUE)
          pred_02<-attr(pred_01,"decision.values")
          pred_02 <- 1/(1+exp((-1)*pred_02[,1]))
          pred_03 <- prediction(pred_02,as.numeric(test_data[,'Class_1']))
          perf_01 <- performance(pred_03,"tpr","fpr")
          auc_01 <- performance(pred_03,"auc")
          auc_02 <- unlist(slot(auc_01,"y.values"))
          AUC_f[i]<-auc_02
        }
        AUC[g-2,c+4]=mean(AUC_f)
      }
    }
    result_22[i_2]<-max(AUC)
  }
  Result_23<-cbind(Result_23,result_22)
}
write.csv(Result_23,'伪氨基酸参数选择.csv',row.names = F)
#########伪氨基酸最优参数###lamada=11，oumiga=0.08#####
############建立训练集测试集的为氨基酸特征#########
paac_fun<-function(data,AAindex,lamda,omiga){
  library(stringr)
  library(data.table)
  Seq_L<-nrow(data)#取蛋白质链数
  AA_L<-nrow(AAindex)#取氨基酸性质个数
  if(ncol(AAindex)>20){AAindex<-AAindex[,-1]}
  F_seq=data
  ABC<-strsplit(as.character(F_seq$seq.content),"",fixed=TRUE)
  LENGTH<-matrix(NA,Seq_L,1)
  for (i in 1:Seq_L) {
    LENGTH[i,]=lengths(ABC[i])
  }
  amino<-matrix(NA,Seq_L,20)
  x<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  tem_E<-as.numeric(length(ABC))
  for (jj in 1:tem_E) {
    tem_EE<-ABC[[jj]]
    for (i in 1:20) {
      amino[jj,i]<-sum(str_count(tem_EE, pattern = x[i]))
    }
  }
  zan_data<-cbind(amino,LENGTH)
  result01<-matrix(NA,Seq_L,20)
  for (i in 1:20) {
    result01[,i]<-zan_data[,i]/zan_data[,21]
  }
  animo_frequency=result01#氨基酸频率 
  rm(tem_E,tem_EE,zan_data)
  #标准化理化性质
  D_LHXZ=AAindex
  result02<-matrix(NA,AA_L,20)
  x_mean<-0
  x_DATA01<-as.matrix(D_LHXZ)
  for (i in 1:AA_L) {
    x_mean[i]<-mean(x_DATA01[i,])
  }
  x_var<-0
  for (i in 1:AA_L) {
    x_he<-0
    for (j in 1:20) {
      x_he[j]<-(x_DATA01[i,j]-x_mean[i])^2
    }
    x_var[i]<-sqrt(sum(x_he)/20)
  }
  for (j in 1:AA_L) {
    for (k in 1:20) {
      result02[j,k]<-(D_LHXZ[j,k]-x_mean[j])/x_var[j]
    }
  }
  AAindex_standar=result02###标准化后的理化性质
  rm(x_mean,x_var,x_he,x_DATA01)
  D_seq<-F_seq$seq.content
  x_seql<-Seq_L
  x_names<-colnames(AAindex)
  Result_21<-NULL
  #提取一个理化性质
  for (i_2 in 1:AA_L) {
    xx_l<-AAindex_standar[i_2,]
    data_t1<-matrix(NA,Seq_L,lamda)
    #提取一条序列
    for (i_1 in 1:x_seql) {
      xx_seq<-D_seq[i_1]
      xx_seq<-strsplit(xx_seq, "", fixed=TRUE)
      xx_seq<-as.vector(xx_seq[[1]])
      xxx_seql<-as.numeric(length(xx_seq))
      for (i in 1:20) {
        xx_seq<-gsub(x_names[i],xx_l[i],xx_seq)
      }
      xx_seq<-gsub("X",0,xx_seq)
      for (i_3 in 1:lamda) {
        var_21<-0
        for (i_4 in 1:(xxx_seql-i)) {
          var_21[i_4]<-(as.numeric(xx_seq[i_4+i_3])-as.numeric(xx_seq[i_4]))^2
        }
        var_22<-sum(var_21)/(xxx_seql-i)
        data_t1[i_1,i_3]<-var_22/7
      }
    }
    Result_21[[i_2]]<-data_t1
  }
  rm(i,i_1,i_2,i_3)
  ####求seita###固定lamda
  Result_22<-matrix(NA,Seq_L,lamda)
  for (i in 1:lamda) {
    Result_22[,i]=Result_21[[1]][,i]+Result_21[[2]][,i]+Result_21[[3]][,i]+
      Result_21[[4]][,i]+Result_21[[5]][,i]+Result_21[[6]][,i]+Result_21[[7]][,i]
  }
  seita=Result_22
  ######将omiga带入#########
  result_21=matrix(NA,Seq_L,20+lamda)
  a=omiga
  tem_A=apply(as.matrix(animo_frequency[,c(1:20)]), 1,sum)#氨基酸频率和
  tem_B<-matrix(a,lamda,1)#欧米伽
  tem_C<-as.matrix(seita)%*%tem_B#欧米伽和
  tem_D<-tem_A+tem_C#总和
  for (i_1 in 1:20) {
    result_21[,i_1]<-animo_frequency[,i_1]/tem_D[,1]
  }
  for (i_2 in 1:lamda) {
    result_21[,20+i_2]<-a*seita[,i_2]/tem_D[,1]
  }
  return(result_21)
}
train_data=read.csv("PDB1075_Train_Data.csv")
test_data=read.csv("PDB186_Test_Data.csv")
AAindex_1=read.csv("8种氨基酸性质.csv",row.names = 1)
train_data_1=paac_fun(train_data,AAindex_1,lamda = 3,omiga = 0.16)
test_data_1=paac_fun(test_data,AAindex_1,lamda = 3,omiga = 0.16)
write.csv(train_data_1,'PDB1075pacc.csv',row.names = F)
write.csv(test_data_1,'PDB186pacc.csv',row.names = F)
#############计算分布###########
amnio_distribut<-function(data,AAindex){
  library(stringr)
  library(data.table)
  Seq_L=nrow(data)
  A_L=nrow(AAindex)
  animo_names=colnames(AAindex)
  result_01=matrix(NA,Seq_L,21)
  result_02=matrix(NA,Seq_L,1)
  F_seq=data
  ABC<-strsplit(as.character(F_seq$seq.content),"",fixed=TRUE)
  LENGTH<-matrix(NA,Seq_L,1)
  for (i in 1:Seq_L) {
    LENGTH[i,]=lengths(ABC[i])
  }
  T_0=rep(0,Seq_L)
  T_25<-round(rep(0.25,Seq_L)*LENGTH)
  T_50<-round(rep(0.5,Seq_L)*LENGTH)
  T_75<-round(rep(0.75,Seq_L)*LENGTH)
  T_1=LENGTH
  T_all=cbind(T_0,T_25,T_50,T_75,T_1)
  ##取一条性质
  for (i_1 in 1:A_L) {
    LHXZ_1=as.matrix(AAindex[i_1,])
    #取一条链#
    for (i_2 in 1:Seq_L) {
      xx_seq=ABC[[i_2]]
      for (i in 1:20) {
        xx_seq<-gsub(animo_names[i],LHXZ_1[i],xx_seq)
      }
      xx_seq<-gsub("X",0,xx_seq)
      var_1=xx_seq
      chain_1=as.numeric(var_1[1:T_all[i_2,2]])
      chain_2=as.numeric(var_1[1:T_all[i_2,3]])
      chain_3=as.numeric(var_1[1:T_all[i_2,4]])
      chain_4=as.numeric(var_1)
      chain_5=as.numeric(var_1[T_all[i_2,2]:T_all[i_2,3]])
      chain_6=as.numeric(var_1[T_all[i_2,3]:T_all[i_2,4]])
      chain_7=as.numeric(var_1[T_all[i_2,4]:T_all[i_2,5]])
      result_01[i_2,1]=mean(chain_1)
      result_01[i_2,2]=sd(chain_1)
      result_01[i_2,3]=sd(chain_1)/mean(chain_1)
      result_01[i_2,4]=mean(chain_2)
      result_01[i_2,5]=sd(chain_2)
      result_01[i_2,6]=sd(chain_2)/mean(chain_2)
      result_01[i_2,7]=mean(chain_3)
      result_01[i_2,8]=sd(chain_3)
      result_01[i_2,9]=sd(chain_3)/mean(chain_3)
      result_01[i_2,10]=mean(chain_4)
      result_01[i_2,11]=sd(chain_4)
      result_01[i_2,12]=sd(chain_4)/mean(chain_4)
      result_01[i_2,13]=mean(chain_5)
      result_01[i_2,14]=sd(chain_5)
      result_01[i_2,15]=sd(chain_5)/mean(chain_5)
      result_01[i_2,16]=mean(chain_6)
      result_01[i_2,17]=sd(chain_6)
      result_01[i_2,18]=sd(chain_6)/mean(chain_6)
      result_01[i_2,19]=mean(chain_7)
      result_01[i_2,20]=sd(chain_7)
      result_01[i_2,21]=sd(chain_7)/mean(chain_7)
    }
    result_02=cbind(result_02,result_01) 
  }
  return(result_02[,-1])
}
train_data_2<-amnio_distribut(train_data,AAindex_1)
test_data_2<-amnio_distribut(test_data,AAindex_1)
write.csv(train_data_2,'PDB1075分布.csv',row.names = F)
write.csv(test_data_2,'PDB186分布.csv',row.names = F)
################转换################
####二肽转换####
Dipeptide_fun<-function(data){
  library(stringr)
  library(data.table)
  F_seq=data
  Seq_L=nrow(F_seq)
  Amino_acid<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  po_ne<-Amino_acid[c(4,7)]
  po_po<-Amino_acid[c(2,9,12)]
  po_uncharged<-Amino_acid[c(3,5,6,8,16,17,19)]
  nonpo<-Amino_acid[c(1,10,11,13,14,15,18,20)]
  tem_list<-list(po_ne,po_po,po_uncharged,nonpo)
  tem_namess<-c("po_ne","po_po","po_uncharged","nonpo")
  d=0
  tem_result1<-NULL
  for (i in 1:4) {
    a<-tem_list[[i]]
    for (j in 1:4) {
      b<-tem_list[[j]]
      num_a<-as.numeric(length(a))
      num_b<-as.numeric(length(b))
      num_ab<-num_a*num_b
      tem_result_1<-0
      for (k in 1:num_a) {
        for (m in 1:num_b) {
          c<-paste(a[k],b[m],sep = "")
          tem_result_1<-c(tem_result_1,c)
        }
      }
      d<-d+1
      tem_result_1<-tem_result_1[-1]
      tem_result1[[d]]<-tem_result_1
    }
  }
  ABC<-strsplit(as.character(F_seq$seq.content), "", fixed=TRUE)
  result_d2_1 <- matrix(NA,Seq_L,16)
  xx_1 <- 0
  for (i in 1:Seq_L) {
    Shuju<-ABC[[i]]
    L <- lengths(ABC)
    x <- vector(length = (L[i]-1))
    for (j in 1:(L[i]-1)) {
      h <- j+1
      x[j] <- paste(Shuju[j],Shuju[h],sep = "")
    }
    for (k_2 in 1:16) {
      ad<-k_2
      a_2 <- as.numeric(length(tem_result1[[ad]]))
      xx_1 <- tem_result1[[ad]]
      aaa_1 <- c(1:a_2)
      for (k_3 in 1:a_2) {
        aaa_1[k_3] <- sum(str_count(x, pattern = xx_1[as.numeric(k_3)]))
      }
      result_d2_1[i,k_2] <- sum(aaa_1)
    }
  }
  Seq_sum<-apply(result_d2_1,1, sum)
  Result_finally=matrix(NA,Seq_L,16)
  for (i_8 in 1:16) {
    Result_finally[,i_8]=result_d2_1[,i_8]/Seq_sum
  }
  return(Result_finally)
}
train_data_3<-Dipeptide_fun(train_data)
test_data_3<-Dipeptide_fun(test_data)
write.csv(train_data_3,'PDB1075转换1.csv',row.names = F)
write.csv(test_data_3,'PDB186转换1.csv',row.names = F)
####三肽转换####
Tripeptide_fun<-function(data){
  library(stringr)
  library(data.table)
  F_seq=data
  Seq_L=nrow(F_seq)
  Amino_acid<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  po_ne<-Amino_acid[c(4,7)]
  po_po<-Amino_acid[c(2,9,12)]
  po_uncharged<-Amino_acid[c(3,5,6,8,16,17,19)]
  nonpo<-Amino_acid[c(1,10,11,13,14,15,18,20)]
  tem_list<-list(po_ne,po_po,po_uncharged,nonpo)
  tem_result_20<-NULL
  tem_result_21<-0
  d<-0
  for (i in 1:4) {
    a<-tem_list[[i]]
    for (j in 1:4) {
      b<-tem_list[[j]]
      for (j_2 in 1:4) {
        c<-tem_list[[j_2]]
        num_a<-as.numeric(length(a))
        num_b<-as.numeric(length(b))
        num_c<-as.numeric(length(c))
        num_abc<-num_a*num_b*num_c
        tem_result_22<-0
        for (k_1 in 1:num_a) {
          for (k_2 in 1:num_b) {
            for (k_3 in 1:num_c) {
              m_1<-as.numeric(k_1)
              m_2<-as.numeric(k_2)
              m_3<-as.numeric(k_3)
              abcc<-paste(a[m_1],b[m_2],c[m_3],sep = "")
              tem_result_22<-c(tem_result_22,abcc)
            }
          }
        }
        d<-d+1
        tem_result_22<-tem_result_22[-1]
        tem_result_20[[d]]<-tem_result_22
      }
    } 
  }
  Result_1=tem_result_20
  result_2 <- matrix(NA,Seq_L,64)
  ABC<-strsplit(as.character(F_seq$seq.content), "", fixed=TRUE)
  xx_1 <- 0
  for (i in 1:Seq_L) {
    Shuju<-ABC[[i]]
    L <- lengths(ABC)
    x <- vector(length = (L[i]-2))
    for (j in 1:(L[i]-2)) {
      h <- j+1
      g <- h+1 
      x[j] <- paste(Shuju[j],Shuju[h],Shuju[g],sep = "")
    }
    for (k_2 in 1:64) {
      ad<-as.numeric(k_2)
      a_2 <- as.numeric(length(tem_result_20[[ad]]))
      xx_1 <- tem_result_20[[ad]]
      aaa_1 <- c(1:a_2)
      for (k_3 in 1:a_2) {
        aaa_1[k_3] <- sum(str_count(x, pattern = xx_1[as.numeric(k_3)]))
      }
      result_2[i,k_2] <- sum(aaa_1)
    }
  }
  Seq_sum<-apply(result_2,1, sum)
  Result_finally=matrix(NA,Seq_L,64)
  for (i_7 in 1:64) {
    Result_finally[,i_7]=result_2[,i_7]/Seq_sum
  }
  return(Result_finally)
}
train_data_4<-Tripeptide_fun(train_data)
test_data_4<-Tripeptide_fun(test_data)
write.csv(train_data_4,'PDB1075转换2.csv',row.names = F)
write.csv(test_data_4,'PDB186转换2.csv',row.names = F)