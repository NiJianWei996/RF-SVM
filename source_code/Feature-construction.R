rm(list = ls())
library(stringr)
library(data.table)
library(pROC)
library(readr)
library(openxlsx)
library(e1071)
library(pROC)
library(ROCR)
data_1=read_csv("tr_posi_rawdata.csv")
data_2=read_csv("tr_nega_rawdata.csv")
train_Class=rep(c(1,0),c(4500,4500))
data_3=rbind(data_1,data_2)
AAindex8_1=read.csv('Physicochemical properties of 8 amino acids.csv')
F_seq=data_3
D_LHXZ<-AAindex8_1
SL<-as.numeric(nrow(F_seq))
ABC<-strsplit(as.character(F_seq$seq.content),"",fixed=TRUE)
LENGTH<-matrix(NA,SL,1)
for (i in 1:SL) {
  LENGTH[i,]=lengths(ABC[i])
}
Amino_num<-matrix(NA,SL,20)
zan_b1<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
b1_1<-as.numeric(length(ABC))
for (i in 1:b1_1) {
  zan_a1<-ABC[[i]]
  for (i_1 in 1:20) {
    Amino_num[i,i_1]<-sum(str_count(zan_a1, pattern = zan_b1[i_1]))
  }
}
zan_a2<-cbind(Amino_num,LENGTH)
zan_result1<-matrix(NA,SL,20)
for (i in 1:20) {
  zan_result1[,i]<-zan_a2[,i]/zan_a2[,21]
}
any(is.na(zan_result1))
Result_1<-as.data.frame(zan_result1)
write_csv(Result_1,"train_data_Amino_Fre.csv")
zan_c1<- as.numeric(nrow(D_LHXZ))
zan_result2<-matrix(NA,zan_c1,20)
zan_mean1=0
zan_var1=0
zan_a3<-as.matrix(D_LHXZ[,-1])
zan_result2=matrix(NA,zan_c1,20)
for (i in 1:zan_c1) {
  zan_result2[i,]=t(scale(t(D_LHXZ[i,-1]),center = TRUE, scale = TRUE))
}
any(is.na(zan_result2))
result_1=cbind(D_LHXZ[,1],data.frame(zan_result2))
colnames(result_1)=c('Index',colnames(D_LHXZ)[-1])
zan_a4<-F_seq$seq.content
zan_c2<-as.numeric(nrow(F_seq))
zan_c3=as.numeric(nrow(D_LHXZ))
zan_b2=c(3,5,7,9,11)
zan_c4=length(zan_b2)
zan_c5=11
zan_b2<-c(colnames(D_LHXZ)[-1],"X","U")
Protein_ID=substr(F_seq$seq.ID, 2, nchar(F_seq$seq.ID))
zan_result3<-matrix(NA,SL,zan_c5)
for (i in 1:SL) {
  tem_a1<-zan_a4[i]
  tem_a2<-strsplit(tem_a1, "", fixed=TRUE)
  tem_a3<-as.vector(tem_a2[[1]])
  tem_c1<-as.numeric(length(tem_a3))
  tem_result1=as.data.frame(tem_a3)
  for (i_1 in 1:zan_c3) {
    tem_a4=c(result_1[i_1,-1],0,0)
    tem_a5=tem_a3
    for (i_2 in 1:22) {
      tem_a5<-gsub(zan_b2[i_2],tem_a4[i_2],tem_a5)
    }
    tem_result1=cbind(tem_result1,as.numeric(tem_a5))
  }
  for (i_3 in 1:zan_c5) {
    tem_result2=matrix(NA,tem_c1-i_3,zan_c3)
    for (i_4 in 1:(tem_c1-i_3)) {
      for (i_5 in 1:zan_c3) {
        tem_result2[i_4,i_5]=(tem_result1[i_4,i_5+1]-tem_result1[i_4+i_3,i_5+1])^2
      }
    }
    tem_a6=0
    for (i_6 in 1:tem_c1-i_3) {
      tem_a6[i_6]=sum(tem_result2[i_6,1:zan_c3])/zan_c3
    }
    zan_result3[i,i_3]=sum(tem_a6)/(tem_c1-i_3)
  }
  
}

any(is.na(zan_result3))
zan_result4=data.frame(zan_result3)
result_2=cbind(substring(F_seq$seq.ID,2,nchar(F_seq$seq.ID)),Result_1,zan_result4)
colnames(result_2)=c("Index",paste("X",c(1:31),sep = ""))
windows1=7
weight1=0.13
data_1=Result_1
data_2=zan_result4
data_3=result_2
Class_1=rep(c(1,-1),c(4500,4500))
zan_a1=cbind(data_1,data_2[,1:windows1]*weight1)
zan_b1=apply(zan_a1,1,sum)
zan_a2=matrix(NA,nrow(zan_a1),ncol(zan_a1))
for (i in 1:ncol(zan_a1)) {
  zan_a2[,i]=zan_a1[,i]/zan_b1
}
result_1=as.data.frame(zan_a2)
any(is.na(result_1))
write_csv(result_1,"PseAAC_feature.csv")



data_1=read_csv("tr_posi_rawdata.csv")
data_2=read_csv("tr_nega_rawdata.csv")
data_3=result_1
data_4=rbind(data_1,data_2)
data_5=data_3[,-1]
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
  for (i_1 in 1:A_L) {
    LHXZ_1=as.matrix(AAindex[i_1,])
    for (i_2 in 1:Seq_L) {
      xx_seq=ABC[[i_2]]
      for (i in 1:20) {
        xx_seq<-gsub(animo_names[i],LHXZ_1[i],xx_seq)
      }
      xx_seq<-gsub("X",0,xx_seq)
      xx_seq<-gsub("U",0,xx_seq)
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
result1<-amnio_distribut(data_4,data_5)
Result1=data.frame(result1)
any(is.na(result1))
write_csv(Result1,"AAD_features.csv")
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
result2=Dipeptide_fun(data_4)
any(is.na(result2))
Result2=data.frame(result2)
write_csv(Result2,"Dipeptide_features.csv")

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
result3=Tripeptide_fun(data_4)
any(is.na(result3))
Result3=data.frame(result3)
write_csv(Result3,"Tripeptide_features.csv")
rm(list = ls())
a1=c("PSSM/")
filenames_1=list.files("PSSM",pattern="*.csv")
filenames_2=substr(filenames_1,1,nchar(filenames_1)-4)数
s_1=length(filenames_1)
library(dplyr)
library(readr)
K_1=2
F_1=3
result_3=0
for (i in 1:s_1) {
  zan_b1=paste(a1,filenames_1[i],sep = "")
  zan_a1=read.csv(zan_b1)
  zan_a2=zan_a1[,2:21]
  s_2=nrow(zan_a2)
  zan_a3=as.matrix(zan_a2)
  s_3=ceiling(s_2/K_1)
  s_4=s_2-s_3*(K_1-1)
  zan_a4=NULL
  for (i_1 in 1:(K_1-1)) {
    zan_a4[[i_1]]=zan_a3[(1+(i_1-1)*s_3):(i_1*s_3),]
  }
  zan_a4[[K_1]]=zan_a3[(1+(K_1-1)*s_3):s_2,]
  result_1=0
  for (i_2 in 1:K_1) {
    zan_a5=zan_a4[[i_2]]
    zan_b2=apply(zan_a5,2,mean)
    result_1=c(result_1,zan_b2)
  }
  
  C_2=0
  for (i_3 in 1:K_1) {
    zan_a6=zan_a4[[i_3]]
    s_5=nrow(zan_a6)
    
    
    for (i_5 in 1:20) {
      
      for (i_4 in 1:F_1) {
        C_1=0
        for (i_6 in 1:(s_5-F_1)) {
          C_1[i_6]=(zan_a6[i_6,i_5]-zan_a6[i_6+i_4,i_5])^2
        }
        C_2=c(C_2,mean(C_1))
        
      }
    }
    
    
  }
  result_2=c(result_1[-1],C_2[-1])
  result_3=rbind(result_3,result_2)
}
Result_DPP=result_3[-1,]
any(is.na(Result_DPP))
Result_DPP_1=data.frame(Result_DPP)
write_csv(Result_DPP_1,"Local_DPP_features.csv")

