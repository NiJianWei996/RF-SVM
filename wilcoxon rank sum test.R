
Class=rep(c(1,0),c(525,550))
tem51_data_11a=which(Class==1)
tem51_data_11b=which(Class==0)
a=top_150[,1]
tem51_train_data_150=Train_SVM[,c('Class',a)]
tem_ww_s=0
tem_ww_p=0
a=which(Class==1)
b=which(Class==0)
for (i in 1:150) {
  x=tem51_train_data_150[a,i+1]
  y=tem51_train_data_150[b,i+1]
  tem_ww=wilcox.test(x,y,alternative = 'less',exact = F,correct = T)
  tem_ww_s[i]=round(tem_ww$statistic,2)
  tem_ww_p[i]=round(tem_ww_p,6)
}
c=colnames(tem51_train_data_150)
c=c[-1]
ww_result01=as.data.frame(cbind(c,tem_ww_s,tem_ww_p))
write.csv(ww_result01,'威尔克森秩和检验结果.csv')
rm(a,b,c)