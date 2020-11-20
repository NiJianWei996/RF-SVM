my_file<-"PDB1075.txt"
file<-my_file
seq <- readLines(file) # 读入序列，每个元素存入一行
seq <- seq[seq != ""] # 去除空行
is.anno <- regexpr(">", seq, perl=T) # 正则匹配（regular expression）注释行,是注释行为1，否则为-1
seq.anno <- seq[ which(is.anno == 1) ] # 注释内容
seq.content <- seq[ which(is.anno == -1) ] # 序列内容
##--计算每条序列内容所占的行数，便于后来拼接--##
start <- which(is.anno == 1) # 注释行行号
end <- start[ 2:length(start) ]-1 # 第二条记录注释行到最后一条记录注释行行号减一，即为每条记录结束行号，这里会统计少一行——最后一行的结束未统计
end <- c(end, length(seq) ) # 末尾添加一行：所有序列结束行
distance <- end - start # 每条记录所占行号
index <- 1:length(start) # 生成一个一到记录总个数的向量
index <- rep(index, distance) # 分组标签
seqs <- tapply(seq.content, index, paste, collapse="") # 拼接每条序列内容，返回一个列表，列表每个元素为一条序列的内容
seq.content<-as.character( seqs ) # 将列表转换为向量，向量每个元素为一条序列的内容
seq.len <- nchar(seq.content) # 获得序列长度
seq.ID <- gsub(">(\\w+\\|){3}([A-Za-z0-9.]+)\\|.*", "\\2", seq.anno, perl = T) # 获取序列的ID
result <- data.frame( seq.ID, seq.anno, seq.len, seq.content ) # 组件结果：ID，长度，注释行，序列内容
LL<-as.numeric(nrow(result))
Index<-c(1:LL)
DDPP<-cbind(Index,result[,c(1,3,4)])
write.csv(DDPP,file = "Train_Data.csv",row.names = F)