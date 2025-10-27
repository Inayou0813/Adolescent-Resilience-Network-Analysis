setwd("E:/[2]postgraduate/【6】论文/心理韧性")
library(foreign)
Data <- read.spss("E:/[2]postgraduate/【6】论文/心理韧性/【01】SPSS_DATA/R_input_data/First layer.sav",use.value.labels = FALSE,to.data.frame = TRUE)
View(Data)
## ---- Import and install packages ----
#install.packages("readxl")
#install.packages("qgraph")
#install.packages("bootnet")
#install.packages("ggplot2")

library(readxl)
library(qgraph)
library(ggplot2)
library(bootnet)
library(ppcor)
library(openxlsx)
library(officer)
library(flextable)
library(foreign)
library(dplyr)
## ---- Import the data ----
##1.Construct the first overall network
Network3 <- Data[, c(
  "Log_H_F1", "Log_H_E1", "Log_H_T1", "Log_H_P1", "Log_H_DHEA", "Log_H_1AG",
  "AC","GP","PT","HS","FS"
)]
Network3_raw <- Network3   # 备份原始数据
View(Network3_raw)
Network3 <- as.data.frame(scale(Network3))  # 转成z分数
View(Network3)
##2.Check for missing values
any_missing <- any(is.na(Network3))
cat("数据集中是否存在缺失值：", any_missing, "\n")
###③查看每一列中的缺失值数量
missing_per_column <- colSums(is.na(Network3))
cat("每一列中的缺失值数量：\n")
print(missing_per_column)

 ##3.Assign variable names to each column based on your own dataset
myname3<-c( "F", "E", "T", "P", "DHEA", "1AG",
           "AC","GP","PT","HS","FS")
colnames(Network3)<-myname3

  ##4.Group the variables according to the attributes of each scale
feature_group3<-list(
  'Hair steroids' = c(1:6), 
  Resilience=c(7:11)
)
 
 ###5.Set label name
labelname3<-c( "Cortisol", "Cortisone", "Testosterone", "Progesterone", "DHEA", "1-AG", "Affect control", "Goal planning","Positive thinking", "Family support", "Help seeking")

## ---- Compute the network ----
pdf("Zlayer3_total.pdf", family = "Times", width = 7, height = 5)
# Compute the partial correlation matrix using cor_auto
pcor_layer3 <- qgraph(cor_auto(Network3),   
                   corMethod = "pcor",          
                   sampleSize = nrow(Network3), 
                   groups = feature_group3,
                   nodeNames = labelname3,
                   minimum = "sig",           
                   layout = "spring",
                   details = TRUE,
                   palette = "colorblind",
                   legend = TRUE,
                   legend.cex = 0.3,
                   mar = c(1.5, 5, 2, 4),    
                   layoutScale = c(1, 0.83),  
                   layoutOffset = c(-0.12, 0), edge.width = TRUE)    
dev.off()

## ---- Extract edge weights ----
head(pcor_layer3$Edgelist)

library(Matrix)
nodes3 <- unique(c(pcor_layer3$Edgelist$from, pcor_layer3$Edgelist$to))
n3 <- length(nodes3)
adj_mat3 <- matrix(0, n3, n3)
for (i in seq_along(pcor_layer3$Edgelist$weight)) {
  from <- pcor_layer3$Edgelist$from[i]
  to <- pcor_layer3$Edgelist$to[i]
  w <- pcor_layer3$Edgelist$weight[i]
  adj_mat3[from, to] <- w
  adj_mat3[to, from] <- w  
}
rownames(adj_mat3) <- colnames(adj_mat3) <- pcor_layer3$Nodes$labels
adj_mat3


rownames(adj_mat3) <- colnames(adj_mat3) <- labelname3
print(adj_mat3)


write.csv(adj_mat3, file = "adjacency_matrix_layer3.csv", row.names = TRUE)

## ---- Calculate centrality indices ----
centrality_auto(pcor_layer3)
centralityPlot(pcor_layer3, include = c("ExpectedInfluence","Strength", "Closeness","Betweenness"),scale = "z-scores",orderBy ="ExpectedInfluence")
centralityPlot(pcor_layer3, include = c("ExpectedInfluence","Strength", "Closeness","Betweenness"),scale = "raw",orderBy ="ExpectedInfluence")
#dev.off()

## ---- Calculate bridge EI indices ----

adj_matrix <- getWmat(pcor_layer3)  


feature_group3<-list(
  'Hair steroids' = c(1:6), 
  Resilience=c(7:11)
)


bridge_results <- bridge(adj_matrix, communities = feature_group3)
EI_1step <- bridge_results$`Bridge Expected Influence (1-step)`
EI_2step <- bridge_results$`Bridge Expected Influence (2-step)`
print(EI_1step)
print(EI_2step)

pdf("Layer3_total_bridge_EI_plot.pdf", width=6, height=6)
plot(bridge_results,
     include = c("Bridge Expected Influence (1-step)","Bridge Expected Influence (2-step)"),
     order = "value",
     zscore = FALSE,
     color = TRUE,
     colpalette = "Set1") 
dev.off()

pdf("Layer3_total_bridge_EI1_plot.pdf", width=3, height=6)
plot(bridge_results,
     include = c("Bridge Expected Influence (1-step)"),
     order = "value",
     zscore = FALSE,
     color = TRUE,
     colpalette = "Set1")    
dev.off()

## ---- Accuracy analysis of network structure ----
library(bootnet)
bpcor_layer11 <- estimateNetwork(Network1, default = "EBICglasso", tuning = 0.13,corMethod="cor_auto")
library(bootnet)

# 1. 原始数据 Network1 是 data.frame 或 matrix
# 2. 估计偏相关网络（非正则化，使用cor_auto自动相关计算）
bpcor_layer3 <- estimateNetwork(Network3,
                                default = "pcor",
                                corMethod = "cor_auto",
                                weighted = TRUE,
                                signed = TRUE,
                                verbose = TRUE)

# Test differences between variables in edge weights and centrality indices.
Results3 <- bootnet(bpcor_layer3, statistics = c("Strength","Closeness","Betweenness","expectedInfluence","edge"), nBoots=1000, nCores=12, caseMin =
                      0.05, caseMax = 0.75, caseN = 10)

plot(Results3, labels = TRUE)


 ##2 Test for differences in edge weights and centrality indices across variables.
plot(Results3, "strength", plot = "difference", family = "serif") #分析节点在强度上的差异检验
plot(Results3, "closeness", plot = "difference") #分析节点在紧密性上的差异检验
plot(Results3, "betweenness", plot = "difference") #分析节点在中介性上的差异检验
plot(Results3, "expectedInfluence", plot = "difference")
plot(Results3, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")#分析边线权重的差异检验
#onlyNonzero用于设定不呈现权重为0的情况
?bootnet

 ##	3.	Test the stability of nodes in centrality indices.
Results3 <- bootnet(bpcor_layer3, statistics=c("ExpectedInfluence","strength","closeness","betweenness"),type = "case", seMin =
                      0.05, caseMax = 0.75, caseN = 10) 
#type=case就是case-bootstrap，减少样本看结果的稳定性
plot(Results3 ,statistics=c("ExpectedInfluence","strength","closeness","betweenness"))#计算节点在中心指标的稳定性
corStability(Results3)#计算CS-coefficients，高于0.25的是好的


## ---- Girl ----
Data1 <- read.spss("E:/[2]postgraduate/【6】论文/心理韧性/【01】SPSS_DATA/R_input_data/Girls_withitems.sav",use.value.labels = FALSE,to.data.frame = TRUE)
View(Data1)

## ---- Import the data ----
##1.构建第一个总网络
GNetwork3 <- Data1[, c(
  "Log_H_F1", "Log_H_E1", "Log_H_T1", "Log_H_P1", "Log_H_DHEA", "Log_H_1AG",
  "AC","GP","PT","HS","FS"
)]
View(GNetwork3)
##2.检查是否有缺失值
any_missing <- any(is.na(GNetwork3))
cat("数据集中是否存在缺失值：", any_missing, "\n")
###③查看每一列中的缺失值数量
missing_per_column <- colSums(is.na(GNetwork3))
cat("每一列中的缺失值数量：\n")
print(missing_per_column)

##3.根据自己的数据给每列值赋予名称
myname3<-c( "F", "E", "T", "P", "DHEA", "1AG",
            "AC","GP","PT","HS","FS")
colnames(GNetwork3)<-myname3

##4.根据每列量表属性进行分组
feature_group3<-list(
  'Hair steroids' = c(1:6), 
  Resilience=c(7:11)
)

###5.设置label name
labelname3<-c( "Cortisol", "Cortisone", "Testosterone", "Progesterone", "DHEA", "1-AG", "Affect control", "Goal planning","Positive thinking", "Family support", "Help seeking")

## ---- Compute the network ----
# 使用cor_auto计算偏相关矩阵（默认处理混合数据类型）
Gpcor_layer3 <- qgraph(cor_auto(GNetwork3),   # 使用cor_auto代替cor
                      corMethod = "pcor",          # 指定偏相关（部分相关）
                      sampleSize = nrow(GNetwork3), # 保持样本量参数
                      
                      # 网络图参数
                      groups = feature_group3,
                      nodeNames = labelname3,
                      minimum = "sig",             # 仅显示显著连接
                      layout = "spring",
                      details = TRUE,
                      palette = "colorblind",
                      legend = TRUE,
                      legend.cex = 0.3,
                      mar = c(1.5, 5, 2, 4),     # 图形边距
                      
                      # 布局微调
                      layoutScale = c(0.95, 0.9),  # 调整网络图宽窄比例
                      layoutOffset = c(-0.12, 0), edge.width = TRUE)    # 调整网络图与图例的偏移距离

## ---- Extract edge weights ----
# 提取边权矩阵
head(Gpcor_layer3$Edgelist)#Edgelist中有边权值


#把输出的weight变成矩阵形式
library(Matrix)
Gnodes3 <- unique(c(Gpcor_layer3$Edgelist$from, Gpcor_layer3$Edgelist$to))
Gn3 <- length(Gnodes3)
Gadj_mat3 <- matrix(0, Gn3, Gn3)
for (i in seq_along(Gpcor_layer3$Edgelist$weight)) {
  from <- Gpcor_layer3$Edgelist$from[i]
  to <- Gpcor_layer3$Edgelist$to[i]
  w <- Gpcor_layer3$Edgelist$weight[i]
  Gadj_mat3[from, to] <- w
  Gadj_mat3[to, from] <- w  # 无向图对称赋值
}
rownames(Gadj_mat3) <- colnames(Gadj_mat3) <- Gpcor_layer3$Nodes$labels
Gadj_mat3

#给行列重新命名
rownames(Gadj_mat3) <- colnames(Gadj_mat3) <- labelname3
print(Gadj_mat3)

#写成Excel
write.csv(Gadj_mat3, file = "Girls_adjacency_matrix_layer3.csv", row.names = TRUE)

## ---- Calculate centrality indices ----
centrality_auto(Gpcor_layer3)
centralityPlot(Gpcor_layer3, include = c("ExpectedInfluence","Strength", "Closeness","Betweenness"),scale = "z-scores",orderBy ="ExpectedInfluence")
centralityPlot(Gpcor_layer3, include = c("ExpectedInfluence","Strength", "Closeness","Betweenness"),scale = "raw",orderBy ="ExpectedInfluence")
#dev.off()
#ExpectedInfluence节点期望影响，由于节点强度(Strength)难以很好的聚合负向的边线(edge)，因此研究者对强度的指标进行了一定的调整而得到，这个指标结果在聚合正向的边线时与强度较为相似；
#centralityPlot(list(cor=qgraph_cor, pcor=qgraph_pcor, glasso= qgraph_glasso), include = c("Strength", "Closeness","Betweenness","ExpectedInfluence"),scale = "z-scores")

## ---- Calculate bridge EI indices ----
# 1. 先提取邻接矩阵（边的权重矩阵）
adj_matrix <- getWmat(Gpcor_layer3)  

# 2. 定义社群分组（你已经写好了）
feature_group3<-list(
  'Hair steroids' = c(1:6), 
  Resilience=c(7:11)
)

# 3. 计算桥中心性（包括桥EI）
bridge_results <- bridge(adj_matrix, communities = feature_group3)
EI_1step <- bridge_results$`Bridge Expected Influence (1-step)`
EI_2step <- bridge_results$`Bridge Expected Influence (2-step)`
print(EI_1step)
print(EI_2step)
# 4. 绘图（画 1&2-step EI，原始值，按值排序，社群上色）
pdf("Layer3_girls_bridge_EI_plot.pdf", width=6, height=6)  # 设置文件名和尺寸
plot(bridge_results,
     include = c("Bridge Expected Influence (1-step)","Bridge Expected Influence (2-step)"),
     order = "value",
     zscore = FALSE,
     color = TRUE,           # 按社群上色
     colpalette = "Set1")    # RColorBrewer 调色板名称
dev.off()

# 5. 绘图（只画 1-step EI，原始值，按值排序，社群上色）
pdf("Layer3_girls_bridge_EI1_plot.pdf", width=3, height=6)  # 设置文件名和尺寸
plot(bridge_results,
     include = c("Bridge Expected Influence (1-step)"),
     order = "value",
     zscore = FALSE,
     color = TRUE,           # 按社群上色
     colpalette = "Set1")    # RColorBrewer 调色板名称
dev.off()

## ---- Accuracy analysis of network structure ----
##1 计算bootstrap区间
library(bootnet)

# 1. 原始数据 Network1 是 data.frame 或 matrix
# 2. 估计偏相关网络（非正则化，使用cor_auto自动相关计算）
Gbpcor_layer3 <- estimateNetwork(GNetwork3,
                                default = "pcor",
                                corMethod = "cor_auto",
                                weighted = TRUE,
                                signed = TRUE,
                                verbose = TRUE)

# 用 bootnet 做bootstrap估计网络稳定性和边置信区间
GResults3.1 <- bootnet(Gbpcor_layer3, statistics = c("Strength","Closeness","Betweenness","expectedInfluence","edge"), nBoots=1000, nCores=12, caseMin =
                       0.05, caseMax = 0.75, caseN = 10)


# 查看结果
plot(GResults3.1, labels = TRUE)

##2 不同变量在边线权重和中心指标上的差异检验
plot(GResults3.1, "strength", plot = "difference", family = "serif") #分析节点在强度上的差异检验
plot(GResults3.1, "closeness", plot = "difference") #分析节点在紧密性上的差异检验
plot(GResults3.1, "betweenness", plot = "difference") #分析节点在中介性上的差异检验
plot(GResults3.1, "expectedInfluence", plot = "difference")
plot(GResults3.1, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")#分析边线权重的差异检验

#onlyNonzero用于设定不呈现权重为0的情况
?bootnet
##3检验节点在中心指标的稳定性

GResults3.1 <- bootnet(Gbpcor_layer3, statistics=c("ExpectedInfluence","strength","closeness","betweenness"),type = "case", seMin =
                       0.05, caseMax = 0.75, caseN = 10) 
#type=case就是case-bootstrap，减少样本看结果的稳定性
plot(GResults3.1 ,statistics=c("ExpectedInfluence","strength","closeness","betweenness"))#计算节点在中心指标的稳定性
corStability(GResults3.1)#计算CS-coefficients，高于0.25的是好的

## ---- Boy ----
Data2 <- read.spss("E:/[2]postgraduate/【6】论文/心理韧性/【01】SPSS_DATA/R_input_data/Boys_withitems.sav",use.value.labels = FALSE,to.data.frame = TRUE)
View(Data2)

## ---- Import the data ----
##1.构建第一个总网络
BNetwork3 <- Data2[, c(
  "Log_H_F1", "Log_H_E1", "Log_H_T1", "Log_H_P1", "Log_H_DHEA", "Log_H_1AG",
  "AC","GP","PT","HS","FS"
)]
View(BNetwork3)
##2.检查是否有缺失值
any_missing <- any(is.na(BNetwork3))
cat("数据集中是否存在缺失值：", any_missing, "\n")
###③查看每一列中的缺失值数量
missing_per_column <- colSums(is.na(BNetwork3))
cat("每一列中的缺失值数量：\n")
print(missing_per_column)

##3.根据自己的数据给每列值赋予名称
myname3<-c( "F", "E", "T", "P", "DHEA", "1AG",
            "AC","GP","PT","HS","FS")
colnames(BNetwork3)<-myname3

##4.根据每列量表属性进行分组
feature_group3<-list(
  'Hair steroids' = c(1:6), 
  Resilience=c(7:11)
)

###5.设置label name
labelname3<-c( "Cortisol", "Cortisone", "Testosterone", "Progesterone", "DHEA", "1-AG", "Affect control", "Goal planning","Positive thinking", "Family support", "Help seeking")

## ---- Compute the network ----
# 使用cor_auto计算偏相关矩阵（默认处理混合数据类型）
Bpcor_layer3 <- qgraph(cor_auto(BNetwork3),   # 使用cor_auto代替cor
                       corMethod = "pcor",          # 指定偏相关（部分相关）
                       sampleSize = nrow(BNetwork3), # 保持样本量参数
                       
                       # 网络图参数
                       groups = feature_group3,
                       nodeNames = labelname3,
                       minimum = "sig",             # 仅显示显著连接
                       layout = "spring",
                       details = TRUE,
                       palette = "colorblind",
                       legend = TRUE,
                       legend.cex = 0.3,
                       mar = c(1.5, 5, 2, 4),     # 图形边距
                       
                       # 布局微调
                       layoutScale = c(0.95, 0.9),  # 调整网络图宽窄比例
                       layoutOffset = c(-0.12, 0), edge.width = TRUE)    # 调整网络图与图例的偏移距离

## ---- Extract edge weights ----
# 提取边权矩阵
head(Gpcor_layer3$Edgelist)#Edgelist中有边权值


#把输出的weight变成矩阵形式
library(Matrix)
Bnodes3 <- unique(c(Bpcor_layer3$Edgelist$from, Bpcor_layer3$Edgelist$to))
Bn3 <- length(Bnodes3)
Badj_mat3 <- matrix(0, Bn3, Bn3)
for (i in seq_along(Bpcor_layer3$Edgelist$weight)) {
  from <- Bpcor_layer3$Edgelist$from[i]
  to <- Bpcor_layer3$Edgelist$to[i]
  w <- Bpcor_layer3$Edgelist$weight[i]
  Badj_mat3[from, to] <- w
  Badj_mat3[to, from] <- w  # 无向图对称赋值
}
rownames(Badj_mat3) <- colnames(Badj_mat3) <- Bpcor_layer3$Nodes$labels
Badj_mat3

#给行列重新命名
rownames(Badj_mat3) <- colnames(Badj_mat3) <- labelname3
print(Badj_mat3)

#写成Excel
write.csv(Badj_mat3, file = "Boys_adjacency_matrix_layer3.csv", row.names = TRUE)

## ---- Calculate centrality indices ----
centrality_auto(Bpcor_layer3)
centralityPlot(Bpcor_layer3, include = c("ExpectedInfluence","Strength", "Closeness","Betweenness"),scale = "z-scores",orderBy ="ExpectedInfluence")
centralityPlot(Bpcor_layer3, include = c("ExpectedInfluence","Strength", "Closeness","Betweenness"),scale = "raw",orderBy ="ExpectedInfluence")
#dev.off()
#ExpectedInfluence节点期望影响，由于节点强度(Strength)难以很好的聚合负向的边线(edge)，因此研究者对强度的指标进行了一定的调整而得到，这个指标结果在聚合正向的边线时与强度较为相似；
#centralityPlot(list(cor=qgraph_cor, pcor=qgraph_pcor, glasso= qgraph_glasso), include = c("Strength", "Closeness","Betweenness","ExpectedInfluence"),scale = "z-scores")

## ---- Calculate bridge EI indices ----
# 1. 先提取邻接矩阵（边的权重矩阵）
adj_matrix <- getWmat(Bpcor_layer3)  

# 2. 定义社群分组（你已经写好了）
feature_group3<-list(
  'Hair steroids' = c(1:6), 
  Resilience=c(7:11)
)

# 3. 计算桥中心性（包括桥EI）
bridge_results <- bridge(adj_matrix, communities = feature_group3)
EI_1step <- bridge_results$`Bridge Expected Influence (1-step)`
EI_2step <- bridge_results$`Bridge Expected Influence (2-step)`
print(EI_1step)
print(EI_2step)
# 4. 绘图（画 1&2-step EI，原始值，按值排序，社群上色）
pdf("Layer3_boys_bridge_EI_plot.pdf", width=6, height=6)  # 设置文件名和尺寸
plot(bridge_results,
     include = c("Bridge Expected Influence (1-step)","Bridge Expected Influence (2-step)"),
     order = "value",
     zscore = FALSE,
     color = TRUE,           # 按社群上色
     colpalette = "Set1")    # RColorBrewer 调色板名称
dev.off()

# 5. 绘图（只画 1-step EI，原始值，按值排序，社群上色）
pdf("Layer3_boys_bridge_EI1_plot.pdf", width=3, height=6)  # 设置文件名和尺寸
plot(bridge_results,
     include = c("Bridge Expected Influence (1-step)"),
     order = "value",
     zscore = FALSE,
     color = TRUE,           # 按社群上色
     colpalette = "Set1")    # RColorBrewer 调色板名称
dev.off()

## ---- Accuracy analysis of network structure ----
##1 计算bootstrap区间
library(bootnet)

# 1. 原始数据 Network1 是 data.frame 或 matrix
# 2. 估计偏相关网络（非正则化，使用cor_auto自动相关计算）
Bbpcor_layer3 <- estimateNetwork(BNetwork3,
                                 default = "pcor",
                                 corMethod = "cor_auto",
                                 weighted = TRUE,
                                 signed = TRUE,
                                 verbose = TRUE)

# 用 bootnet 做bootstrap估计网络稳定性和边置信区间
BResults3 <- bootnet(Bbpcor_layer3, statistics = c("Strength","Closeness","Betweenness","edge"), nBoots=2500, nCores=12, caseMin =
                       0, caseMax = 0.1, caseN = 20)
BResults3.1 <- bootnet(Bbpcor_layer3, statistics = c("Strength","Closeness","Betweenness","expectedInfluence","edge"), nBoots=1000, nCores=12, caseMin =
                       0.05, caseMax = 0.75, caseN = 10)

# 查看结果
plot(BResults3.1, labels = TRUE)



##2 不同变量在边线权重和中心指标上的差异检验
plot(BResults3.1, "strength", plot = "difference", family = "serif") #分析节点在强度上的差异检验
plot(BResults3.1, "closeness", plot = "difference") #分析节点在紧密性上的差异检验
plot(BResults3.1, "betweenness", plot = "difference") #分析节点在中介性上的差异检验
plot(BResults3.1, "expectedInfluence", plot = "difference")
plot(BResults3.1, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")#分析边线权重的差异检验

##3检验节点在中心指标的稳定性
BResults3.1 <- bootnet(Bbpcor_layer3, statistics=c("ExpectedInfluence","strength","closeness","betweenness"),type = "case", seMin =
                       0.05, caseMax = 0.75, caseN = 10) 
#type=case就是case-bootstrap，减少样本看结果的稳定性
plot(BResults3.1 ,statistics=c("ExpectedInfluence","strength","closeness","betweenness"))#计算节点在中心指标的稳定性
corStability(BResults3.1)#计算CS-coefficients，高于0.25的是好的

## ---- NCT ----
source("NCT_sources.R")
result <- nct_network_compare(
  data1 = BNetwork3,
  data2 = GNetwork3,
  n_perms = 1000,              # 你想的置换次数
  seed = 123,                  # 随机种子，保证结果可复现
  network_method = "partial", # 或 "pearson" 或 "partial"
  p_adjust_method = "BH",      # 多重检验校正方法
  progress = TRUE
)

#输出结果
# 查看整体差异的p值
print(result$global_strength_p)

# 查看边差异的p值矩阵
print(result$edge_pvals)

# 查看调整后的p值矩阵
print(result$edge_pvals_adj)

#结构差异的 p 值
print(result$structure_p)

# 查看网络方法
print(result$network_method)

str(result)  # 显示结果的结构和内容
