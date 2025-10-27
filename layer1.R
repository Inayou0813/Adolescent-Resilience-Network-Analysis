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
library(networktools)
## ---- Import the data ----
##1.Construct the first overall network
Network1 <- Data[, c(
  "Log_H_F1", "Log_H_E1", "Log_H_T1", "Log_H_P1", "Log_H_DHEA", "Log_H_1AG",
  "TR"
)]
View(Network1)

##2.Check for missing values
any_missing <- any(is.na(Network1))
cat("数据集中是否存在缺失值：", any_missing, "\n")
###③查看每一列中的缺失值数量
missing_per_column <- colSums(is.na(Network1))
cat("每一列中的缺失值数量：\n")
print(missing_per_column)

 ##3.Assign variable names to each column based on your own dataset
myname<-c( "F", "E", "T", "P", "DHEA", "1AG",
           "TR")
colnames(Network1)<-myname

  ##4.Group the variables according to the attributes of each scale
feature_group<-list(
  'Hair steroids' = c(1:6), 
  Resilience=c(7)
)
 
 ###5.Set label names
labelname1<-c( "Cortisol", "Cortisone", "Testosterone", "Progesterone", "DHEA", "1-AG", "Resilience")

## ---- Compute the network ----

pdf("zlayer1_total.pdf", family = "Times", width = 7, height = 5)
par(family = "serif") 
#Compute the partial correlation matrix using cor_auto
pcor_layer1 <- qgraph(cor_auto(Network1),   
                   corMethod = "pcor",          
                   sampleSize = nrow(Network1), 
                   groups = feature_group,
                   nodeNames = labelname1,
                   minimum = "sig",            
                   layout = "spring",
                   details = TRUE,
                   palette = "colorblind",
                   legend = TRUE,
                   legend.cex = 0.5,
                   mar = c(1.5, 5, 2, 4),     
                   layoutScale = c(1, 0.83),  
                   layoutOffset = c(-0.12, 0), edge.width = TRUE)
dev.off()

## ---- Extract edge weights ----
head(pcor_layer1$Edgelist)
pcor_layer1$Edgelist
library(Matrix)
nodes <- unique(c(pcor_layer1$Edgelist$from, pcor_layer1$Edgelist$to))
n <- length(nodes)
adj_mat <- matrix(0, n, n)
for (i in seq_along(pcor_layer1$Edgelist$weight)) {
  from <- pcor_layer1$Edgelist$from[i]
  to <- pcor_layer1$Edgelist$to[i]
  w <- pcor_layer1$Edgelist$weight[i]
  adj_mat[from, to] <- w
  adj_mat[to, from] <- w  
}
rownames(adj_mat) <- colnames(adj_mat) <- pcor_layer1$Nodes$labels
adj_mat


rownames(adj_mat) <- colnames(adj_mat) <- labelname1
print(adj_mat)


write.csv(adj_mat, file = "zlayer1_total_adjacency_matrix.csv", row.names = TRUE)

## ---- Calculate centrality indices ----
centrality_auto(pcor_layer1)
centralityPlot(pcor_layer1, include = c("ExpectedInfluence","Strength", "Closeness","Betweenness"),scale = "z-scores",orderBy ="ExpectedInfluence")
centralityPlot(pcor_layer1, include = c("ExpectedInfluence","Strength", "Closeness","Betweenness"),scale = "raw",orderBy ="ExpectedInfluence")
#dev.off()

## ---- Calculate bridge EI indices ----
library(qgraph)
library(networktools)
adj_matrix <- getWmat(pcor_layer1)  

feature_group <- list(
  'Hair steroids' = c(1:6), 
  'Resilience' = c(7)
)

bridge_results <- bridge(adj_matrix, communities = feature_group)
EI_1step <- bridge_results$`Bridge Expected Influence (1-step)`
EI_2step <- bridge_results$`Bridge Expected Influence (2-step)`
print(EI_1step)
print(EI_2step)

pdf("ZLayer1_total_BEI.pdf",family = "Times",width=6, height=6)  
plot(bridge_results,
     include = c("Bridge Expected Influence (1-step)","Bridge Expected Influence (2-step)"),
     order = "value",
     zscore = FALSE,
     color = TRUE,           
     colpalette = "Set1")    
dev.off()

pdf("ZLayer1_total_BEI1.pdf",family = "Times", width=3, height=6)
plot(bridge_results,
     include = c("Bridge Expected Influence (1-step)"),
     order = "value",
     zscore = FALSE,
     color = TRUE,           
     colpalette = "Set1")    
dev.off()

## ---- Accuracy analysis of network structure ----
library(bootnet)

bpcor_layer1 <- estimateNetwork(Network1,
                               default = "pcor",
                               corMethod = "cor_auto",
                               weighted = TRUE,
                               signed = TRUE,
                               verbose = TRUE)

# Test differences between variables in edge weights and centrality indices.
Results1 <- bootnet(bpcor_layer1, statistics = c("Strength","Closeness","Betweenness","ExpectedInfluence","edge"), nBoots=1000, nCores=12, caseMin =
                      0.05, caseMax = 0.75, caseN = 10)
pdf("Layer1_total_bootstrap.pdf",family = "Times", width=7, height=5)  # 设置文件名和尺寸
plot(Results1, labels = TRUE)
dev.off()


##2 Test for differences in edge weights and centrality indices across variables.
plot(Results1, "strength", plot = "difference", family = "serif") #分析节点在强度上的差异检验
plot(Results1, "closeness", plot = "difference") #分析节点在紧密性上的差异检验
plot(Results1, "betweenness", plot = "difference") #分析节点在中介性上的差异检验
plot(Results1, "expectedInfluence", plot = "difference")
plot(Results1, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")#分析边线权重的差异检验
#onlyNonzero用于设定不呈现权重为0的情况

 ##	3.	Test the stability of nodes in centrality indices.
Results2 <- bootnet(bpcor_layer1, statistics=c("ExpectedInfluence","strength","closeness","betweenness"),type = "case", caseMin =
                      0.05, caseMax = 0.75, caseN = 10) 
#type=case就是case-bootstrap，减少样本看结果的稳定性
pdf("ZLayer1_total_CS.pdf",family = "Times", width=7, height=5)  # 设置文件名和尺寸
plot(Results2.1 ,statistics=c("ExpectedInfluence","strength","closeness","betweenness"))#计算节点在中心指标的稳定性
dev.off()
corStability(Results2.1)#计算CS-coefficients，高于0.25的是好的

## ---- Girl ----
Data1 <- read.spss("E:/[2]postgraduate/【6】论文/心理韧性/【01】SPSS_DATA/R_input_data/Girls_withitems.sav",use.value.labels = FALSE,to.data.frame = TRUE)
View(Data1)
## ---- Import the data ----
GNetwork1 <- Data1[, c(
  "Log_H_F1", "Log_H_E1", "Log_H_T1", "Log_H_P1", "Log_H_DHEA", "Log_H_1AG",
  "TR"
)]

GNetwork1_raw <- GNetwork1
View(GNetwork1_raw)

##2.Check for missing values
any_missing <- any(is.na(GNetwork1))
cat("数据集中是否存在缺失值：", any_missing, "\n")
###③查看每一列中的缺失值数量
missing_per_column <- colSums(is.na(GNetwork1))
cat("每一列中的缺失值数量：\n")
print(missing_per_column)

##3.Assign variable names to each column based on your own dataset
myname<-c( "F", "E", "T", "P", "DHEA", "1AG",
           "TR")
colnames(GNetwork1)<-myname
##4.Group the variables according to the attributes of each scale
feature_group<-list(
  'Hair steroids' = c(1:6), 
  Resilience=c(7)
)

###5.Set label name
labelname1<-c( "Cortisol", "Cortisone", "Testosterone", "Progesterone", "DHEA", "1-AG", "Resilience")

## ---- Compute the network ----
# 使用cor_auto计算偏相关矩阵（默认处理混合数据类型）
pdf("zlayer1_girls.pdf", family = "Times", width = 7, height = 5)
Gpcor_layer1 <- qgraph(cor_auto(GNetwork1),   # 使用cor_auto代替cor
                      corMethod = "pcor",          # 指定偏相关（部分相关）
                      sampleSize = nrow(GNetwork1), # 保持样本量参数
                      
                      # 网络图参数
                      groups = feature_group,
                      nodeNames = labelname1,
                      minimum = "sig",             # 仅显示显著连接
                      layout = "spring",
                      details = TRUE,
                      palette = "colorblind",
                      legend = TRUE,
                      legend.cex = 0.3,
                      mar = c(1.5, 5, 2, 4),     # 图形边距
                      
                      # 布局微调
                      layoutScale = c(1, 0.83),  # 调整网络图宽窄比例
                      layoutOffset = c(-0.12, 0), edge.width = TRUE)    # 调整网络图与图例的偏移距离
# 关闭设备，保存PDF
dev.off()

## ---- Extract edge weights ----
# 提取边权矩阵
head(Gpcor_layer1$Edgelist)#Edgelist中有边权值
Gpcor_layer1$Edgelist

#把输出的weight变成矩阵形式
library(Matrix)
nodes <- unique(c(Gpcor_layer1$Edgelist$from, Gpcor_layer1$Edgelist$to))
n <- length(nodes)
adj_mat <- matrix(0, n, n)
for (i in seq_along(Gpcor_layer1$Edgelist$weight)) {
  from <- Gpcor_layer1$Edgelist$from[i]
  to <- Gpcor_layer1$Edgelist$to[i]
  w <- Gpcor_layer1$Edgelist$weight[i]
  adj_mat[from, to] <- w
  adj_mat[to, from] <- w  # 无向图对称赋值
}
rownames(adj_mat) <- colnames(adj_mat) <- Gpcor_layer1$Nodes$labels
adj_mat

#给行列重新命名
rownames(adj_mat) <- colnames(adj_mat) <- labelname1
print(adj_mat)

#写成Excel
write.csv(adj_mat, file = "zlayer1_girls_adjacency_matrix.csv", row.names = TRUE)

## ---- Calculate centrality EI indices ----
centrality_auto(Gpcor_layer1)
centralityPlot(Gpcor_layer1, include = c("ExpectedInfluence","Strength", "Closeness","Betweenness"),scale = "z-scores",orderBy ="ExpectedInfluence")
centralityPlot(Gpcor_layer1, include = c("ExpectedInfluence","Strength", "Closeness","Betweenness"),scale = "raw",orderBy ="ExpectedInfluence")
#dev.off()
#ExpectedInfluence节点期望影响，由于节点强度(Strength)难以很好的聚合负向的边线(edge)，因此研究者对强度的指标进行了一定的调整而得到，这个指标结果在聚合正向的边线时与强度较为相似；
#centralityPlot(list(cor=qgraph_cor, pcor=qgraph_pcor, glasso= qgraph_glasso), include = c("Strength", "Closeness","Betweenness","ExpectedInfluence"),scale = "z-scores")

## ---- Calculate bridge EI indices ----
# 1. 先提取邻接矩阵（边的权重矩阵）
adj_matrix <- getWmat(Gpcor_layer1)  

# 2. 定义社群分组（你已经写好了）
feature_group <- list(
  'Hair steroids' = c(1:6), 
  'Resilience' = c(7)
)

# 3. 计算桥中心性（包括桥EI）
bridge_results <- bridge(adj_matrix, communities = feature_group)
EI_1step <- bridge_results$`Bridge Expected Influence (1-step)`
EI_2step <- bridge_results$`Bridge Expected Influence (2-step)`
print(EI_1step)
print(EI_2step)
# 4. 绘图（画 1&2-step EI，原始值，按值排序，社群上色）
pdf("ZLayer1_girls_BEI_plot.pdf", family = "Times", width=6, height=6)  # 设置文件名和尺寸
plot(bridge_results,
     include = c("Bridge Expected Influence (1-step)","Bridge Expected Influence (2-step)"),
     order = "value",
     zscore = FALSE,
     color = TRUE,           # 按社群上色
     colpalette = "Set1")    # RColorBrewer 调色板名称
dev.off()

# 5. 绘图（只画 1-step EI，原始值，按值排序，社群上色）
pdf("ZLayer1_girls_BEI1_plot.pdf", family = "Times", width=3, height=6)  # 设置文件名和尺寸
plot(bridge_results,
     include = c("Bridge Expected Influence (1-step)"),
     order = "value",
     zscore = FALSE,
     color = TRUE,           # 按社群上色
     colpalette = "Set1")    # RColorBrewer 调色板名称
dev.off()

## ---- Accuracy analysis of network structure ----
##1 计算bootstrap区间
Gbpcor_layer1 <- estimateNetwork(GNetwork1,
                                default = "pcor",
                                corMethod = "cor_auto",
                                weighted = TRUE,
                                signed = TRUE,
                                verbose = TRUE
                              )

# Use bootnet to perform bootstrap estimation of network stability and edge confidence intervals.
GResults1 <- bootnet(Gbpcor_layer1, statistics = c("Strength","Closeness","Betweenness","ExpectedInfluence","edge"), nBoots=1000, nCores=12, caseMin =
                       0.05, caseMax = 0.75, caseN = 10)
# 查看结果
pdf("ZLayer1_girls_bootstrap.pdf", family = "Times", width=7, height=5)  # 设置文件名和尺寸
plot(GResults1, labels = TRUE)
dev.off()
edge_summary_layer1G <- summary(GResults1)


##2 Test for differences in edge weights and centrality indices across variables.
plot(GResults1, "strength", plot = "difference", family = "serif") #分析节点在强度上的差异检验
plot(GResults1, "closeness", plot = "difference") #分析节点在紧密性上的差异检验
plot(GResults1, "betweenness", plot = "difference") #分析节点在中介性上的差异检验
plot(GResults1, "expectedInfluence", plot = "difference")
plot(GResults1, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")#分析边线权重的差异检验
#onlyNonzero用于设定不呈现权重为0的情况

##	3.	Test the stability of nodes in centrality indices.
GResults2 <- bootnet(Gbpcor_layer1, statistics=c("ExpectedInfluence","strength","closeness","betweenness"),type = "case", caseMin = 0.05,
                     caseMax = 0.20, caseN = 20, nBoots = 2500) 
GResults2 <- bootnet(Gbpcor_layer1, statistics=c("ExpectedInfluence","strength","closeness","betweenness"),type = "case", seMin =
                       0.05, caseMax = 0.75, caseN = 10) 
#type=case就是case-bootstrap，减少样本看结果的稳定性
pdf("ZLayer1_girls_CS.pdf", family = "Times", width=7, height=5)  # 设置文件名和尺寸
plot(GResults2 ,statistics=c("ExpectedInfluence","strength","closeness","betweenness"))#计算节点在中心指标的稳定性
dev.off()
corStability(GResults2)#计算CS-coefficients，高于0.25的是好的

GResults2 <- bootnet(Gbpcor_layer1, statistics=c("ExpectedInfluence","strength","closeness","betweenness"),type = "case", nBoots=2500, nCores=12, caseMin =
                      0, caseMax = 0.1, caseN = 20)



## ---- Boy ----
Data2 <- read.spss("E:/[2]postgraduate/【6】论文/心理韧性/【01】SPSS_DATA/R_input_data/Boys_withitems.sav",use.value.labels = FALSE,to.data.frame = TRUE)
View(Data2)
## ---- Import the data ----
##1.构建第一个总网络
BNetwork1 <- Data2[, c(
  "Log_H_F1", "Log_H_E1", "Log_H_T1", "Log_H_P1", "Log_H_DHEA", "Log_H_1AG",
  "TR"
)]

##2.检查是否有缺失值
any_missing <- any(is.na(BNetwork1))
cat("数据集中是否存在缺失值：", any_missing, "\n")
###③查看每一列中的缺失值数量
missing_per_column <- colSums(is.na(BNetwork1))
cat("每一列中的缺失值数量：\n")
print(missing_per_column)

##3.根据自己的数据给每列值赋予名称
myname<-c( "F", "E", "T", "P", "DHEA", "1AG",
           "TR")
colnames(BNetwork1)<-myname
##4.根据每列量表属性进行分组
feature_group<-list(
  'Hair steroids' = c(1:6), 
  Resilience=c(7)
)

###5.设置label name
labelname1<-c( "Cortisol", "Cortisone", "Testosterone", "Progesterone", "DHEA", "1-AG", "Resilience")

## ---- Compute the network ----
# 使用cor_auto计算偏相关矩阵（默认处理混合数据类型）
pdf("zlayer1_boys.pdf", family = "Times", width = 7, height = 5)

# 使用cor_auto计算偏相关矩阵（默认处理混合数据类型）
Bpcor_layer1 <- qgraph(cor_auto(BNetwork1),   # 使用cor_auto代替cor
                       corMethod = "pcor",          # 指定偏相关（部分相关）
                       sampleSize = nrow(BNetwork1), # 保持样本量参数
                       
                       # 网络图参数
                       groups = feature_group,
                       nodeNames = labelname1,
                       minimum = "sig",             # 仅显示显著连接
                       layout = "spring",
                       details = TRUE,
                       palette = "colorblind",
                       legend = TRUE,
                       legend.cex = 0.3,
                       mar = c(1.5, 5, 2, 4),     # 图形边距
                       
                       # 布局微调
                       layoutScale = c(1, 0.83),  # 调整网络图宽窄比例
                       layoutOffset = c(-0.12, 0), edge.width = TRUE)    # 调整网络图与图例的偏移距离
# 关闭设备，保存PDF
dev.off()
## ---- Extract edge weights ----
# 提取边权矩阵
head(Bpcor_layer1$Edgelist)#Edgelist中有边权值
Bpcor_layer1$Edgelist

#把输出的weight变成矩阵形式
library(Matrix)
Bnodes <- unique(c(Bpcor_layer1$Edgelist$from, Bpcor_layer1$Edgelist$to))
Bn <- length(Bnodes)
Badj_mat <- matrix(0, Bn, Bn)
for (i in seq_along(Bpcor_layer1$Edgelist$weight)) {
  from <- Bpcor_layer1$Edgelist$from[i]
  to <- Bpcor_layer1$Edgelist$to[i]
  w <- Bpcor_layer1$Edgelist$weight[i]
  Badj_mat[from, to] <- w
  Badj_mat[to, from] <- w  # 无向图对称赋值
}
rownames(Badj_mat) <- colnames(Badj_mat) <- Bpcor_layer1$Nodes$labels
Badj_mat

#给行列重新命名
rownames(Badj_mat) <- colnames(Badj_mat) <- labelname1
print(Badj_mat)

#写成Excel
write.csv(Badj_mat, file = "zlayer1_boys_adjacency_matrix.csv", row.names = TRUE)

## ---- Calculate centrality EI indices ----
centrality_auto(Bpcor_layer1)
centralityPlot(Bpcor_layer1, include = c("ExpectedInfluence","Strength", "Closeness","Betweenness"),scale = "z-scores",orderBy ="ExpectedInfluence")
centralityPlot(Bpcor_layer1, include = c("ExpectedInfluence","Strength", "Closeness","Betweenness"),scale = "raw",orderBy ="ExpectedInfluence")
#dev.off()
#ExpectedInfluence节点期望影响，由于节点强度(Strength)难以很好的聚合负向的边线(edge)，因此研究者对强度的指标进行了一定的调整而得到，这个指标结果在聚合正向的边线时与强度较为相似；
#centralityPlot(list(cor=qgraph_cor, pcor=qgraph_pcor, glasso= qgraph_glasso), include = c("Strength", "Closeness","Betweenness","ExpectedInfluence"),scale = "z-scores")

## ---- Calculate bridge EI indices ----
# 1. 先提取邻接矩阵（边的权重矩阵）
adj_matrix <- getWmat(Bpcor_layer1)  

# 2. 定义社群分组（你已经写好了）
feature_group <- list(
  'Hair steroids' = c(1:6), 
  'Resilience' = c(7)
)

# 3. 计算桥中心性（包括桥EI）
bridge_results <- bridge(adj_matrix, communities = feature_group)
EI_1step <- bridge_results$`Bridge Expected Influence (1-step)`
EI_2step <- bridge_results$`Bridge Expected Influence (2-step)`
print(EI_1step)
print(EI_2step)
# 4. 绘图（画 1&2-step EI，原始值，按值排序，社群上色）
pdf("ZLayer1_boys_BEI.pdf", family = "Times",width=6, height=6)  # 设置文件名和尺寸
plot(bridge_results,
     include = c("Bridge Expected Influence (1-step)","Bridge Expected Influence (2-step)"),
     order = "value",
     zscore = FALSE,
     color = TRUE,           # 按社群上色
     colpalette = "Set1")    # RColorBrewer 调色板名称
dev.off()

# 5. 绘图（只画 1-step EI，原始值，按值排序，社群上色）
pdf("ZLayer1_boys_BEI1.1.pdf", family = "Times", width=3, height=6)  # 设置文件名和尺寸
plot(bridge_results,
     include = c("Bridge Expected Influence (1-step)"),
     order = "value",
     zscore = FALSE,
     color = TRUE,           # 按社群上色
     colpalette = "Set1")    # RColorBrewer 调色板名称
dev.off()

## ---- Accuracy analysis of network structure ----

library(bootnet)

Bbpcor_layer1 <- estimateNetwork(BNetwork1,
                                 default = "pcor",
                                 corMethod = "cor_auto",
                                 weighted = TRUE,
                                 signed = TRUE,
                                 verbose = TRUE)

# 1.Use bootnet to perform bootstrap estimation of network stability and edge confidence intervals.
BResults1 <- bootnet(Bbpcor_layer1, statistics = c("Strength","Closeness","Betweenness","edge"), nBoots=2500, nCores=12, caseMin =
                       0, caseMax = 0.1, caseN = 20)
BResults1 <- bootnet(Bbpcor_layer1, statistics = c("Strength","Closeness","Betweenness","ExpectedInfluence","edge"), nBoots=1000, nCores=12, caseMin =
                       0.05, caseMax = 0.75, caseN = 10)
# 查看结果
pdf("ZLayer1_boys_bootstrap.pdf", family = "Times", width=7, height=5)  # 设置文件名和尺寸
plot(BResults1, labels = TRUE)
dev.off()

# 查看边权重和置信区间
summary(Results1)


##2 Test for differences in edge weights and centrality indices across variables.
plot(BResults1, "strength", plot = "difference", family = "serif") #分析节点在强度上的差异检验
plot(BResults1, "closeness", plot = "difference") #分析节点在紧密性上的差异检验
plot(BResults1, "betweenness", plot = "difference") #分析节点在中介性上的差异检验
plot(BResults1, "expectedInfluence", plot = "difference")
plot(BResults1, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")#分析边线权重的差异检验
#onlyNonzero用于设定不呈现权重为0的情况
?bootnet

##	3.	Test the stability of nodes in centrality indices.

BResults2 <- bootnet(Bbpcor_layer1, statistics=c("ExpectedInfluence","strength","closeness","betweenness"),type = "case", seMin =
                       0.05, caseMax = 0.75, caseN = 10) 
#type=case就是case-bootstrap，减少样本看结果的稳定性
pdf("ZLayer1_boys_CS.pdf", family = "Times", width=7, height=5)  # 设置文件名和尺寸
plot(BResults2 ,statistics=c("ExpectedInfluence","strength","closeness","betweenness"))#计算节点在中心指标的稳定性
dev.off()

corStability(BResults2)#计算CS-coefficients，高于0.25的是好的


## ---- NCT ----
source("NCT_sources.R")
result <- nct_network_compare(
  data1 = BNetwork1,
  data2 = GNetwork1,
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
