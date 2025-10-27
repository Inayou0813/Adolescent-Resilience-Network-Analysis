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
library(officer)
## ---- Import the data ----
##1.Construct the first overall network
Network4 <- Data[, c(
  "Log_H_F1", "Log_H_E1", "Log_H_T1", "Log_H_P1", "Log_H_DHEA", "Log_H_1AG",
  "GP_3", "GP_4", "rHS_6", "HS_7", 
  "rHS_9", "PT_10", "GP_11", "rHS_12", "PT_13", "PT_14",
  "HS_18", "GP_20", 
  "GP_24", "PT_25", "rHS_26"
)]
View(Network4)
##2.Check for missing values
any_missing <- any(is.na(Network4))
cat("数据集中是否存在缺失值：", any_missing, "\n")
###③查看每一列中的缺失值数量
missing_per_column <- colSums(is.na(Network4))
cat("每一列中的缺失值数量：\n")
print(missing_per_column)

 ##3.Assign variable names to each column based on your own dataset
myname4 <- c("F", "E", "T", "P", "DHEA", "1AG",
             "R3", "R4","R6", "R7",
             "R9", "R10", "R11", "R12", "R13", "R14", 
             "R18", "R20",
             "R24", "R25", "R26")
colnames(Network4)<-myname4

  ##4.Group the variables according to the attributes of each scale
feature_group4<-list(
  'Hair steroids' = c(1:6), 
  Resilience=c(7:21)
)
 
 ###5.Set label name
labelname4 <- c(
  "Cortisol", "Cortisone", "Testosterone", "Progesterone", "DHEA", "1-AG",
  "Life goals", "Maturity from setbacks", "No confidant", "Friend support",
  "Unsure who to seek", "Process for growth", "Problem-solving plan", "Keep emotions inside", "Adversity motivation", "Growth from adversity",
  "Seek listener", "Focused effort",
  "Goal setting", "Positive outlook", "Reluctant to share"
)


## ---- Compute the network ----
par(family = "serif")  
# Compute the partial correlation matrix using cor_auto
pcor_layer4 <- qgraph(cor_auto(Network4),   
                   corMethod = "pcor",         
                   sampleSize = nrow(Network4), 
                   groups = feature_group4,
                   nodeNames = labelname4,
                   minimum = "sig",          
                   layout = "spring",
                   repulsion = 25,
                   aspect = FALSE,
                   rotation = 4,
                   layout.par = list(repulse.rad = 200), 
                   area = 10,  
                   details = TRUE,
                   palette = "colorblind",
                   legend = TRUE,
                   legend.cex = 0.3,
                   mar = c(1.5, 5, 2, 4),     
                   details = FALSE, 
                   layoutScale = c(1, 0.8),  
                   layoutOffset = c(-0.12, 0), edge.width = TRUE)  




## ---- Extract edge weights ----
head(pcor_layer4$Edgelist)

library(Matrix)
nodes4 <- unique(c(pcor_layer4$Edgelist$from, pcor_layer4$Edgelist$to))
n4 <- length(nodes4)
adj_mat4 <- matrix(0, n4, n4)
for (i in seq_along(pcor_layer4$Edgelist$weight)) {
  from <- pcor_layer4$Edgelist$from[i]
  to <- pcor_layer4$Edgelist$to[i]
  w <- pcor_layer4$Edgelist$weight[i]
  adj_mat4[from, to] <- w
  adj_mat4[to, from] <- w
}
rownames(adj_mat4) <- colnames(adj_mat4) <- pcor_layer4$Nodes$labels
adj_mat4


rownames(adj_mat4) <- colnames(adj_mat4) <- labelname4
print(adj_mat4)


write.csv(adj_mat4, file = "adjacency_matrix_layer4.csv", row.names = TRUE)

## ---- Calculate centrality indices ----
centrality_auto(pcor_layer4)
centralityPlot(pcor_layer4, include = c("ExpectedInfluence","Strength", "Closeness","Betweenness"),scale = "z-scores",orderBy ="ExpectedInfluence")
centralityPlot(pcor_layer4, include = c("ExpectedInfluence","Strength", "Closeness","Betweenness"),scale = "raw",orderBy ="ExpectedInfluence")
#dev.off()

## ---- Calculate bridge EI indices ----
adj_matrix <- getWmat(pcor_layer4)  

feature_group4<-list(
  'Hair steroids' = c(1:6), 
  Resilience=c(7:21)
)

bridge_results <- bridge(adj_matrix, communities = feature_group4)
EI_1step <- bridge_results$`Bridge Expected Influence (1-step)`
EI_2step <- bridge_results$`Bridge Expected Influence (2-step)`
print(EI_1step)
print(EI_2step)

pdf("Layer4_total_bridge_EI_plot.pdf", width=6, height=6)
plot(bridge_results,
     include = c("Bridge Expected Influence (1-step)","Bridge Expected Influence (2-step)"),
     order = "value",
     zscore = FALSE,
     color = TRUE,           
     colpalette = "Set1")    
dev.off()

pdf("Layer4_total_bridge_EI1_plot.pdf", width=3, height=6)
plot(bridge_results,
     include = c("Bridge Expected Influence (1-step)"),
     order = "value",
     zscore = FALSE,
     color = TRUE,           
     colpalette = "Set1") 
dev.off()

## ---- Accuracy analysis of network structure. ----
library(bootnet)

bpcor_layer4 <- estimateNetwork(Network4,
                               default = "pcor",
                               corMethod = "cor_auto",
                               weighted = TRUE,
                               signed = TRUE,
                               verbose = TRUE)

# Use bootnet to perform bootstrap estimation of network stability and edge confidence intervals.
Results4 <- bootnet(bpcor_layer4, statistics = c("Strength","Closeness","Betweenness","edge","expectedInfluence"), nBoots=1000, nCores=12, caseMin =
                      0.05, caseMax = 0.75, caseN = 10)

plot(Results4, labels = TRUE)

 ##2 Test for differences in edge weights and centrality indices across variables.
plot(Results4, "strength", plot = "difference", family = "serif") #分析节点在强度上的差异检验
plot(Results4, "closeness", plot = "difference") #分析节点在紧密性上的差异检验
plot(Results4, "betweenness", plot = "difference") #分析节点在中介性上的差异检验
plot(Results4, "expectedInfluence", plot = "difference")
plot(Results4, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")#分析边线权重的差异检验
#onlyNonzero用于设定不呈现权重为0的情况
?bootnet

 ##	3.	Test the stability of nodes in centrality indices.
Results3 <- bootnet(bpcor_layer3, statistics=c("ExpectedInfluence","strength","closeness","betweenness"),type = "case", seMin =
                      0.05, caseMax = 0.75, caseN = 10) 
plot(Results3 ,statistics=c("ExpectedInfluence","strength","closeness","betweenness"))#计算节点在中心指标的稳定性
corStability(Results3)#计算CS-coefficients，高于0.25的是好的
