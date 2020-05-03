rm(list=ls())

source('RobustDEA.runPackages.R')
source('RobustDEA.runCombine.R')
source('RobustDEA.calcWeights.R')
source('RobustDEA.mainFun.R')

load('maqc_count.rda')

count.matrix <- as.matrix(count.matrix)
count.table <- count.matrix[rowSums(count.matrix)>0,]
colnames(count.table) <-NULL

conds <- rep(1:2, each=7)

p.valeus <- RobustDEA(count.table, replicates = conds, topGenes = 500)
