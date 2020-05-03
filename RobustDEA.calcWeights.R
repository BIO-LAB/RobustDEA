

stat.venn.counts <- function(pvals, deg.ranks,topGene=2000){
  
  pvalues.matrix = pvals
  deg.ranks = deg.ranks
  
  deg.list <- c(as.character(deg.ranks[1:topGene,1]),as.character(deg.ranks[1:topGene,2]),as.character(deg.ranks[1:topGene,3]),as.character(deg.ranks[1:topGene,4]))
  deg.list <- unique(deg.list)
  top.geneRanks <- deg.ranks[1:topGene,]
  deg.pvals <- pvalues.matrix[deg.list,]
  
  venn.matrix = matrix(rep(0,4*length(deg.list)),nrow = length(deg.list),ncol = 4)
  rownames(venn.matrix) = deg.list
  
  for(i in seq(1:length(deg.list))){
    geneName = deg.list[i]
    indexList = apply(top.geneRanks,2,FUN = function(char){grep(geneName,char)})
    venn.matrix[i,] = sapply(indexList, FUN = function(index){length(index)})
    
  }
  
  return(venn.matrix)
 
}





calculate.weights <- function(deg.pvals,venn.matrix){
  
  require(ROCR)
  
  groupList = list(c(1,2,3),c(1,2,4),c(1,3,4),c(2,3,4))
  classLabelG1 = rowSums(venn.matrix[,groupList[[1]]])==3
  classLabelG2 = rowSums(venn.matrix[,groupList[[2]]])==3
  classLabelG3 = rowSums(venn.matrix[,groupList[[3]]])==3
  classLabelG4 = rowSums(venn.matrix[,groupList[[4]]])==3
  
  deg.pvals <- deg.pvals[names(classLabelG1),]
  aucG1 <- calculate.auc(deg.pvals,classLabelG1)
  

  deg.pvals <- deg.pvals[names(classLabelG2),]
  aucG2 <- calculate.auc(deg.pvals,classLabelG2)

  deg.pvals <- deg.pvals[names(classLabelG3),]
  aucG3 <- calculate.auc(deg.pvals,classLabelG3)
 
  deg.pvals <- deg.pvals[names(classLabelG4),]
  aucG4 <- calculate.auc(deg.pvals,classLabelG4)
  
  
  aucList <- rev(c(aucG1[4],aucG2[3],aucG3[2],aucG4[1]))
  aucWeights <- aucList/sum(aucList)
  print(aucWeights)
  
  return(aucWeights)
  
}


calculate.auc <- function(pvals, classLabels){
  
  weights <- rep(0,length(pvals))
  
  for(i in seq(length(pvals))){
    p.vect <- pvals[,i]
    roc.pred <- prediction(1-p.vect,classLabels)
    roc.auc <- performance(roc.pred, 'auc',fpr.stop=1)@y.values[[1]]
    weights[i] <- roc.auc
  }
  
  return(weights)
} 






