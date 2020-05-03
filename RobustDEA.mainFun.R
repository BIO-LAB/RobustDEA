

RobustDEA <- function(countsTable, replicates,topGenes = 2000 ){
  

  

  
#  -------------------  run four DE methods------------------------------------  
  
  res.DESeq2 <- run.DESeq2(count.table,replicates)  #results: baseMean,log2FoldChange,lfcSE,stat, pvalue,padj
  res.edgeR <- run.edgeR(count.table,replicates)    # results: logFC,logCPM,PValue
  res.voom <- run.voom(count.table,replicates)       # results: pvalue,padj
  res.NOISeq <- run.NOISeq(count.table,replicates)   #results: 1_mean,2_mean,theta,prob,log2FC. eg.probability of differential expression ("prob")

  
  
  geneNames <- rownames(count.table)
  methodNames <- c('DESeq2','edgeR','Voom','NOISeq')
  
  combGenePvalues <- data.frame(res.DESeq2$pvalue,res.edgeR$PValue,res.voom$pvalue,res.NOISeq$prob)
  rownames(combGenePvalues) <- geneNames  
  colnames(combGenePvalues) <- methodNames

  combDEGRanks<- data.frame(rownames(res.DESeq2[order(res.DESeq2$pvalue),]),rownames(res.edgeR[order(res.edgeR$PValue),]),
                            rownames(res.voom[order(res.voom$pvalue),]),rownames(res.NOISeq[order(res.NOISeq$prob),]))
  colnames(combDEGRanks) <- methodNames
  
  combData <- list(pvalues.matrix = combGenePvalues,  deg.ranks = combDEGRanks )

  
  #  -------------------  combine four DE methods------------------------------------  
  
  
  pvalues.matrix = combData$pvalues.matrix
  deg.ranks = combData$deg.ranks
  geneNames = rownames(pvalues.matrix)
  
  venn.matrix <- stat.venn.counts(pvalues.matrix, deg.ranks, topGenes)
  weights <- calculate.weights(pvalues.matrix, venn.matrix)
  
  p.robust <- combine.robust(pvalues.matrix, weights)
  p.robust <- p.robust$sum
  names(p.robust) <- geneNames
  
  
  return(p.robust)
  
  
}