
# edgeR  ------------------------------------------------------------------------#

run.edgeR <- function(countsTable,conds){
  
  library(edgeR)
  

  libSizes <- colSums(countsTable)
  
  dg1 <- DGEList(counts=countsTable,group=conds,lib.size=libSizes)
  dg1 <- calcNormFactors(dg1)
  dg1 <- estimateCommonDisp(dg1)
  dg1 <- estimateTagwiseDisp(dg1)
  res.edgeR <- exactTest(dg1,dispersion='auto')
  
  return(res.edgeR$table)
  
}



# voom ----------------------------------------------------------------------------------#

run.voom <- function(countsTable,conds){
  
  require(limma)
  require(edgeR)
  

  nGenes <- dim(countsTable)[1]
  nSamples <-dim(countsTable)[2]
  libSizes <- colSums(countsTable)
  
  nf = calcNormFactors(countsTable, method = "TMM")
  voom.data = voom(countsTable, design = model.matrix(~factor(conds)), lib.size = libSizes * nf)
  voom.data$genes = rownames(countsTable)
  voom.fitlimma = lmFit(voom.data, design = model.matrix(~factor(conds)))
  voom.fitbayes = eBayes(voom.fitlimma)
  voom.pvalues = voom.fitbayes$p.value[, 2]
  voom.adjpvalues = p.adjust(voom.pvalues, method = "BH")
  voom.results = cbind(voom.pvalues,voom.adjpvalues)
  colnames(voom.results) <- c('pvalue','padj')
  
  return(data.frame(voom.results))
}


# NOISeq - ----------------------------------------------------------------------------------#

 run.NOISeq <- function(countsTable,conds){
   
   require(edgeR)
   require(NOISeq)
   
   normFactor <- calcNormFactors(countsTable)
   libSizes <- colSums(countsTable)
   commLibSizes <- prod(libSizes^(1/length(libSizes)))
   normFactor <- normFactor*libSizes/commLibSizes
   
   normCountsTable <- sweep(as.matrix(countsTable),2,normFactor,'/')
   
   convertTable <- readData(normCountsTable,factors = data.frame(grp=as.character(conds)))
   
 #  NOISeq.test <- noiseq(convertTable, factor = 'grp',replicates = 'technical',conditions = c('1','2'), k=0.5, norm='uqua')
   
   NOISeq.test <- noiseqbio(convertTable, norm = "uqua", factor = "grp", random.seed = 123, nclust = 15, lc = 0, r = 50,
                            adj = 1.5, a0per = 0.9, filter = 0)
   
  NOISeq.test@results[[1]]$prob <- 1 - NOISeq.test@results[[1]]$prob
   
   return(NOISeq.test@results[[1]])
#   NOISeq.Prob <- NOISeq.test@results[[1]]$prob
#   return(NOISeq.Prob)
 }


# DESeq2 - ----------------------------------------------------------------------------------#
 
run.DESeq2 <- function(countsTable,conds){
  require(DESeq2)
  
  grp <- factor(conds)
  
  DESeq.cds = DESeqDataSetFromMatrix(countData = countsTable, colData = DataFrame(grp), design = ~ grp)
  DESeq.cds = DESeq(DESeq.cds)
  DESeq.res = DESeq2::results(DESeq.cds)
  
  return(DESeq.res)
  
}


