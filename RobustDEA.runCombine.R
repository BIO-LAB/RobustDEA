
combine.robust <- function(pvals,weights){
  pvals <- as.matrix(pvals)
  len <- dim(pvals)[2]
  sum.pvals <- data.frame(rowSums(pvals*weights))
  colnames(sum.pvals) <- c('sum')
  return(sum.pvals)
}

combine.fisher <- function(pvals) {
  pvals[pvals == 0] <- 0.00001
  fisher.sums <- data.frame(do.call(rbind,apply(pvals,1,fisher.sum,
                                                  zero.sub=zero.sub,na.rm=na.rm)))
  rownames(fisher.sums) <- rownames(pvals)
  fisher.sums$p.value <- 1-pchisq(fisher.sums$S,df=2*fisher.sums$num.p)
  fisher.sums$p.adj <- bonferroni = p.adjust(fisher.sums$p.value,"bonferroni")
  return(fisher.sums)
}

fisher.sum <- function(p) {
  p[p==0] <- 0.00001
  if (na.rm)
    p <- p[!is.na(p)]
  S = -2*sum(log(p))
  res <- data.frame(S=S,num.p=length(p))
  return(res)
}

combine.stouffer <- function (pvals) 
{
    weight <- rep(1, k)
    z <- qnorm(pvals, lower.tail = FALSE)
    cp <- pnorm(sum(weight * z)/sqrt(sum(weight^2)), lower.tail = FALSE)
  return(cp)
}