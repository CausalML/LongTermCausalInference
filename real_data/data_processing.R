set.seed(2)

data.list <- list()
totite <- 1000

for(idx in 1:totite){
  #get the index of selected data
  X <- cbind(quartly$grde911, quartly$grade12, quartly$grd1315, quartly$grade16, quartly$grd1720)
  X <- X %*% c(1,2,3,3,3)
  idxs <- 1:length(X)
  obssel <- runif(length(X)) <= 1
  X <- X[quartly$county == 2 & obssel]
  T <- quartly$e[quartly$county == 2 & obssel]
  trueY <- quartly[[outcomelist[1]]][quartly$county == 2 & obssel]
  idxs <- idxs[quartly$county == 2 & obssel]
  truetau <- mean(trueY[T == 1]) - mean(trueY[T == 0])
  prob0 <- c(1, 1 - offset / 3, 1 - offset / 3 * 2, max(1 - offset, 0.2))
  prob1 <- 1 - (prob0 - min(prob0)) * sum(T == 0) / sum(T == 1)
  X <- X + 1
  prob <- sapply(1:length(X), function(idx) if(T[idx] == 1) prob1[X[idx]] else prob0[X[idx]])
  sel <- runif(length(X)) <= prob
  idxs <- idxs[sel]
  
  #get the observational data
  # obsX <- cbind(to.matrix(covlist, idxs), countycov[idxs])
  obsX <- to.matrix(covlist, idxs)
  obsA <- to.matrix(elist, idxs)
  obsY <- to.matrix(outcomelist, idxs)
  obsS1 <- to.matrix(surro1list, idxs)
  obsS2 <- to.matrix(surro2list, idxs)
  obsS3 <- to.matrix(surro3list, idxs)
  
  #get the experimental data
  expidxs <- which(quartly$county == 3 & (obssel))
  # expX <- cbind(to.matrix(covlist, expidxs), countycov[expidxs])
  expX <- to.matrix(covlist, expidxs)
  expA <- to.matrix(elist, expidxs)
  expY <- to.matrix(outcomelist, expidxs)
  expS1 <- to.matrix(surro1list, expidxs)
  expS2 <- to.matrix(surro2list, expidxs)
  expS3 <- to.matrix(surro3list, expidxs)
  
  data.list[[idx]] <- mget(c("obsX", "obsA", "obsY", "obsS1", "obsS2", "obsS3", 
                             "expX", "expA", "expY", "expS1", "expS2", "expS3", "truetau"))
}
