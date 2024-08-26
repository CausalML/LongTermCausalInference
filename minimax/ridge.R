set.seed(1)

library("glmnet")

outcome.ridge.bridge <- function(Y, S1, X, S3, lambda){
  adj <- cbind(S1, X)
  adj <- scale(adj, scale=F)
  cond <- cbind(S3, X)
  condmean <- colMeans(cond)
  cond <- scale(cond, scale=F)
  vY <- t(adj) %*% Y / length(Y)
  cov <- t(adj) %*% cond / length(Y)
  b <- solve(t(cov) %*% cov + lambda * diag(1 / (length(Y)), nrow=length(vY)), t(cov) %*% vY)
  D <- function(X){X <- X - condmean; return(sum(X * b) + mean(Y))}
  return(D)
}

selection.ridge.bridge <- function(obsS1, obsX, obsS3, expS1, expX, expS3, lambda){
  #determine the index of observational data
  obsindicator <- c(rep(TRUE, dim(obsS1)[1]), rep(FALSE, dim(expS1)[1]))
  
  #get initial parameter
  W <- rbind(cbind(obsS3, obsX), cbind(expS3, expX))
  W[!obsindicator,] <- t(t(W[!obsindicator,]) - colMeans(W[obsindicator,]))
  W[obsindicator,] <- scale(W[obsindicator,], scale = FALSE)
  Z <- cbind(obsS1, obsX)
  Zmean <- colMeans(Z)
  Z <- scale(Z, scale=F)
  
  #get more beyond the initial parameter
  ZZ <- cbind(1, Z)
  gmmobj <- function(theta){
    eta <- exp(ZZ %*% theta)
    eta <- sapply(eta, function(t) min(t, 10))
    avg.moment <- c(t(W[obsindicator,]) %*% eta) / sum(obsindicator) - colMeans(W[!obsindicator,])
    avg.moment <- c(avg.moment, mean(eta) - 1)
    sum(avg.moment^2) + lambda * 1 / (sum(obsindicator)) * sum(theta[-1]^2)
  }
  expcoef <- optim(rep(0, ncol(ZZ)), gmmobj, control=list(reltol=1e-19, maxit=20000), method="L-BFGS-B")$par
  
  D <- function(X){X <- X - Zmean; X <- c(1, X); return(min(exp(sum(X * expcoef)), 10))}
  return(D)
}

#
#Change the dimension of X and S according to the setting 
#specified in each cell of the table; for example, if you 
#would like to replicate the simulation in the bottom-left 
#cell of Table 1 in the main text, then you have to set 
#dimX = 5 and dimS = 5.
#
dimX = 10
dimS = 5
nsim = 200

dimX = dimX + 1 #here "dimX" refers to the dimension of X plus the dimension of S2

for(idx in 1:nsim){
  #the order of random variables in the csv file is: X, S2, S1, S3, Y, A
  print(idx)
  
  ###get experimental and observational data
  data = read.csv(paste("tmp/obs_", idx, ".csv", sep=""), header=F)
  obsX = data[, 1:dimX]##this includes both X and S2
  obsS1 = data[, (dimX+1):(dimX + dimS)]
  obsS3 = data[, (dimX + dimS+1):(dimX + 2 * dimS)]
  obsY = data[, dimX + 2 * dimS + 1]
  obsA = data[, dimX + 2 * dimS + 2]
  
  data = read.csv(paste("tmp/exp_", idx, ".csv", sep=""), header=F)
  expX = data[, 1:dimX]##this includes both X and S2
  expS1 = data[, (dimX+1):(dimX + dimS)]
  expS3 = data[, (dimX + dimS+1):(dimX + 2 * dimS)]
  expY = data[, dimX + 2 * dimS + 1]
  expA = data[, dimX + 2 * dimS + 2]
  
  #estimate the outcome bridge function with A == 1
  h1 <- outcome.ridge.bridge(obsY[obsA == 1], obsS1[obsA == 1,], obsX[obsA == 1,], obsS3[obsA == 1,], 1)
  exph1 <- apply(cbind(expS3[expA == 1,], expX[expA == 1,]), 1, h1)
  ate_or <- mean(exph1)
  
  #selection bridge function
  q1 <- selection.ridge.bridge(obsS1[obsA == 1,], obsX[obsA == 1,], obsS3[obsA == 1,],
                               expS1[expA == 1,], expX[expA == 1,], expS3[expA == 1,], 1)
  
  #doubly robust estimation
  obsh1 <- apply(cbind(obsS3[obsA == 1,], obsX[obsA == 1,]), 1, h1)
  obsq1 <- apply(cbind(obsS1[obsA == 1,], obsX[obsA == 1,]), 1, q1)
  ate_dr <- ate_or + mean(obsq1 * (obsY[obsA == 1] - obsh1))
  
  sigma = mean((exph1 - ate_dr)^2) / sum(expA) + mean(obsq1^2 * (obsY[obsA == 1] - obsh1)^2) / sum(obsA)
  sd = sqrt(sigma)
  
  write.table(c(ate_dr, ate_or, sd), file=paste("tmp/result_ridge_", idx, ".csv", sep=""), row.names=F, col.names=F)
}
