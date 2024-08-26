#load("data.rda")

#all the inputs are n times dim data matrix
outcome.bridge <- function(Y, S2, S1, X, S3){
  adj <- cbind(S1, S2, X)
  adj <- scale(adj, scale=F)
  cond <- cbind(S3, S2, X)
  condmean <- colMeans(cond)
  cond <- scale(cond, scale=F)
  vY <- t(adj) %*% Y
  cov <- t(adj) %*% cond
  b <- solve(cov, vY)
  D <- function(X){X <- X - condmean; return(sum(X * b) + mean(Y))}
  return(D)
}

# selection.bridge <- function(obsS2, obsS1, obsX, obsS3, expS2, expS1, expX, expS3){
#   adj <- rbind(cbind(obsS3, obsS2, obsX), cbind(expS3, expS2, expX))
#   adj <- scale(adj)
#   cond <- rbind(cbind(obsS1, obsS2, obsX), cbind(expS1, expS2, expX))
#   condmean <- colMeans(cond)
#   cond <- scale(cond)
#   obsindicator <- c(rep(1, dim(obsS1)[1]), rep(0, dim(expS1)[1]))
#   vY <- t(adj) %*% (!obsindicator) / mean(!obsindicator)
#   cov <- t(adj) %*% (obsindicator * cond) / mean(obsindicator)
#   b <- solve(cov, vY)
#   D <- function(X){X <- X - condmean; return(sum(X * b) + 1)}
#   return(D)
# }

# selection.bridge <- function(obsS2, obsS1, obsX, obsS3, expS2, expS1, expX, expS3){
#   adj <- rbind(cbind(obsS3, obsS2, obsX), cbind(expS3, expS2, expX))
#   adj <- scale(adj, scale=F)
#   cond <- cbind(obsS1, obsS2, obsX)
#   condmean <- colMeans(cond)
#   cond <- scale(cond, scale=F)
#   obsindicator <- c(rep(1, dim(obsS1)[1]), rep(0, dim(expS1)[1]))
#   vY <- t(adj) %*% (!obsindicator) / mean(!obsindicator)
#   cov <- t(adj[1:(dim(obsS1)[1]), ]) %*% cond / mean(obsindicator)
#   b <- solve(cov, vY)
#   D <- function(X){X <- X - condmean; return(sum(X * b) + 1)}
#   return(D)
# }
selection.bridge <- function(obsS2, obsS1, obsX, obsS3, expS2, expS1, expX, expS3){
  #determine the index of observational data
  obsindicator <- c(rep(TRUE, dim(obsS1)[1]), rep(FALSE, dim(expS1)[1]))
  
  #get initial parameter
  W <- rbind(cbind(obsS3, obsS2, obsX), cbind(expS3, expS2, expX))
  W[!obsindicator,] <- t(t(W[!obsindicator,]) - colMeans(W[obsindicator,]))
  W[obsindicator,] <- scale(W[obsindicator,], scale = FALSE)
  Z <- cbind(obsS1, obsS2, obsX)
  Zmean <- colMeans(Z)
  Z <- scale(Z, scale=F)
  
  #get more beyond the initial parameter
  ZZ <- cbind(1, Z)
  gmmobj <- function(theta){
    eta <- exp(ZZ %*% theta)
    eta <- sapply(eta, function(t) min(t, 10))
    avg.moment <- c(t(W[obsindicator,]) %*% eta) / sum(obsindicator) - colMeans(W[!obsindicator,])
    avg.moment <- c(avg.moment, mean(eta) - 1)
    sum(avg.moment^2)
  }
  theta <- optim(rep(0, ncol(ZZ)), gmmobj, control=list(reltol=1e-19, maxit=20000), method="L-BFGS-B")$par
  
  D <- function(X){X <- X - Zmean; X <- c(1, X); return(min(exp(sum(X * theta)), 10))}
  return(D)
}

#the one form Athey et al.
athey <- function(Y, S2, S1, X, S3){
  adj <- cbind(S3, S2, S1, X)
  adjmean <- colMeans(adj)
  adj <- scale(adj, scale=F)
  vY <- t(adj) %*% Y
  cov <- t(adj) %*% adj
  b <- solve(cov, vY)
  D <- function(X){X <- X - adjmean; return(sum(X * b) + mean(Y))}
  return(D)
}

res.list <- list()
nuisance.list <- list()

for(idx in 1:length(data.list)){
  print(idx)
  data <- data.list[[idx]]
  list2env(data, .GlobalEnv)
  
  #if you would like to generate Table 9 in the Supplementary Material, you need to uncomment the following 
  #two lines of codes
  # obsX <- matrix(nrow=nrow(obsX), ncol=0)
  # expX <- matrix(nrow=nrow(expX), ncol=0)
  
  #estimate the outcome bridge function with A == 1
  h1 <- outcome.bridge(obsY[obsA == 1], obsS2[obsA == 1,], obsS1[obsA == 1,], obsX[obsA == 1,], obsS3[obsA == 1,])
  #estimate the outcome bridge function with A == 0
  h0 <- outcome.bridge(obsY[obsA == 0], obsS2[obsA == 0,], obsS1[obsA == 0,], obsX[obsA == 0,], obsS3[obsA == 0,])
  tauest <- mean(apply(cbind(expS3[expA == 1,], expS2[expA == 1,], expX[expA == 1,]), 1, h1))
  tauest <- tauest - mean(apply(cbind(expS3[expA == 0,], expS2[expA == 0,], expX[expA == 0,]), 1, h0))

  #selection bridge function
  q1 <- selection.bridge(obsS2[obsA == 1,], obsS1[obsA == 1,], obsX[obsA == 1,], obsS3[obsA == 1,],
                         expS2[expA == 1,], expS1[expA == 1,], expX[expA == 1,], expS3[expA == 1,])
  q0 <- selection.bridge(obsS2[obsA == 0,], obsS1[obsA == 0,], obsX[obsA == 0,], obsS3[obsA == 0,],
                         expS2[expA == 0,], expS1[expA == 0,], expX[expA == 0,], expS3[expA == 0,])
  obsq1 <- apply(cbind(obsS1[obsA == 1,], obsS2[obsA == 1,], obsX[obsA == 1,]), 1, q1)
  tausel <- mean(obsq1 * obsY[obsA == 1])
  obsq0 <- apply(cbind(obsS1[obsA == 0,], obsS2[obsA == 0,], obsX[obsA == 0,]), 1, q0)
  tausel <- tausel - mean(obsq0 * obsY[obsA == 0])
  
  
  obsh1 <- apply(cbind(obsS3[obsA == 1,], obsS2[obsA == 1,], obsX[obsA == 1,]), 1, h1)
  obsq1 <- apply(cbind(obsS1[obsA == 1,], obsS2[obsA == 1,], obsX[obsA == 1,]), 1, q1)
  taudr <- tauest + mean(obsq1 * (obsY[obsA == 1] - obsh1))
  obsh0 <- apply(cbind(obsS3[obsA == 0,], obsS2[obsA == 0,], obsX[obsA == 0,]), 1, h0)
  obsq0 <- apply(cbind(obsS1[obsA == 0,], obsS2[obsA == 0,], obsX[obsA == 0,]), 1, q0)
  taudr <- taudr - mean(obsq0 * (obsY[obsA == 0] - obsh0))

  #the one from Athey et al.
  h1 <- athey(obsY[obsA == 1], obsS2[obsA == 1,], obsS1[obsA == 1,], obsX[obsA == 1,], obsS3[obsA == 1,])
  h0 <- athey(obsY[obsA == 0], obsS2[obsA == 0,], obsS1[obsA == 0,], obsX[obsA == 0,], obsS3[obsA == 0,])
  tauathey <- mean(apply(cbind(expS3[expA == 1,], expS2[expA == 1,], expS1[expA == 1,], expX[expA == 1,]), 1, h1))
  tauathey <- tauathey - mean(apply(cbind(expS3[expA == 0,], expS2[expA == 0,], expS1[expA == 0,], expX[expA == 0,]), 1, h0))

  #naive estimation
  taunaive <- mean(obsY[obsA == 1]) - mean(obsY[obsA == 0])
  
  #ground truth
  #tautrue <- 0.3 * (mean(expY[expA == 1]) - mean(expY[expA == 0])) + 0.7 * truetau
  tautrue <- mean(expY[expA == 1]) - mean(expY[expA == 0])
  
  # res.list[[idx]] <- mget(c("tauest", "tauathey", "taunaive", "tautrue", "truetau", "taudr"))
  res.list[[idx]] <- mget(c("tauathey", "taunaive", "tautrue", "truetau", "tauest", "taudr", "tausel"))
  nuisance.list[[idx]] <- mget(c("obsh1", "obsq1", "obsh0", "obsq0"))
}
