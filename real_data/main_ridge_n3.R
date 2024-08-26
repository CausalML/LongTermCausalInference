library("glmnet")

outcome.ridge.bridge <- function(Y, S2, S1, X, S3, lambda){
  adj <- cbind(S1, S2, X)
  adj <- scale(adj, scale=F)
  cond <- cbind(S3, S2, X)
  condmean <- colMeans(cond)
  cond <- scale(cond, scale=F)
  vY <- t(adj) %*% Y / length(Y)
  cov <- t(adj) %*% cond / length(Y)
  b <- solve(t(cov) %*% cov + lambda * diag(1 / (length(Y)), nrow=length(vY)), t(cov) %*% vY)
  D <- function(X){X <- X - condmean; return(sum(X * b) + mean(Y))}
  return(D)
}

# selection.ridge.bridge <- function(obsS2, obsS1, obsX, obsS3, expS2, expS1, expX, expS3, lambda){
#   adj <- rbind(cbind(obsS3, obsS2, obsX), cbind(expS3, expS2, expX))
#   adj <- scale(adj)
#   cond <- rbind(cbind(obsS1, obsS2, obsX), cbind(expS1, expS2, expX))
#   condmean <- colMeans(cond)
#   cond <- scale(cond)
#   obsindicator <- c(rep(1, dim(obsS1)[1]), rep(0, dim(expS1)[1]))
#   vY <- t(adj) %*% (!obsindicator) / mean(!obsindicator) / sum(obsindicator)
#   cov <- t(adj) %*% (obsindicator * cond) / mean(obsindicator) / sum(obsindicator)
#   b <- solve(t(cov) %*% cov + lambda * diag(1 / (sum(obsindicator)), nrow=length(vY)), t(cov) %*% vY)
#   D <- function(X){X <- X - condmean; return(sum(X * b) + 1)}
#   return(D)
# }

# selection.ridge.bridge <- function(obsS2, obsS1, obsX, obsS3, expS2, expS1, expX, expS3, lambda){
#   adj <- rbind(cbind(obsS3, obsS2, obsX), cbind(expS3, expS2, expX))
#   adj <- scale(adj, scale=F)
#   cond <- cbind(obsS1, obsS2, obsX)
#   condmean <- colMeans(cond)
#   cond <- scale(cond, scale=F)
#   obsindicator <- c(rep(1, dim(obsS1)[1]), rep(0, dim(expS1)[1]))
#   vY <- t(adj) %*% (!obsindicator) / mean(!obsindicator) / sum(obsindicator)
#   cov <- t(adj[1:(dim(obsS1)[1]), ]) %*% cond / mean(obsindicator) / sum(obsindicator)
#   b <- solve(t(cov) %*% cov + lambda * diag(1 / (sum(obsindicator)), nrow=length(vY)), t(cov) %*% vY)
#   D <- function(X){X <- X - condmean; return(sum(X * b) + 1)}
#   return(D)
# }
selection.ridge.bridge <- function(obsS2, obsS1, obsX, obsS3, expS2, expS1, expX, expS3, lambda){
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
    sum(avg.moment^2) + lambda * 1 / (sum(obsindicator)) * sum(theta[-1]^2)
  }
  expcoef <- optim(rep(0, ncol(ZZ)), gmmobj, control=list(reltol=1e-19, maxit=20000), method="L-BFGS-B")$par
  
  D <- function(X){X <- X - Zmean; X <- c(1, X); return(min(exp(sum(X * expcoef)), 10))}
  return(D)
}

athey.ridge <- function(Y, S2, S1, X, S3, lambda){
  adj <- cbind(S3, S2, S1, X)
  adjmean <- colMeans(adj)
  adj <- scale(adj, scale=F)
  vY <- t(adj) %*% Y / length(Y)
  cov <- t(adj) %*% adj / length(Y)
  b <- solve(cov + lambda * diag(1 / (length(Y)), nrow=length(vY)), vY)
  D <- function(X){X <- X - adjmean; return(sum(X * b) + mean(Y))}
  return(D)
}

athey.cv <- function(Y, S2, S1, X, S3){
  adj <- cbind(S3, S2, S1, X)
  cv.fit <- cv.glmnet(adj, Y, alpha=0, family="gaussian")
  D <- function(X){return(predict(cv.fit, newx=X))}
  return(D)
}

res.ridge.list4 <- list()
nuisance.ridge.list4 <- list()

for(idx in 1:length(data.list)){
  print(idx)
  data <- data.list[[idx]]
  list2env(data, .GlobalEnv)
  
  #if you would like to generate Table 9 in the Supplementary Material, you need to uncomment the following 
  #two lines of codes
  # obsX <- matrix(nrow=nrow(obsX), ncol=0)
  # expX <- matrix(nrow=nrow(expX), ncol=0)
  
  #estimate the outcome bridge function with A == 1
  h1 <- outcome.ridge.bridge(obsY[obsA == 1], obsS2[obsA == 1,], obsS1[obsA == 1,], obsX[obsA == 1,], obsS3[obsA == 1,], 1)
  #estimate the outcome bridge function with A == 0
  h0 <- outcome.ridge.bridge(obsY[obsA == 0], obsS2[obsA == 0,], obsS1[obsA == 0,], obsX[obsA == 0,], obsS3[obsA == 0,], 1)
  tauridge1 <- mean(apply(cbind(expS3[expA == 1,], expS2[expA == 1,], expX[expA == 1,]), 1, h1))
  tauridge1 <- tauridge1 - mean(apply(cbind(expS3[expA == 0,], expS2[expA == 0,], expX[expA == 0,]), 1, h0))
  
  #selection bridge function
  q1 <- selection.ridge.bridge(obsS2[obsA == 1,], obsS1[obsA == 1,], obsX[obsA == 1,], obsS3[obsA == 1,],
                               expS2[expA == 1,], expS1[expA == 1,], expX[expA == 1,], expS3[expA == 1,], 1)
  q0 <- selection.ridge.bridge(obsS2[obsA == 0,], obsS1[obsA == 0,], obsX[obsA == 0,], obsS3[obsA == 0,],
                               expS2[expA == 0,], expS1[expA == 0,], expX[expA == 0,], expS3[expA == 0,], 1)
  obsq1 <- apply(cbind(obsS1[obsA == 1,], obsS2[obsA == 1,], obsX[obsA == 1,]), 1, q1)
  tauselridge1 <- mean(obsq1 * obsY[obsA == 1])
  obsq0 <- apply(cbind(obsS1[obsA == 0,], obsS2[obsA == 0,], obsX[obsA == 0,]), 1, q0)
  tauselridge1 <- tauselridge1 - mean(obsq0 * obsY[obsA == 0])
  
  #doubly robust estimation
  obsh1 <- apply(cbind(obsS3[obsA == 1,], obsS2[obsA == 1,], obsX[obsA == 1,]), 1, h1)
  obsq1 <- apply(cbind(obsS1[obsA == 1,], obsS2[obsA == 1,], obsX[obsA == 1,]), 1, q1)
  taudrridge1 <- tauridge1 + mean(obsq1 * (obsY[obsA == 1] - obsh1))
  obsh0 <- apply(cbind(obsS3[obsA == 0,], obsS2[obsA == 0,], obsX[obsA == 0,]), 1, h0)
  obsq0 <- apply(cbind(obsS1[obsA == 0,], obsS2[obsA == 0,], obsX[obsA == 0,]), 1, q0)
  taudrridge1 <- taudrridge1 - mean(obsq0 * (obsY[obsA == 0] - obsh0))
  
  #estimate the imputation function by Athey et al. with ridge regularization
  h1 <- athey.ridge(obsY[obsA == 1], obsS2[obsA == 1,], obsS1[obsA == 1,], obsX[obsA == 1,], obsS3[obsA == 1,], 1)
  h0 <- athey.ridge(obsY[obsA == 0], obsS2[obsA == 0,], obsS1[obsA == 0,], obsX[obsA == 0,], obsS3[obsA == 0,], 1)
  tauatheyridge <- mean(apply(cbind(expS3[expA == 1,], expS2[expA == 1,], expS1[expA == 1,], expX[expA == 1,]), 1, h1))
  tauatheyridge <- tauatheyridge - mean(apply(cbind(expS3[expA == 0,], expS2[expA == 0,], expS1[expA == 0,], expX[expA == 0,]), 1, h0))
  
  #estimate by Athey et al. via cross validation
  h1 <- athey.cv(obsY[obsA == 1], obsS2[obsA == 1,], obsS1[obsA == 1,], obsX[obsA == 1,], obsS3[obsA == 1,])
  h0 <- athey.cv(obsY[obsA == 0], obsS2[obsA == 0,], obsS1[obsA == 0,], obsX[obsA == 0,], obsS3[obsA == 0,])
  tauatheycv <- mean(h1(cbind(expS3[expA == 1,], expS2[expA == 1,], expS1[expA == 1,], expX[expA == 1,])))
  tauatheycv <- tauatheycv - mean(h0(cbind(expS3[expA == 0,], expS2[expA == 0,], expS1[expA == 0,], expX[expA == 0,])))
  
  #ground truth
  #tautrue <- 0.3 * (mean(expY[expA == 1]) - mean(expY[expA == 0])) + 0.7 * truetau
  tautrue <- mean(expY[expA == 1]) - mean(expY[expA == 0])
  # res.list[[idx]] <- mget(c("tauest", "tauathey", "tauatheypar", "taunaive", "tautrue", "truetau"))
  # res.ridge.list4[[idx]] <- mget(c("tauridge0", "tauridge1", "tauridge2", "tauridge3", "tauridge4", "tauridge5",
  #                                  "taudrridge0", "taudrridge1", "taudrridge2", "taudrridge3", "taudrridge4", "taudrridge5",
  #                                 "tautrue", "truetau"))
  res.ridge.list4[[idx]] <- mget(c("tauridge1", "tauselridge1", "taudrridge1", "tauatheyridge", "tauatheycv", "tautrue", "truetau"))
  nuisance.ridge.list4[[idx]] <- mget(c("obsh1", "obsq1", "obsh0", "obsq0"))
}
