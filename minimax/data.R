library(MASS)

#set up dimensions
nsim <- 200

#
#Change the dimension of X, S, U as well as the level
#of nonlinearity q according to the setting specified
#in each cell of the table; for example, if you would 
#like to replicate the simulation in the bottom-left 
#cell, then you have to set dimX <- 5; dimS <- dimU <- 5
#and set q <- 2
#
dimX <- 10
dimS <- dimU <- 5
q <- 2

nobs <- 2000
nexp <- 2000

set.seed(1)

g <- function(lindat) {
	dat <- matrix(0, nrow=nrow(lindat), ncol=ncol(lindat))
	for(i in 1:length(dat)){
		t <- lindat[i]
		dat[i] <- sign(t) * abs(t)^q
	}
	dat
}

rnd_norm_coef <- function(dim1, dim2, scale, min=0, max=1){

	coef <- matrix(0, nrow=dim1, ncol=dim2)

	for(j in 1:dim2){
		vec <-  runif(dim1, min=min, max=max)
		vec <- scale * vec / sqrt(sum(vec^2))
		coef[,j] <- vec
	}
	coef
}

rnd_data <- function(n, p){
	mvrnorm(n, mu=rep(0, p), Sigma=diag(p))
}

res.list <- c()

for(i in 1:nsim){
	print(i)

	#set up for system parameters
	kappaU <- rnd_norm_coef(dimU, 1, sqrt(0.5))
	kappaX <- rnd_norm_coef(dimX, 1, sqrt(0.5))
	tau1 <- rnd_norm_coef(1, dimS, sqrt(0.5)); beta1 <- rnd_norm_coef(dimX, dimS, sqrt(0.5)); gamma1 <- rnd_norm_coef(dimU, dimS, sqrt(0.5))
	tau2 <- sqrt(0.5); alpha2 <- rnd_norm_coef(dimS, 1, sqrt(0.5)); beta2 <- rnd_norm_coef(dimX, 1, sqrt(0.5)); gamma2 <- rnd_norm_coef(dimU, 1, sqrt(0.5))
	tau3 <- rnd_norm_coef(1, dimS, sqrt(0.5)); alpha3 <- rnd_norm_coef(1, dimS, sqrt(0.5)); beta3 <- rnd_norm_coef(dimX, dimS, sqrt(0.5)); gamma3 <- rnd_norm_coef(dimU, dimS, sqrt(0.5))
	tauy <- sqrt(0.5); alphay <- rnd_norm_coef(dimS, 1, sqrt(0.5)); betay <- rnd_norm_coef(dimX, 1, sqrt(0.5)); gammay <- rnd_norm_coef(dimU, 1, sqrt(0.5))

	#set up for observational data
	X <- rnd_data(nobs, dimX)
	U <- rnd_data(nobs, dimU)
	prob <- 1 / (1 + exp(X %*% kappaX + U %*% kappaU))
	A <- matrix(runif(nobs) <= prob, ncol=1)
	S1 <- A %*% tau1 + X %*% beta1 + U %*% gamma1 + sqrt(0.5) * rnd_data(nobs, dimS)
	S2 <- A * tau2 + S1 %*% alpha2 + X %*% beta2 + U %*% gamma2 + sqrt(0.5) * rnd_data(nobs, 1)
	S3 <- A %*% tau3 + S2 %*% alpha3 + X %*% beta3 + U %*% gamma3 + sqrt(0.5) * rnd_data(nobs, dimS)
	Y <- A * tauy + S3 %*% alphay + X %*% betay + U %*% gammay + sqrt(0.5) * rnd_data(nobs, 1)

	write.table(cbind(g(X), g(S2), g(S1), g(S3), Y, A), paste("tmp/obs_", i, ".csv", sep=""), row.names=FALSE, col.names=FALSE, sep=",")
	
	#set up for experimental data
	X <- rnd_data(nexp, dimX)
	U <- rnd_data(nexp, dimU)
	A <- matrix(runif(nexp) <= 1/2, ncol=1)
	S1 <- A %*% tau1 + X %*% beta1 + U %*% gamma1 + sqrt(0.5) * rnd_data(nexp, dimS)
	S2 <- A * tau2 + S1 %*% alpha2 + X %*% beta2 + U %*% gamma2 + sqrt(0.5) * rnd_data(nexp, 1)
	S3 <- A %*% tau3 + S2 %*% alpha3 + X %*% beta3 + U %*% gamma3 + sqrt(0.5) * rnd_data(nexp, dimS)
	Y <- A * tauy + S3 %*% alphay + X %*% betay + U %*% gammay + sqrt(0.5) * rnd_data(nexp, 1)

	write.table(cbind(g(X), g(S2), g(S1), g(S3), Y, A), paste("tmp/exp_", i, ".csv", sep=""), row.names=FALSE, col.names=FALSE, sep=",")
	S1 <- tau1; S2 <- tau2 + S1 %*% alpha2; S3 <- tau3 + S2 %*% alpha3
	Y <- tauy + S3 %*% alphay

	res.list <- c(res.list, Y)
}

save(res.list, file="tmp/result.rda")
