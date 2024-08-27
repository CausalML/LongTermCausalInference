library(MASS)

#set up dimensions
nsim <- 200
alpha <- 0.05

folder <- "tmp/"

drerr <- function(i, resname){
	t <- read.csv(paste(folder, resname, "_", i, ".csv", sep=""), header=F)$V1[1]
	load(paste(folder, "result.rda", sep=""))
	tt <- res.list[i]
	t - tt
}

cp <- function(i, resname){
	t <- read.csv(paste(folder, resname, "_", i, ".csv", sep=""), header=F)$V1
	load(paste(folder, "result.rda", sep=""))
	tt <- res.list[i]
	abs(t[1] - tt) <= qnorm(1 - alpha / 2) * t[3]
}

cilen <- function(i, resname){
	t <- read.csv(paste(folder, resname, "_", i, ".csv", sep=""), header=F)$V1
	load(paste(folder, "result.rda", sep=""))
	2 * qnorm(1 - alpha / 2) * t[3]
}

bias_ratio <- function(a, b){sqrt(a^2 / b)}

print(mean(sapply(1:nsim, function(idx) cp(idx, "result"))))

print(mean(sapply(1:nsim, function(idx )cilen(idx, "result"))))

print(sqrt(mean(sapply(1:nsim, function(idx) drerr(idx, "result"))^2)))

print(mean(sapply(1:nsim, function(idx) drerr(idx, "result"))))

#ridge
print(mean(sapply(1:nsim, function(idx) cp(idx, "result_ridge"))))

print(mean(sapply(1:nsim, function(idx) cilen(idx, "result_ridge"))))

print(sqrt(mean(sapply(1:nsim, function(idx) drerr(idx, "result_ridge"))^2)))

print(mean(sapply(1:nsim, function(idx) drerr(idx, "result_ridge"))))
