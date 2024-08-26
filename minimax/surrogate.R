library(MASS)

#set up dimensions
args <- commandArgs(trailingOnly=TRUE)
start <- as.integer(args[1])
start <- start * 8
nsim <- 8

res.list <- c()

for(i in 1:nsim){
	print(i)

	#set up for system parameters
	system(paste("python3 surrogate.py", start + i, sep= " "))
}
