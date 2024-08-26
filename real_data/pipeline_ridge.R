args <- commandArgs(trailingOnly=TRUE)
offset <- args[1]

load(paste("data_", offset, ".rda", sep=""))
#load(paste("result_", offset, ".rda", sep=""))

#all the inputs are n times dim data matrix

source("main_ridge_n1.R")
save(res.ridge.list, file=paste("result_ridge_", offset, ".rda", sep=""))
