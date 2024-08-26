args <- commandArgs(trailingOnly=TRUE)
offset <- args[1]

load(paste("data_", offset, ".rda", sep=""))

#all the inputs are n times dim data matrix

source("main.R")
source("main_ridge_n3.R")
save(res.list, nuisance.list, res.ridge.list4, nuisance.ridge.list4, file=paste("result_", offset, ".rda", sep=""))
