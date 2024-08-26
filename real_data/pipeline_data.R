library("R.matlab")
#setwd("/Users/yuhaow/Downloads/simulation") #feel free to set your own directory
quartly <- readMat("quarterly.mat")

#
#list of surrogates for S1, S2 and S3. The default option is listed below, which sets
#S1 as 1-2 quarters, S2 as 3-4 quarters, S3 as 5-6 quarters. To generate the results in
#Tables 5-8, you may change the list below according to the choice of quarters prescribed
#in each table. For example, to generate Table 5, you need to set 
#surro1list <- c("ptcedd1", "ptcedd2", "ptcedd3", "ptcedd4"); surro2list <- c("ptcedd5", "ptcedd6")
#surro3list <- c("ptcedd7", "ptcedd8", "ptcedd9", "ptcedd10")
#

surro1list <- c("ptcedd1", "ptcedd2") #list for the 1st surrogate
surro2list <- c("ptcedd3", "ptcedd4") #list for the 2nd surrogate
surro3list <- c("ptcedd5", "ptcedd6") #list for the 3rd surrogate
outcomelist <- c("ptcedd20")
elist <- c("e")

covlist <- c("grew1", "gepop1", "xsexf", "xhsdip", "x1chld",
             "single", "dumkids", "xchld05", "white", "hisp", "black", "age")

to.matrix <- function(namelist, idxs){
  A <- matrix(0, nrow =length(quartly$county), ncol=length(namelist))
  for(idx in 1:length(namelist)){
    A[, idx] <- quartly[[namelist[idx]]]
  }
  return(A[idxs, ])
}

offset.list <- list(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6)

for(offset in offset.list){
  
  source("data_processing.R")
  
  save(data.list, file=paste("data_", offset, ".rda", sep=""))
}
