#offset.list <- list(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
# offset.list <- list(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4)
#offset.list <- list(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4)
offset.list <- list(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6)
sc <- 1
digi <- 3

trim <- function(est){
  if(est <= 0)
    0
  else{
    if(est > 1)
      1
    else
      est
  }
}

getmean <- function(listval, keyword){
  mean(sapply(listval, function(res) {abs(res[[keyword]] - res$tautrue)}))
}

getmed <- function(listval, keyword){
  median(sapply(listval, function(res) {abs(res[[keyword]] - res$tautrue)}))
}

formatting <- function(vec, digi){
	sprintf(round(vec, digi), fmt = paste("%#.", digi, "f", sep=""))
}

formatperc <- function(mtd, gnd, digi){
	as.character(round((1 - mtd / gnd) * 100))
}

totpair <- 15
totoffset <- length(offset.list)

resmaetab <- matrix(nrow=totpair, ncol=totoffset)
resmedtab <- matrix(nrow=totpair, ncol=totoffset)
resdigimaetab <- matrix(nrow=totpair, ncol=totoffset)
resdigimedtab <- matrix(nrow=totpair, ncol=totoffset)

for(oidx in 1:length(offset.list)){
  offset <- offset.list[oidx]

  load(paste("result_", offset, ".rda", sep=""))
  load(paste("result_ridge_", offset, ".rda", sep=""))
  print(paste("offset:", offset))
  
  for(idx in 1:length(res.ridge.list)) {res.ridge.list[[idx]]$tautrue <- res.list[[idx]]$tautrue}
  
  pairs <- list(list(res.list, "tauest"), list(res.ridge.list, "tauridge1"), list(res.ridge.list, "tauridge2"),
                list(res.ridge.list4, "tauridge1"), list(res.list, "tausel"), list(res.ridge.list, "tauselridge1"),
                list(res.ridge.list, "tauselridge2"), list(res.ridge.list4, "tauselridge1"), list(res.list, "taudr"),
                list(res.ridge.list, "taudrridge1"), list(res.ridge.list, "taudrridge2"), 
                list(res.ridge.list4, "taudrridge1"), list(res.list, "tauathey"),
                list(res.ridge.list4, "tauatheycv"), list(res.list, "taunaive"))
  
  # for(i in 1:length(res.list)){
  #   res.list[[i]]$tauest <- trim(res.list[[i]]$tauest)
  #   res.list[[i]]$tausel <- trim(res.list[[i]]$tausel)
  #   res.list[[i]]$taudr <- trim(res.list[[i]]$taudr)
  # }

  for(idx in 1:totpair) {resdigimaetab[idx, oidx] <- getmean(pairs[[idx]][[1]], pairs[[idx]][[2]])}
  
  print(resdigimaetab[1, oidx])
  
  for(idx in 1:totpair) {resdigimedtab[idx, oidx] <- getmed(pairs[[idx]][[1]], pairs[[idx]][[2]])}

  for(idx in 1:(totpair-1)) 
    {resmaetab[idx, oidx] <- formatperc(resdigimaetab[idx, oidx], resdigimaetab[totpair, oidx], 0)}
  resmaetab[totpair, oidx] <- formatting(resdigimaetab[totpair, oidx], 3)
  
  for(idx in 1:(totpair-1)) 
    {resmedtab[idx, oidx] <- formatperc(resdigimedtab[idx, oidx], resdigimedtab[totpair, oidx], 0)}
  resmedtab[totpair, oidx] <- formatting(resdigimedtab[totpair, oidx], 3)
  
  
  print(paste(resmaetab[,oidx], collapse=" & "))
  print(paste(resmedtab[,oidx], collapse=" & "))
}
