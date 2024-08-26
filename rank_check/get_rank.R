library("R.matlab")
#setwd("/Users/yuhaow/Downloads/simulation") #feel free to set your own directory
quartly <- readMat("quarterly.mat")

#
#list of surrogates for S1, S2 and S3. The default option is listed below, which sets
#S1 as 1-2 quarters, S2 as 3-4 quarters, S3 as 5-6 quarters. To generate the results in
#Tables 1 and 3-4 in the Supplementary Material, you may change the list below according 
#to the choice of quarters prescribed in each cell of those tables. For example, to generate 
#data in the second row of Table 1 in the supplement, you need to set 
#surro1list <- c("ptcedd1", "ptcedd2", "ptcedd3", "ptcedd4"); 
#surro2list <- c("ptcedd5", "ptcedd6");
#surro3list <- c("ptcedd7", "ptcedd8", "ptcedd9", "ptcedd10")
#
surro1list <- c("ptcedd1", "ptcedd2") #list for the 1st surrogate
surro2list <- c("ptcedd3", "ptcedd4") #list for the 2nd surrogate
surro3list <- c("ptcedd5", "ptcedd6") #list for the 3rd surrogate

U <- cbind(quartly$grde911, quartly$grade12, quartly$grd1315, quartly$grade16, quartly$grd1720)
U <- U %*% c(1,2,3,3,3)
U <- U[quartly$county == 2]

binary.to.digit <- function(binvec){
	sum(sapply(1:length(binvec), function(idx) binvec[idx] * 2^(idx - 1)))
}

to.surro <- function(namelist){
  A <- matrix(0, nrow =length(quartly$county), ncol=length(namelist))
  for(idx in 1:length(namelist)){
    A[, idx] <- quartly[[namelist[idx]]]
  }
  A <- A[quartly$county == 2, ]

  B <- apply(A, 1, binary.to.digit)

  return(B)
}

get_indep_svd <- function(a, b, nsim){
	svs <- rep(0, nsim)
	for(i in 1:nsim){
		bs <- sample(b, size=length(b), replace=T)
		table <- table(a, bs)
		for(idx in 1:4) table[idx, ] <- table[idx,] / sum(table[idx,])
		svs[i] <- min(svd(table)$d)
	}
	return(svs)
}

S1 <- to.surro(surro1list)
S2 <- to.surro(surro2list)
S3 <- to.surro(surro3list)
Act <- quartly$e[quartly$county == 2]

for(act in 0:1){

print("SVD value for S1")
for(s2choice in 0:(2^length(surro2list) - 1)){
	print(paste("S2: ", s2choice))
	print(paste("A: ", act))
	tabl <- as.matrix(table(U[S2 == s2choice & Act == act], S1[S2 == s2choice & Act == act]))
	for(idx in 1:nrow(tabl)) tabl[idx, ] <- tabl[idx,] / sum(tabl[idx,])
	print(min(svd(tabl)$d))
	#print(table)
}

print("SVD value for S3")
for(s2choice in 0:(2^length(surro2list) - 1)){
	print(paste("S2 value: ", s2choice))
	print(paste("A: ", act))
	tabl <- as.matrix(table(U[S2 == s2choice & Act == act], S3[S2 == s2choice & Act == act]))
	for(idx in 1:nrow(tabl)) tabl[idx, ] <- tabl[idx,] / sum(tabl[idx,])
	print(min(svd(tabl)$d))
	#print(table)
}}
