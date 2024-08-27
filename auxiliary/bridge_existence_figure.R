library(tidyverse)
library(latex2exp)

rep_num = 10000
res1 = rep(NA, rep_num)
res2 = rep(NA, rep_num)

for (i in seq_along(res1)){
  alpha2 = matrix(rnorm(4), 2, 2)
  gamma1 = matrix(rnorm(4), 2, 2)
  gamma2 = matrix(rnorm(4), 2, 2)
  gamma3 = matrix(rnorm(4), 2, 2)
  sigma1 = 0.5
  sigma2 = 0.5
  I = diag(2)
  
  M = sigma1 * alpha2 %*% t(alpha2) + I
  temp = solve(M, alpha2)
  final = (I - sigma1 * t(alpha2) %*% temp) %*% gamma1 -
    sigma1 * t(temp) %*% gamma2
  
  res1[i] = min(svd(final)$d)
  res2[i] = min(svd(gamma3)$d)
}

range(res1)
hist(res1)
hist(res2)

data.frame(singular_value = res1) %>% ggplot(aes(x = singular_value)) +
  geom_histogram(bins = 45, fill = "white", color = "black") + 
  geom_vline(aes(xintercept = 0.1), color = "red") +  
  xlab(TeX(r"( Smallest Singular Value of $\gamma_3$)"))
data.frame(singular_value = res2) %>% ggplot(aes(x = singular_value)) +
  geom_histogram(bins = 45, fill = "white", color = "black") + 
  geom_vline(aes(xintercept = 0.1), color = "red") +  
  xlab(TeX(r"( Smallest Singular Value of $\lambda_4$)"))


