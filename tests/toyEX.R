# Simulation setup

if(!interactive()){
  args <- as.numeric(commandArgs(trailingOnly = TRUE))
}else{
  args <- c(10)
}

timenow1 = Sys.time()
timenow0 = gsub(' ', '_', gsub('[-:]', '', timenow1))
timenow = paste(timenow0, ".txt", sep = "")

# install.packages( "MatrixModels", type="win.binary" )  

library(np)
# library(ggplot2)
library(nleqslv)
library(knitr)
# suppressMessages(library(CVXR))
# suppressMessages(library(kableExtra))
suppressMessages(library(dplyr))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))

set.seed(10)
SIMNUM = args[1]

if(!interactive()){
  dir.create(timenow0)
  setwd(timenow0)
  
  sink(timenow, append=TRUE)
}
# setup parallel backend to use many processors
cores = min(detectCores() - 3, 101)
print(paste("cores =", cores))
# cl <- makeCluster(cores, outfile = timenow) #not to overload your computer
cl <- makeCluster(cores)
registerDoParallel(cl)

N = 10000
x1= rnorm(N, 2, 1)
x1 = ifelse(x1 < 0, 0, x1)
x1 = ifelse(x1 > 4, 4, x1)
# x = cbind(rnorm(N, 2, 1), runif(N, 0, 4), runif(N, 0, 4))
x = cbind(x1, runif(N, 0, 4))
e = rnorm(N, 0, 0.2)
# y = x[,1] + x[,2] + e
y = x[,1] + x[,2]^2 / 8 * 3 + e

pi = 1 / (1 + exp(-(0.2 * y - 1.5))) # Poisson
mean(pi)
summary(pi)

theta = mean(y)

delta = rbinom(N, 1, pi)
Index_S = (delta == 1)
pi_S = pi[Index_S]

data = data.frame(x)
data_S = cbind(pi, data)[Index_S,]
bw1 <- np::npregbw(reformulate(colnames(data_S[,-1,drop = F]), response = "1 / pi_S * (1 / pi_S - 1)"),
                   regtype = "lc", data = data_S[,-1,drop = F])

final_res <- foreach(
  simnum = 1:SIMNUM, 
  .packages = c("nleqslv", "np"), 
  .errorhandling="pass") %dopar% {
    # set.seed(simnum) # To be removed
    
    # Poisson sampling
    delta = rbinom(N, 1, pi)
    Index_S = (delta == 1)
    n = sum(Index_S); #print(n)
    pi_S = pi[Index_S]
    d_S0 = 1 / pi_S
    pimat_S = diag(d_S0^2 - d_S0) / N^2 # 1 / pi_i * (1 - 1 / pi_i)
    
    x_S = x[Index_S,,drop = F]
    y_S = y[Index_S] # plot(x_S, y_S)
    
    Z_S = cbind(1, x_S)
    Z_St = t(Z_S)
    Zbar = colMeans(cbind(1, x))
    
    gammahat = solve(Z_St %*% (Z_S * d_S0), Z_St %*% (y_S * d_S0))
    theta1 = drop(t(Zbar) %*% gammahat) + sum((y_S - drop(Z_S %*% gammahat)) * d_S0) / N
    
    yvec = npksum(txdat=drop(x_S), tydat= d_S0 * (d_S0 - 1) * y_S , exdat = drop(x), bws=bw1$bw)$ksum/
      npksum(txdat=drop(x_S), tydat= d_S0 * (d_S0 - 1), exdat = drop(x), bws=bw1$bw)$ksum
    
    theta2 = mean(yvec) + sum((y_S  - yvec[Index_S]) / pi_S) / sum(1 / pi_S)
    
    V1 = crossprod(y_S - drop(Z_S %*% gammahat), pimat_S %*% (y_S - drop(Z_S %*% gammahat)))
    
    V2 = crossprod(y_S - yvec[Index_S], pimat_S %*% (y_S - yvec[Index_S]))
    
    theta_vec = c(theta1, theta2) - theta
    var_vec = c(V1, V2)
    
    CR_vec = ifelse(abs(theta_vec) < qnorm(0.975) * sqrt(var_vec), 1, 0)
    
    list(theta = theta_vec, 
         Var = var_vec,
         CR = CR_vec)
  }

final_res1 = lapply(final_res, function(x) x[[1]])
final_res2 = lapply(final_res, function(x) x[[2]])
final_res3 = lapply(final_res, function(x) x[[3]])

stopCluster(cl)
timenow2 = Sys.time()
print(timenow2 - timenow1)
# if(sum(!sapply(final_res1, function(x) is.numeric(unlist(x)))) != 0) stop(paste(final_res1))
paste("# of failure:", sum(!sapply(final_res1, function(x) is.numeric(unlist(x))))); final_res0 = final_res1
final_res0[!sapply(final_res0, function(x) is.numeric(unlist(x)))]
final_res1 = final_res1[sapply(final_res1, function(x) is.numeric(unlist(x)))]
res1 = do.call("rbind", final_res1)

print(res1[apply(res1, 1, function(x) any(is.na(x) | is.infinite(x))),])
colnames(res1)

Bias = colMeans(res1, na.rm = T)
SE = apply(res1, 2, function(x) sqrt(var(x, na.rm = T) * (sum(!is.na(x))-1)/sum(!is.na(x)) ))
RMSE = apply(res1, 2, function(x) sqrt(mean(x^2, na.rm = T)))

res2 = do.call("rbind", final_res2)
# colMeans(sqrt(res2), na.rm = T)
RB = (colMeans(res2, na.rm = T) - SE^2) / SE^2

res3 = do.call("rbind", final_res3)
CR = colMeans(sqrt(res3), na.rm = T)

round(cbind(Bias, SE, RMSE, RB, CR), 4)
# round(cbind(Bias, SE, RMSE, RB, CR), 6)

xtable::xtable(cbind(Bias, SE, RMSE, RB, CR), digits = c(1,4,4,4,2,3))
xtable::xtable(cbind(cbind(Bias, SE, RMSE) * 1e2, RB, CR), digits = c(1,3,3,3,2,3))
