# Simulation setup

if(!interactive()){
  args <- as.numeric(commandArgs(trailingOnly = TRUE))
}else{
  args <- c(5)
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

set.seed(11)
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

# N = 100000
# x = matrix(rnorm(N, 1, 1), nc = 1)
# e = rnorm(N, 0, 1)
# beta = c(1, 1)
# y = beta[1] + x[,1] * beta[2] + e
# # y = beta[1] + x[,1]^2 * beta[2] + e
# # y = (x[,1] - 1)^2 + e
# # y = beta[1] + 1 / (x[,1]^2 + 1)* beta[2] + e
# 
# En = 200
# pi = 1 / ((x + y + 10) / 100)
# pi = pi / sum(pi) * En

# pi = 1 / (1 + exp(-(-6 + .2 * x[,1]^2)))
# pi = 1 / (1 + exp(-(-4 + .2 * x[,1]^2)))
# pi = 1 / (1 + exp(-(-6 -  0.2 * x[,1]^2)))
# pi = 1 / (1 + exp(-(-6 -  0.2 * x[,1]^2)))
# pi = 1 / (1 + exp(-(rep(-6, N))))
# pi = 1 / (1 + exp(-(-7 + .2 * x[,1] + 0.2 * y)))
# pi = 1 / (1 + exp(-(-6 + 0.2 * e)))
# pi = 1 / (1 + exp(-(-0.2 * e)))
# pi = 1 / (1 + exp(-(-6 + 0.2 * sin(e))))
# pi[order(pi)[1:100]] <- pi[order(pi)[1:100]] / 100
# summary(pi)

#### Simulation setup from a note on weight smoothing in survey sampling; kim2023note
#### They used poisson sampling
# N = 10000
# En = 500
# x = matrix(runif(N, 0, 2), nc = 1)
# e = rnorm(N, 0, 0.5)
# y = 0.5 + 0.5 * x[,1] + e
# 
# mpi = 1 / (1 + exp(-(x + y)))
# mpi = mpi / sum(mpi) * En
# phi = 100
# pi = rbeta(N, mpi * phi, (1 - mpi) * phi)

#### Simulation setup from kim2010exponential
#### They used PPSWR sampling.
# N = 10000
# En = 500
# x = matrix(rexp(N, 1) + 1, nc = 1)
# e = rnorm(N, 0, 1)
# # y = 3 + x + x * e
# y = (5 - 1 / sqrt(8)) + 1 / sqrt(8) * (x - 2)^2 + e
# z = rchisq(N, 1) + abs(y)
# # z = rchisq(N, 1)
# prob = z / sum(z)
# pi = 1 - (1 - prob)^En

# #### Simulation setup from Wu & Rao 2009; wu2009pseudo
# #### They used PPSWOR sampling; We use Case (ii)
# N = 2000
# En = 200
# x = matrix(rexp(N), nc = 1)
# z = matrix(rexp(N), nc = 1)
# e = rchisq(N, 1) - 1
# y = 1 + z + x + 2 * e
# # cor(y, 1 + z + x)
# pi = (z / sum(z)) * En
# # max(pi)
# # prob = z / sum(z)
# # pi = 1 - (1 - prob)^En

#### Simulation setup from Qin et. al. 2002; qin2002estimation
#### They used Poisson sampling (under missing data setup)
# N = 3000
# x = matrix(rchisq(N, 6), nc = 1) / 2
# # x = matrix(rchisq(N * 2, 6), nc = 2) / 2
# 
# # e = rnorm(N, 0, 1); # e = rnorm(N)
# e = rnorm(N, 0, 0.5);
# # y = x + e * sqrt(x)                    # Model 1
# # y = x + 0.05 * x^2 + e * sqrt(x)       # Model 2
# # y = 1.5 + x + e * sqrt(x)              # Model 3
# # y = 3 + x - 0.05 * x^2 + e * sqrt(x)   # Model 4
# # y = -2 + x[,1] + 0.2 * x[,1]^2 + e * sqrt(x[,1])   # Works well
# # y = -0.5 + x + 0.1 * x^2 + e * sqrt(x)   # Works well
# y = 1.5 + x + e
# #### a_qin = -3, -2, -1, 1
# a_qin = -2
# z = 0.2 * y + a_qin
# 
# pi = 1 / (1 + exp(-z)) # Poisson
# mean(pi)
# summary(pi)

print("Sum g(d_i) is unknown")
#### New simulation setup with two covariates
N = 10000
# N = 2000
x1= rnorm(N, 2, 1)
# x = cbind(x1, runif(N, 0, 4), runif(N, 0, 4))
x = cbind(x1, runif(N, 0, 4))
e = rnorm(N, 0, 1)
# e = rnorm(N, 0, 0.2)
print(paste("sd(e) =", sd(e)))
# z = runif(N, -2, 2)
z = rnorm(N, 0, 1)
# y = x[,1] + x[,2] + z + e; print("linear, z")
# y = x[,1] + x[,2]^2 / 8 * 3 + z + e; print("nonlinear, z")

# y = x[,1] + x[,2] + z^2 + e; print("linear, z^2")
# y = x[,1] + x[,2]^2 / 8 * 3 + z^2 + e; print("nonlinear, z^2")

y = x[,1] + x[,2]^2 / 8 * 3 + e; print("noninformative, nonlinear, z")

# pi = pt(-z / 1.5 - 2, 3) # Works well
# pi = pt(-z / 2 - 2, 3) # Works well
pi = pt(-z - 2, 3) # Works well
pi = ifelse(pi >.7, .7, pi)

mean(1 / pi) / sd(1 / pi); print("CV") # Coefficient of Variation

# pi = pt(-z / 5 - 2, 3) # Works well

# pi = pt(-z/5 - x[,1] * x[,2], 3)

# pi = pt(-z * x[,1] - 2, 3)

# pi = pt(- z + x[,1] - 4.5, 3)

# pi = pt(- z + x[,1] - x[,2] - 3, 3)

# pi = pt(- z + x[,1] - x[,2] + 1, 3)

# pi = pt(x1^2 / 5 - x[,2]- 1.5, 3)

# pi = pt(x1^2 / 5 + x[,2]- 5.5, 3)

# pi = pt(-2 * z - 2, 3)
# pi = ifelse(pi >.7, .7, pi)

# pi = pt(-3 * z - 2, 3)
# pi = ifelse(pi >.7, .7, pi)

# pi = ifelse(pi < 1e-3, runif(N, pi, 1e-3), pi)
# pi = ifelse(pi >1 -  1e-3, runif(N, 1 - 1e-3, pi), pi)

# pi = ifelse(pi < 1e-3, runif(N, ifelse(pi > 1e-4, pi, 1e-4), 1e-3), pi)
# pi = ifelse(pi >1 -  1e-3, runif(N, 1 - 1e-3, pi), pi)

# pi = ifelse(pi < 1e-3, runif(N, pi, 1e-3), pi)
# pi = ifelse(pi > 0.9, runif(N, 0, 0.9), pi)

# pi = ifelse(pi < 1e-1, runif(N, pi, 1e-1), pi)
# pi = ifelse(pi >1 -  1e-1, runif(N, 1 - 1e-1, pi), pi)

# pi = 1 / (1 + exp(-(0.2 * x1 - 1.5))) # Poisson
mean(pi)
summary(pi)
# hist(pi)

# if(type == "Hb") Hcnt = Hcnt + 1
#### used for Huber loss. ####
# summary(1 / pi)
# if(Hcnt == 1){
#   del = quantile(1 / pi, 0.25)
# }else if(Hcnt == 2){
#   del = quantile(1 / pi, 0.50)
# }else if(Hcnt == 3){
#   del = quantile(1 / pi, 0.75)        
# }
del = quantile(1 / pi, 0.80)  

# En = 350 # PPS
# prob = 1 / (1 + exp(-z))
# prob = prob / sum(prob)
# pi = 1 - (1 - prob)^En
# mean(pi)
# summary(pi)

# y = (y < 7)
# z = y + 5
# pi = z / sum(z) * 400
# mean(pi)


# We need two codntions to see outperformance of the proposed method:
# y is non-linear ftn
# pi and x are related in a functional form.

theta = mean(y)
# hist(pi)

# type_vec = c( "SL", "EL", "ET", "CE", "HD")
type_vec = c("EL", "ET", "CE", "HD", "PH")
# type_vec = c("PH")
# type_vec = c("HD")
# cal_vec = c("DS", "GEC1", "GEC2", "Knl2", "GEC0")
# cal_vec = c("DS", "GEC1", "GEC2")
# cal_vec = c("DS")

# Used only for bandwidth finding
# data = data.frame(x)
# bw1 <- np::npregbw(reformulate(colnames(data), response = "1 / pi"),
#                    regtype = "lc", data = data)

# delta = rbinom(N, 1, pi)
# Index_S = (delta == 1)
# pi_S = pi[Index_S]
# data = data.frame(x)
# data_S = cbind(pi, data)[Index_S,]
# bw1 <- np::npregbw(reformulate(colnames(data_S[,-1,drop = F]), response = "1 / pi_S"),
#                    regtype = "lc", data = data_S[,-1,drop = F])

# plot(bw1)
# points(x[Index_S], 1 / pi_S * (1 / pi_S - 1), col = "red")

# bw0 <- np::npregbw(reformulate(colnames(data), response = "1 / pi"),
#                    regtype = "lc", data = data)

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
    
    # PPS-WR sampling
    # Index_S = unique(sample(1:N, size = En, replace = T, prob = prob))
    # n = length(Index_S)
    # pi_S = pi[Index_S]
    # d_S0 = 1 / pi_S
    # prob_S = prob[Index_S]
    # pimat_S = outer(pi_S, pi_S, "+") - (1 - (1 - outer(prob_S, prob_S, "+"))^En)
    # diag(pimat_S) = pi_S  # pi_{ij}
    # pimat_S = (outer(d_S0, d_S0, "*") - 1 / pimat_S) / N^2
    
    x_S = x[Index_S,,drop = F]
    y_S = y[Index_S] # plot(x_S, y_S)
    data = data.frame(x)
    data_S = cbind(pi, data)[Index_S,]
    
    #Hajek estimator
    theta_Hajek = sum(y_S * d_S0) / sum(d_S0)
    theta_HT = sum(y_S * d_S0) / N
    
    Var_HT = crossprod(y_S, pimat_S %*% y_S)
    Var_Hajek = crossprod(y_S - theta_Hajek, pimat_S %*% (y_S - theta_Hajek)) 
    
    G = function(x, type, del){
      switch(type,
             SL = x^2/2,
             EL = -log(x),
             ET = x * (log(x) - 1),
             CE = (x-1) * log(x-1) - x * log(x),
             HD = -2 * sqrt(x),
             PH = del^2 * sqrt(1 + (x / del)^2))
    }
    
    # g = function(x, type){
    #   switch(type,
    #          SL = x,
    #          EL = -1 / x,
    #          ET = log(x),
    #          CE = log(1 - 1 / x))
    # }
    
    # ginv = function(x, type){
    #   switch(type,
    #          SL = x,
    #          EL = -1 / x,
    #          ET = exp(x),
    #          CE = 1 / (1 - exp(x)))
    # }
    
    # ginvprime = function(x, type){
    #   switch(type,
    #          SL = rep(1, length(x)),
    #          EL = 1 / x^2,
    #          ET = exp(x),
    #          CE = x * (x - 1))
    # }
    
    gprime1 = function(x, type, del){ # Inverse of gprime(d)
      switch(type,
             SL = rep(1, length(x)),
             EL = x^2,
             ET = x,
             CE = x * (x - 1),
             HD = 2 * x^(1.5),
             Hb = ifelse(abs(x) < del, 1, NA), # ???
             PH = (1 + (x / del)^2)^(1.5))
    }
    
    # Fint = function(lambda, d_S, Z_S, Zbar, type, ..., returnw = F){
    #   Zlambda = Z_S %*% lambda
    #   if(type == "SL"){
    #     return(sum(d_S * (drop(Zlambda))^2 / 2) / N - sum(Zbar %*% lambda))
    #   }else if(type == "EL"){
    #     if(any(Zlambda >= 0)) return(.Machine$double.xmax)
    #     return(sum(d_S * log(-drop(Zlambda))) / N - sum(Zbar %*% lambda))
    #   }else if(type == "ET"){
    #     return(sum(d_S * exp(drop(Zlambda))) / N - sum(Zbar %*% lambda))
    #   }else if(type == "CE"){
    #     if(any(Zlambda >= 0)) return(.Machine$double.xmax)
    #     return(sum(d_S * (drop(Zlambda) - log(1 - exp(drop(Zlambda)) ))) / N - sum(Zbar %*% lambda))
    #   }
    # }
    
    f = function(lambda, d_S, Z_S, Zbar, type, del, ..., returnw = F){
      # w_S = d_S * ginv(drop(Z_S %*% lambda), type = type)
      if(type == "HD" & any(Z_S %*% lambda >= 0)) return(rep(Inf, length(lambda)))
      if(type == "PH" & any(abs(Z_S %*% lambda) >= del)) return(rep(Inf, length(lambda)))
      
      if(type == "SL"){
        w_S = d_S * drop(Z_S %*% lambda)
      }else if(type == "EL"){
        w_S = -d_S / drop(Z_S %*% lambda)
      }else if(type == "ET"){
        w_S = d_S * exp(drop(Z_S %*% lambda))
      }else if(type == "CE"){
        w_S = d_S / (1 - exp(drop(Z_S %*% lambda)))
      }else if(type == "HD"){
        w_S = d_S / drop(Z_S %*% lambda)^2
      }else if(type == "PH"){
        w_S = d_S / sqrt(1 / drop(Z_S %*% lambda)^2 - 1 / del^2)
      }
      if(type != "SL" & any(w_S <= 0)) return(rep(Inf, length(lambda)))
      if(type == "CE" & any(w_S <= 1)) return(rep(Inf, length(lambda)))
      if(returnw == T){
        return(w_S)
      }else{
        return(colSums(Z_S * w_S) / N - Zbar)
      }
    }
    
    h = function(lambda, d_S, Z_S, Z_St, Zbar, type, del){
      # return(Z_St %*% (Z_S * d_S * ginvprime(drop(Z_S %*% lambda), type = type) / N))
      if(type == "SL"){
        return(Z_St %*% (Z_S * d_S) / N)
      }else if(type == "EL"){
        w_S = -d_S / drop(Z_S %*% lambda)
        return(Z_St %*% (Z_S * (w_S^2 / d_S)) / N)
      }else if(type == "ET"){
        w_S = d_S * exp(drop(Z_S %*% lambda))
        return(Z_St %*% (Z_S * w_S) / N)
      }else if(type == "CE"){
        p_Stmp = 1 / (1 - exp(drop(Z_S %*% lambda)))
        return(Z_St %*% (Z_S * d_S * p_Stmp * (p_Stmp - 1)) / N)
      }else if(type == "HD"){
        return(Z_St %*% (Z_S * (-2 * d_S / drop(Z_S %*% lambda)^3)) / N)
      }else if(type == "PH"){
        return(Z_St %*% (Z_S * (d_S * (1 - (drop(Z_S %*% lambda) / del)^2)^(-1.5))) / N)
      }
    }
    
    targetftn = function(W, d_S, Z_S, Z_St, init, Zbar, type, cal, del,..., returnw = F){
      # d_S = rep(1, n) / n
      # Z_S = cbind(1, x_S, log(pi_S)); 
      # Z_St = t(Z_S)
      # init = c(log(n / N), 0, -1)
      # Zbar = c(1, colMeans(x), W)
      if(type == "ET" | type == "SL" | type == "PH"){
        if(W < 0) return(.Machine$double.xmax)
      }else if(type == "EL" | type == "CE" | type == "HD"){
        if(W > 0) return(.Machine$double.xmax)
      }
      alphaHT = Zbar[length(Zbar)]
      Zbar[length(Zbar)] <- W
      nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, Z_S = Z_S, 
                            Z_St = Z_St, Zbar = Zbar, type = type, del = del,
                            method = "Newton", control = list(maxit = 1e5, allowSingular = T))
      # if(nleqslv_res$termcd != 1 & nleqslv_res$termcd != 2){
      if(nleqslv_res$termcd != 1){
        if(max(abs(f(nleqslv_res$x, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del))) > 1e-5)
          return(.Machine$double.xmax)
      }
      w_S = f(nleqslv_res$x, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, returnw = T)
      # if(any(is.infinite(w_S))) return(.Machine$double.xmax)
      
      if(returnw == F){
        if(cal == "GEC1") return(sum(G(w_S, type = type, del = del)) - N * W)
        # else if(cal == "GEC2") return(sum(G(w_S, type = type)) - N * alphaHT * log(abs(W)))
        else if(cal == "GEC2") return(sum(G(w_S, type = type, del = del)) - N * (alphaHT + 1) * log(abs(W + 1)))
        # else if(cal == "GEC2") return(sum(G(w_S, type = type)) - N / g(alphaHT + 1, type = type) * G(abs(W + 1), type = type)) # Not working well
      }else{
        return(w_S)
      }
      
      # if(returnw == F){
      #   if(type == "EL"){
      #     return(-sum(log(w_S)) - N * W)
      #     # return(-sum(log(w_S)) - N^3 / sum(d_S0)^2 * W)
      #     # return(-sum(log(w_S)) - sum(d_S0) * W)
      #     # return(-sum(log(w_S)) - sum(Z_S[,2] * d_S0) / Zbar[2] * W)
      #     # return(-sum(log(w_S)))
      #     # return(-sum(log(w_S)) + n  * log(-W))
      #     # if(W > -1){
      #     #   return(-sum(log(w_S)) - (N - n)  * log(1 + W))
      #     # }else{
      #     #   return(.Machine$double.xmax)
      #     # }
      #   }else if(type == "SL"){
      #     return(sum(w_S^2) / 2 - N * W)
      #     # return(sum(w_S^2) / 2 - sum(d_S0) * W)
      #     # return(sum(w_S^2) / 2 - N^2 / sum(d_S0^2) * W^2 / 2)
      #     # if(W > -1){
      #     #   return(sum(w_S^2) / 2 - N * (1 + sum((d_S0)^2) / N) * log(1 + W))
      #     # }else{
      #     #   return(.Machine$double.xmax)
      #     # }
      #   }else if(type == "ET"){
      #     return(sum((w_S) * (log(w_S) - 1)) - N * W)
      #     # return(sum((w_S) * (log(w_S) - 1)) - N * (1 + sum(d_S0 * log(d_S0)) / N) * log(1 + W))
      #     # if(W > 0){
      #     #   return(sum((w_S) * (log(w_S) - 1)) - N * (1 + sum(d_S0 * log(d_S0)) / N) * log(1 + W))
      #     # }else{
      #     #   return(.Machine$double.xmax)
      #     # }
      #   }
      #   else{
      #     return(-sum(log(w_S)))
      #   }
      # }else{
      #   return(w_S)
      # }
    }
    
    theta_res = NULL
    Var_res = NULL
    
    cal_vec = c("DS", "GEC"); forcntvec = 1:(length(type_vec) * length(cal_vec))
    
    for(forcnt in forcntvec){
      cal = cal_vec[length(cal_vec) - (length(cal_vec) - forcnt) %% length(cal_vec)]
      type = type_vec[(forcnt-1) %/% length(cal_vec) + 1]
      # for(type in type_vec){
      #   for(cal in cal_vec){
      # if(cal == "DS" & type %in% c("CE")) next # CE_Div too unstable...
      # if(type == "SL" & cal != "DS") next # Skip SL
      # if(type == "SL") next # Skip SL
      varnames = paste(type, cal, sep = "_")
      # if(cal == "DS" & type == "SL"){
      #   varnames = "Reg"
      # }
      
      if(cal == "DS" & type == "Hb") next
      
      # if(type == "Hb") Hcnt = Hcnt + 1
      
      if(cal == "DS"){
        d_S = d_S0
      }else{
        d_S = rep(1, n)
      }
      
      if(type == "EL"){
        u_vec = -pi
      }else if(type == "ET"){
        u_vec = -log(pi)
      }else if(type == "SL"){
        u_vec = 1 / pi
      }else if(type == "CE"){
        u_vec = log(1 - pi)
      }else if(type == "HD"){
        u_vec = -sqrt(pi);
      }else if(type == "Hb"){
        u_vec = ifelse(1 / pi > del, del, 1 / pi)
      }else if(type == "PH"){
        u_vec = 1 / pi / sqrt(1 + (1 / pi / del)^2)
      }
      u_vec_S = u_vec[Index_S]; Uhat = mean(u_vec); 
      Z_S = cbind(1, x_S, u_vec_S); Zbar = c(1, colMeans(x), Uhat)
      Z_St = t(Z_S)
      
      init = rep(0, length(Zbar))
      if(cal != "DS"){
        init[length(init)] = 1
        # init = c(0, 0, 1)
      }else{
        if(type == "EL"){
          init[1] = -1
          # init = c(-1, 0)
        }else if(type == "ET"){
          # init = c(0, 0)
        }else if(type == "SL"){
          init[1] = 1
          # init = c(1, 0)
        }else if(type == "CE"){
          init[1] = -1
          # init = c(-1, 0)
        }else if(type == "HD"){
          init[1] = -1
          # init = c(-1, 0)
        }else if(type == "PH"){
          init[1] = 1 / sqrt(1 + 1 / del^2)
          # init = c(1, 0)
        }
        # init = rep(0, length(Zbar))
      }
      # if(type %in% c("Hb")){
      #   w = CVXR::Variable(length(d_S0))
      #   constraints <- list(Z_St %*% w / N == Zbar)
      #   
      #   # if(cal == "DS"){
      #   #   Phi_R <- CVXR::Minimize(sum(d_S * huber(w / d_S - 1, M = del)))
      #   # }else if(cal == "GEC"){
      #   #   Phi_R <- CVXR::Minimize(sum(huber(w, M = del)))
      #   # }
      #   
      #   if(type == "Hb"){
      #     Phi_R <- CVXR::Minimize(sum(huber(w, M = del)))          
      #   # }else if(type == "PH"){
      #   #   Phi_R <- CVXR::Minimize(sum(cvxr_norm(hstack(rep(1, length(d_S0)), w / del), p = 2, axis = 1)))          
      #   }
      # 
      #   p <- CVXR::Problem(Phi_R, constraints)
      #   res <- CVXR::solve(p, solver = "ECOS_BB")
      #   
      #   w_S = drop(res$getValue(w))
      #   # drop(Z_St %*% res$getValue(w) / N - Zbar)
      #   # plot(res$getValue(w), d_S0)
      #   # summary(res$getValue(w))
      # }else 
      if(cal != "DS" & cal != "GEC"){
        # if(type == "SL") stop(Zbar[length(Zbar)])
        nlmres= nlm(targetftn, Zbar[length(Zbar)], d_S = d_S, Z_S = Z_S, Z_St = Z_St, 
                    init = init, Zbar = Zbar, type = type, cal = cal, del = del)
        # if(nlmres$code != 1 & nlmres$code != 2 & nlmres$code != 3) print(nlmres$estimate)
        if(nlmres$code != 1 & nlmres$code != 2 & nlmres$code != 3) stop(nlmres$code)
        W = nlmres$estimate
        if(nlmres$minimum >= .Machine$double.xmax){
          # stop(targetftn(Zbar[length(Zbar)], d_S = d_S, Z_S = Z_S, Z_St = Z_St, init = init, Zbar = Zbar, type = type))
          w_S = NA
        }else{
          w_S = targetftn(W, d_S, Z_S, Z_St, init, Zbar, type, returnw = T, cal = cal, del = del)
        }
        # if(type == "ET") stop(W)
        # if(any(is.infinite(sum(y_S * w_S)))) stop(summary(w_S))
        # if(any(is.infinite(sum(y_S * w_S)))) stop(nlmres)
      }else if(cal != "DS" | !(type %in% c("CE"))){
        # Fint(init, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type, control = list(maxit = 10000))
        # optim(init, Fint, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type, control = list(maxit = 10000), method = "L-BFGS-B")
        # f(init, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type)
        # if(any(!is.finite(Zbar))) stop(paste(cal, type))
        
        # library(CVXR)
        # w = Variable(length(d_S0))
        # constraints <- list(Z_St %*% w / N == Zbar, w >= 1)
        # # Phi_R <- Maximize(sum(entr(w) - entr(w - 1)))
        # # Phi_R <- Minimize(sum(kl_div(w-1, w) - log(w)))
        # Phi_R <- Minimize(sum(d_S0 * (kl_div(w / d_S0-1, w / d_S0) - log(w / d_S0))))
        # p <- Problem(Phi_R, constraints)          
        # res <- CVXR::solve(p, solver = "ECOS_BB")
        
        nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, Z_S = Z_S, 
                              Z_St = Z_St, Zbar = Zbar, type = type, del = del,
                              method = "Newton", control = list(maxit = 1e5, allowSingular = T))
        # print(type); print(cal); print(nleqslv_res)
        # if(cal == "GEC" & type == "ET") print(nleqslv_res$x)
        if(nleqslv_res$termcd != 1){
          if(max(abs(f(nleqslv_res$x, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del))) > 1e-5)
            w_S = NA
          # stop(c(paste(type, cal, sep = "_"), nleqslv_res))
        }else{
          # if(type == "CE" & cal == "DS") stop(c(paste(type, cal, sep = "_"), nleqslv_res))
          w_S = f(nleqslv_res$x, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, returnw = T)             
        }
      }
      
      if(cal == "DS" & (type %in% c("CE")) ){
        w_S = drop(d_S0 + (Z_S * d_S0) %*% solve(Z_St %*% (Z_S * d_S0), Zbar * N - t(Z_S) %*% d_S0))
      }
      
      theta_res = c(theta_res, setNames(sum(y_S * w_S) / N, varnames))
      
      pimat_S = diag(1 - pi_S) * (w_S / N)^2 # 1 / pi_i * (1 - 1 / pi_i)
      # pimat_S = diag(1 - pi_S) * (w_S / sum(d_S0))^2 # 1 / pi_i * (1 - 1 / pi_i)
      # if(type != "EL" & type != "CE") pimat_S = diag(w_S - 1) * w_S / N^2 # 1 / pi_i * (1 - 1 / pi_i)
      
      if(cal == "DS"){
        gammahat = solve(Z_St %*% (Z_S * d_S0), Z_St %*% (y_S * d_S0))
      }else{
        alphaHT = sum(u_vec_S / pi_S) / N
        
        hatSigmazz = Z_St %*% (Z_S * gprime1(d_S0, type, del = del)) / N
        xcolums = 1:(ncol(Z_S) - 1); gcolums = ncol(Z_S)
        hatSigmaxx = hatSigmazz[xcolums, xcolums]
        hatSigmagx = hatSigmazz[gcolums, xcolums, drop = F]
        hatSigmagg = hatSigmazz[gcolums, gcolums, drop = F]
        hatSigmagg_x = drop(hatSigmagg - hatSigmagx %*% solve(hatSigmaxx, t(hatSigmagx)))
        
        gammahat = solve(hatSigmazz * N, Z_St %*% (y_S * gprime1(d_S0, type, del = del)))
      }
      
      # if(cal == "DS" & type == "CE"){
      #   theta_res[length(theta_res)] = drop(t(Zbar) %*% gammahat) + sum((y_S - drop(Z_S %*% gammahat)) * d_S0) / N
      # }
      
      if(cal == "DS"){
        Varhat = crossprod(y_S - drop(Z_S %*% gammahat), pimat_S %*% (y_S - drop(Z_S %*% gammahat))) 
        
      }else if(cal == "GEC"){
        Varhat = crossprod(y_S - drop(Z_S %*% gammahat), pimat_S %*% (y_S - drop(Z_S %*% gammahat))) 
      }else if(cal == "GEC1"){
        Z_S2 = Z_S; Z_S2[,ncol(Z_S2)] <- c(hatSigmagx %*% solve(hatSigmaxx, Z_St[xcolums,]))
        Varhat = crossprod(y_S - drop(Z_S2 %*% gammahat), pimat_S %*% (y_S - drop(Z_S2 %*% gammahat)))
      }else if(cal == "GEC2"){
        # Z_S2 = Z_S; Z_S2[,3] <- (1 / hatSigmagg_x / (1 / alphaHT + 1 / hatSigmagg_x)) * c(hatSigmagx %*% solve(hatSigmaxx, Z_St[xcolums,]))
        Z_S2 = Z_S; Z_S2[,ncol(Z_S2)] <- (1 / hatSigmagg_x / (1 / (alphaHT + 1) + 1 / hatSigmagg_x)) * c(hatSigmagx %*% solve(hatSigmaxx, Z_St[xcolums,]))
        Varhat = crossprod(y_S - drop(Z_S2 %*% gammahat), pimat_S %*% (y_S - drop(Z_S2 %*% gammahat)))
      }
      Var_res = c(Var_res, setNames(Varhat, varnames))
      # }
    }

    cal_vec = c("DS", "GEC1", "GEC2"); forcntvec = 1:(length(type_vec) * length(cal_vec))
    
    for(forcnt in forcntvec){
      cal = cal_vec[length(cal_vec) - (length(cal_vec) - forcnt) %% length(cal_vec)]
      type = type_vec[(forcnt-1) %/% length(cal_vec) + 1]
      # type = type_vec[length(type_vec) - (length(type_vec) - forcnt) %% length(type_vec)]
      # cal = cal_vec[(forcnt-1) %/% length(type_vec) + 1]
      # for(type in type_vec){
      #   for(cal in cal_vec){
      # if(cal == "DS" & type %in% c("ET", "CE", "HD")) next # CE_DS too unstable...
      if(type == "SL" & cal != "DS") next # Skip SL
      varnames = paste(type, cal, sep = "_")
      # if(cal != "DS"){
      #   varnames = paste(type, cal, sep = "_")
      # }else if(type == "EL"){
      #   varnames = "PEL"
      # }else if(type == "SL"){
      #   varnames = "Reg"
      # }
      
      if(cal == "DS" & type == "Hb") next
      
      if(cal == "DS"){
        d_S = d_S0; Z_S = cbind(1, x_S); Z_St = t(Z_S); Zbar = c(1, colMeans(x))
      }else{
        d_S = rep(1, n)
        if(type == "EL"){
          u_vec = -pi
        }else if(type == "ET"){
          u_vec = -log(pi)
        }else if(type == "SL"){
          u_vec = 1 / pi
        }else if(type == "CE"){
          u_vec = log(1 - pi)
        }else if(type == "HD"){
          u_vec = -sqrt(pi);
        }else if(type == "Hb"){
          u_vec = ifelse(1 / pi > del, del, 1 / pi)
        }else if(type == "PH"){
          u_vec = 1 / pi / sqrt(1 + (1 / pi / del)^2)
        }
        u_vec_S = u_vec[Index_S]; Uhat = mean(u_vec); 
        if(cal == "Knl"){
          # bw1 <- np::npregbw(reformulate(colnames(data_S[,-1,drop = F]), response = "d_S0"),
          #                    regtype = "lc", data = data_S[,-1,drop = F])
          
          bw1 <- np::npregbw(reformulate(colnames(data_S[,-1,drop = F]), response = "d_S0 * u_vec_S"),
                             regtype = "lc", data = data_S[,-1,drop = F])
          
          # kre1 <- np::npreg(bws = bw1)
          # uhat1 = predict(kre1, newdata = data); # mean(1 / uhat1)
          uvec = npksum(txdat=drop(x_S), tydat= d_S0 * u_vec_S , exdat = drop(x), bws=bw1$bw)$ksum/
            npksum(txdat=drop(x_S), tydat= d_S0, exdat = drop(x), bws=bw1$bw)$ksum
          # Uhat = mean(uvec[!is.nan(uvec)])
          # uvec[Index_S] <- u_vec_S
          # Uhat = mean(uvec)
          Uhat = mean(uvec) + sum((u_vec_S  - uvec[Index_S]) / pi_S) / sum(1 / pi_S)
          # # u_vec_S = uvec[Index_S]
          
          # Uhat = sum(u_vec_S * d_S0) / N# To be removed
        }else if(cal == "Knl2"){
          bw1 <- np::npregbw(reformulate(colnames(data_S[,-1,drop = F]), response = "1 / pi_S * (1 / pi_S - 1)"),
                             regtype = "lc", data = data_S[,-1,drop = F])
          # plot(bw1)
          # kre1 <- np::npreg(bws = bw1)
          # uhat1 = predict(kre1, newdata = data); # mean(1 / uhat1)
          uvec = npksum(txdat=drop(x_S), tydat= d_S0 * (d_S0 - 1) * u_vec_S , exdat = drop(x), bws=bw1$bw)$ksum/
            npksum(txdat=drop(x_S), tydat= d_S0 * (d_S0 - 1), exdat = drop(x), bws=bw1$bw)$ksum
          
          # plot(uvec, u_vec)
          # plot(x_S[,2], uvec[Index_S])
          # plot(x_S[,2], u_vec[Index_S])
          
          # Uhat = mean(uvec[!is.nan(uvec)]) + sum((u_vec_S  - uvec[Index_S]) / pi_S) / N
          # uvec[Index_S] <- u_vec_S
          # Uhat = mean(uvec) + sum((u_vec_S  - uvec[Index_S]) / pi_S) / N
          Uhat = mean(uvec) + sum((u_vec_S  - uvec[Index_S]) / pi_S) / sum(1 / pi_S)
          # u_vec_S = uvec[Index_S]
          
          # uvec0 = npksum(txdat=drop(x), tydat= (1 / pi - 1) * u_vec , exdat = drop(x), bws=bw0$bw)$ksum/
          #   npksum(txdat=drop(x), tydat= (1 / pi - 1), exdat = drop(x), bws=bw0$bw)$ksum
        }
        Z_S = cbind(1, x_S, u_vec_S); Zbar = c(1, colMeans(x), Uhat)
        Z_St = t(Z_S)
      }
      init = rep(0, length(Zbar))
      if(cal != "DS"){
        init[length(init)] = 1
        # init = c(0, 0, 1)
      }else{
        if(type == "EL"){
          init[1] = -1
          # init = c(-1, 0)
        }else if(type == "ET"){
          # init = c(0, 0)
        }else if(type == "SL"){
          init[1] = 1
          # init = c(1, 0)
        }else if(type == "CE"){
          init[1] = -1
          # init = c(-1, 0)
        }else if(type == "HD"){
          init[1] = -1
          # init = c(-1, 0)
        }else if(type == "PH"){
          init[1] = 1 / sqrt(1 + 1 / del^2)
          # init = c(1, 0)
        }
        # init = rep(0, length(Zbar))
      }
      #   if(type == "Hb"){
      #     Phi_R <- CVXR::Minimize(sum(huber(w, M = del)))          
      #   # }else if(type == "PH"){
      #   #   Phi_R <- CVXR::Minimize(sum(cvxr_norm(hstack(rep(1, length(d_S0)), w / del), p = 2, axis = 1)))          
      #   }
      # 
      #   p <- CVXR::Problem(Phi_R, constraints)
      #   res <- CVXR::solve(p, solver = "ECOS_BB")
      #   
      #   w_S = drop(res$getValue(w))
      #   # drop(Z_St %*% res$getValue(w) / N - Zbar)
      #   # plot(res$getValue(w), d_S0)
      #   # summary(res$getValue(w))
      # }else 
      if(cal != "DS" & cal != "GEC0" & cal != "Knl" & cal != "Knl2"){
        # if(type == "SL") stop(Zbar[length(Zbar)])
        nlmres= nlm(targetftn, Zbar[length(Zbar)], d_S = d_S, Z_S = Z_S, Z_St = Z_St, 
                    init = init, Zbar = Zbar, type = type, cal = cal, del = del)
        # if(nlmres$code != 1 & nlmres$code != 2 & nlmres$code != 3) print(nlmres$estimate)
        if(nlmres$code != 1 & nlmres$code != 2 & nlmres$code != 3) stop(nlmres$code)
        W = nlmres$estimate
        if(nlmres$minimum >= .Machine$double.xmax){
          # stop(targetftn(Zbar[length(Zbar)], d_S = d_S, Z_S = Z_S, Z_St = Z_St, init = init, Zbar = Zbar, type = type))
          w_S = NA
        }else{
          w_S = targetftn(W, d_S, Z_S, Z_St, init, Zbar, type, returnw = T, cal = cal, del = del)
        }
        # if(type == "ET") stop(W)
        # if(any(is.infinite(sum(y_S * w_S)))) stop(summary(w_S))
        # if(any(is.infinite(sum(y_S * w_S)))) stop(nlmres)
      }else if(cal != "DS" | !(type %in% c("CE"))){
        # Fint(init, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type, control = list(maxit = 10000))
        # optim(init, Fint, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type, control = list(maxit = 10000), method = "L-BFGS-B")
        # f(init, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type)
        # if(any(!is.finite(Zbar))) stop(paste(cal, type))
        
        # library(CVXR)
        # w = Variable(length(d_S0))
        # constraints <- list(Z_St %*% w / N == Zbar, w >= 1)
        # # Phi_R <- Maximize(sum(entr(w) - entr(w - 1)))
        # # Phi_R <- Minimize(sum(kl_div(w-1, w) - log(w)))
        # Phi_R <- Minimize(sum(d_S0 * (kl_div(w / d_S0-1, w / d_S0) - log(w / d_S0))))
        # p <- Problem(Phi_R, constraints)          
        # res <- CVXR::solve(p, solver = "ECOS_BB")
        
        nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, Z_S = Z_S, 
                              Z_St = Z_St, Zbar = Zbar, type = type, del = del,
                              method = "Newton", control = list(maxit = 1e5, allowSingular = T))
        # print(type); print(cal); print(nleqslv_res)
        # if(cal == "GEC0" & type == "ET") print(nleqslv_res$x)
        if(nleqslv_res$termcd != 1){
          if(max(abs(f(nleqslv_res$x, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del))) > 1e-5)
            w_S = NA
          # stop(c(paste(type, cal, sep = "_"), nleqslv_res))
        }else{
          # if(type == "CE" & cal == "DS") stop(c(paste(type, cal, sep = "_"), nleqslv_res))
          w_S = f(nleqslv_res$x, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, returnw = T)               
        }
      }
      
      if(cal == "DS" & (type %in% c("CE")) ){
        w_S = drop(d_S0 + (Z_S * d_S0) %*% solve(Z_St %*% (Z_S * d_S0), Zbar * N - t(Z_S) %*% d_S0))
      }
      
      theta_res = c(theta_res, setNames(sum(y_S * w_S) / N, varnames))
      
      pimat_S = diag(1 - pi_S) * (w_S / N)^2 # 1 / pi_i * (1 - 1 / pi_i)
      # pimat_S = diag(1 - pi_S) * (w_S / sum(d_S0))^2 # 1 / pi_i * (1 - 1 / pi_i)
      # if(type != "EL" & type != "CE") pimat_S = diag(w_S - 1) * w_S / N^2 # 1 / pi_i * (1 - 1 / pi_i)
      
      if(cal %in% c("DS")){
        gammahat = solve(Z_St %*% (Z_S * d_S0), Z_St %*% (y_S * d_S0))
      }else{
        if(cal != "DS" & cal != "GEC0" & cal != "Knl" & cal != "Knl2"){
          alphaHT = W
        }else{
          alphaHT = sum(u_vec_S / pi_S) / N
        }
        
        hatSigmazz = Z_St %*% (Z_S * gprime1(d_S0, type, del = del)) / N
        # Z = cbind(1, x, u_vec)   # True m, To be changed
        # hatSigmazz = t(Z) %*% (Z * gprime1(1 / pi, type, del = del) * pi) / N # To be changed
        
        xcolums = 1:(ncol(Z_S) - 1); gcolums = ncol(Z_S)
        hatSigmaxx = hatSigmazz[xcolums, xcolums]
        hatSigmagx = hatSigmazz[gcolums, xcolums, drop = F]
        hatSigmagg = hatSigmazz[gcolums, gcolums, drop = F]
        hatSigmagg_x = drop(hatSigmagg - hatSigmagx %*% solve(hatSigmaxx, t(hatSigmagx)))
        
        gammahat = solve(hatSigmazz * N, Z_St %*% (y_S * gprime1(d_S0, type, del = del)))
        # gammahat = solve(hatSigmazz * N, t(Z) %*% (y * gprime1(1 / pi, type, del = del) * pi)) # To be changed
      }
      
      # if(cal == "DS" & type == "CE"){
      #   theta_res[length(theta_res)] = drop(t(Zbar) %*% gammahat) + sum((y_S - drop(Z_S %*% gammahat)) * d_S0) / N
      # }
      
      # if(cal == "GEC2") stop((alphaHT + 1) / (hatSigmagg_x + alphaHT + 1))
      
      if(cal %in% c("DS", "GEC0")){
        Varhat = crossprod(y_S - drop(Z_S %*% gammahat), pimat_S %*% (y_S - drop(Z_S %*% gammahat))) 
      }else if(cal == "GEC1"){
        Z_S2 = Z_S; Z_S2[,ncol(Z_S2)] <- c(hatSigmagx %*% solve(hatSigmaxx, Z_St[xcolums,]))
        Varhat = crossprod(y_S - drop(Z_S2 %*% gammahat), pimat_S %*% (y_S - drop(Z_S2 %*% gammahat)))
      }else if(cal == "GEC2"){
        # Z_S2 = Z_S; Z_S2[,ncol(Z_S2)] <- (1 / hatSigmagg_x / (1 / (alphaHT + 1) + 1 / hatSigmagg_x)) * c(hatSigmagx %*% solve(hatSigmaxx, Z_St[xcolums,]))
        Z_S2 = Z_S; Z_S2[,ncol(Z_S2)] <- (alphaHT + 1) / (hatSigmagg_x + alphaHT + 1) * c(hatSigmagx %*% solve(hatSigmaxx, Z_St[xcolums,]))
        Varhat = crossprod(y_S - drop(Z_S2 %*% gammahat), pimat_S %*% (y_S - drop(Z_S2 %*% gammahat)))
      }else if(cal == "Knl" | cal == "Knl2"){
        # Varhat = NA
        Z_S2 = Z_S; Z_S2[,ncol(Z_S2)] <- uvec[Index_S]
        # Z_S2 = Z_S; Z_S2[,ncol(Z_S2)] <- uvec0[Index_S] # Not valid
        Varhat = crossprod(y_S - drop(Z_S2 %*% gammahat), pimat_S %*% (y_S - drop(Z_S2 %*% gammahat)))
        
        # Varhat = crossprod(y_S - drop(Z_S %*% gammahat), pimat_S %*% (y_S - drop(Z_S %*% gammahat))) +
        #   crossprod(u_vec_S - uvec[Index_S], pimat_S %*% (u_vec_S - uvec[Index_S])) * gammahat[length(gammahat)]^2
        
        # Variance estimation using eta ####
        # Z = cbind(1, x, uvec)  # estimated m
        # # Z = cbind(1, x, u_vec)   # True m
        # eta = drop(Z %*% gammahat)
        # eta[Index_S] = eta[Index_S] + (y_S - drop(Z[Index_S,] %*% gammahat)) * d_S0
        # # eta[Index_S] = eta[Index_S] + (y_S - drop(Z[Index_S,] %*% gammahat)) * w_S
        # Varhat = var(eta) / N
        
        # crossprod(y_S - drop(cbind(cbind(1, x)[Index_S,], u_vec_S) %*% gammahat), pimat_S %*% (y_S - drop(cbind(cbind(1, x)[Index_S,], u_vec_S) %*% gammahat)))
        
        # eta2 = uvec
        # eta2[Index_S] = eta2[Index_S] + (u_vec_S - uvec[Index_S]) / d_S0
        # Varhat = crossprod(y_S - drop(Z_S %*% gammahat), pimat_S %*% (y_S - drop(Z_S %*% gammahat))) +
        #   var(eta2) * gammahat[length(gammahat)]^2 / N
        
      }
      
      # if(cal == "GEC1" & type == "PH"){ # To be corrected
      #   Z_S2 = Z_S; Z_S2[,ncol(Z_S2)] <- (alphaHT + 1) / (hatSigmagg_x + alphaHT + 1) * c(hatSigmagx %*% solve(hatSigmaxx, Z_St[xcolums,]))
      #   theta_res[length(theta_res)] = drop(t(Zbar) %*% gammahat) + sum((y_S - drop(Z_S %*% gammahat)) * d_S0) / N
      # }
      
      # if(cal == "Knl"){ # To be removed
      #   # Varhat = NA
      #   Z_S2 = Z_S; Z_S2[,ncol(Z_S2)] <- 0
      #   Varhat = crossprod(y_S - drop(Z_S2 %*% gammahat), pimat_S %*% (y_S - drop(Z_S2 %*% gammahat))) 
      # }
      Var_res = c(Var_res, setNames(Varhat, varnames))
      # }
      
      
      
      # if(cal == "Knl2"){ # To be removed
      #   Zbar[length(Zbar)] = mean(uvec)
      #   theta_res[length(theta_res)] = drop(t(Zbar) %*% gammahat) + sum((y_S - drop(Z_S2 %*% gammahat)) * d_S0) / N
      # }
      
    }
    
    theta_vec = c("HT" = theta_HT, "Hajek" = theta_Hajek, theta_res) - theta
    var_vec = c("HT" = Var_HT, "Hajek" = Var_Hajek, Var_res)
    CR_vec = ifelse(abs(theta_vec) < qnorm(0.975) * sqrt(var_vec), 1, 0)
    
    list(theta = theta_vec, 
         Var = var_vec,
         CR = CR_vec)
  }
final_res1 = lapply(final_res, function(x) x[[1]])

stopCluster(cl)
timenow2 = Sys.time()
print(timenow2 - timenow1)
# if(sum(!sapply(final_res1, function(x) is.numeric(unlist(x)))) != 0) stop(paste(final_res1))
paste("# of failure:", sum(!sapply(final_res1, function(x) is.numeric(unlist(x))))); final_res0 = final_res1
print(final_res0[!sapply(final_res0, function(x) is.numeric(unlist(x)))])
final_res1 = final_res1[sapply(final_res1, function(x) is.numeric(unlist(x)))]
res1 = do.call("rbind", final_res1)

print(res1[apply(res1, 1, function(x) any(is.na(x) | is.infinite(x))),])
colnames(res1)

Bias = colMeans(res1, na.rm = T)
SE = apply(res1, 2, function(x) sqrt(var(x, na.rm = T) * (sum(!is.na(x))-1)/sum(!is.na(x)) ))
RMSE = apply(res1, 2, function(x) sqrt(mean(x^2, na.rm = T)))

round(cbind(Bias, SE, RMSE), 4)

final_res2 = lapply(final_res[sapply(final_res0, function(x) is.numeric(unlist(x)))], function(x) x[[2]])
final_res3 = lapply(final_res[sapply(final_res0, function(x) is.numeric(unlist(x)))], function(x) x[[3]])

res2 = do.call("rbind", final_res2)
print(res2[apply(res2, 1, function(x) any(is.na(x) | is.infinite(x))),])
# colMeans(sqrt(res2), na.rm = T)
RB = (colMeans(res2, na.rm = T) - SE^2) / SE^2

res3 = do.call("rbind", final_res3)
print(res3[apply(res3, 1, function(x) any(is.na(x) | is.infinite(x))),])
CR = colMeans(sqrt(res3), na.rm = T)

round(cbind(Bias, SE, RMSE, RB, CR), 4)
# round(cbind(Bias, SE, RMSE, RB, CR), 6)

# xtable::xtable(cbind(Bias, SE, RMSE, RB, CR), digits = c(1,4,4,4,2,3))
xtable::xtable(cbind(cbind(Bias, SE, RMSE) * 1e2, RB, CR), digits = c(1,3,3,3,2,3))

cbind(`SB(%)` = Bias[-1] / RMSE[-1] * 100, `R-RMSE` = RMSE[-1] / RMSE[2] * 100, `CR(%)` = CR[-1] * 100)
xtable::xtable(cbind(`SB(%)` = Bias[-1] / RMSE[-1] * 100, `R-RMSE` = RMSE[-1] / RMSE[2] * 100, `CR(%)` = CR[-1] * 100),
               digits = c(1,0,0,1))

# png(filename=paste(timenow0, ".png", sep = ""), width = 640, height = 360)
boxplot(res1[,-c(1,2)] + theta, ylab = "Estimated population mean", xaxt = "n")
# abline(h = theta, col = "red", lty = 3)
# abline(v = c(0.5, 2.5, 6.5, 10.5, 14.5), col = "blue")
# 
# axis(side = 3, at = c(1, 2, 4.5, 8.5, 12.5), labels = c("PEL", "Reg", cal_vec[-1]))
# label_tmp = c(rep(" ", 2), rep(type_vec[-1], length(cal_vec[-1])))
# axis(side = 1, at = 1:length(label_tmp), labels = label_tmp) 
# dev.off()
