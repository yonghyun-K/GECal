# Model-assisted calibration using 
# Generalized entropy calibration in survey sampling

# Simulation setup

if(!interactive()){
  args <- as.numeric(commandArgs(trailingOnly = TRUE))
}else{
  args <- c(50)
}

timenow1 = Sys.time()
timenow0 = gsub(' ', '_', gsub('[-:]', '', timenow1))
timenow = paste(timenow0, ".txt", sep = "")

# install.packages( "MatrixModels", type="win.binary" )  

library(nleqslv)
library(CVXR)
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
library(caret)
library(tidyverse)
library(kableExtra)

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

# modelcases = expand.grid(c(T,F),c(T,F),c(T,F))
# modelcases = expand.grid(c(T,F),c(T,F))
# modelcases = expand.grid(c(T,F))
# rnames <- apply(apply(modelcases, 2, ifelse, "C", "M"), 1, paste, collapse = "")

modelcases = expand.grid(c(T))
rnames <- c("C")

# # Kang & Schafer, 2007
# n = 1000
# pcol = 4; # pcol = 4 is True
# x = matrix(rnorm(n * pcol, 0, 1), nc= pcol)
# 
# pi = 1 / (1 + exp(-(-x[,1] + 0.5 * x[,2] -
#                       0.25 * x[,3] - 0.1 * x[,4]))) # True
# # pi = runif(n, 0, 1) # To be removed
# # pi = rep(0.3, n) # To be removed
# 
# 
# # vx = (x[,1]^2 + exp(x[,3]) / 1.65) * 80
# vx = exp(x[,1] / 2 + x[,3]) * 70
# # vx = rep(1, n) # To be removed
# 
# e = rnorm(n, 0, sqrt(vx))
# 
# y = 210 + 27.4 * x[,1] + 13.7 * x[,2] +
#   13.7 * x[,3]+ 13.7 * x[,4] + e
# 
# theta = sum(y)
# 
# # pi = 1 / (1 + exp(-(-x[,1] + 0.5 * x[,2] -
# #                       0.25 * x[,3] - 0.1 * x[,4] - e / 40)))

# Additional simulation study
n = 5000
pcol = 3; # pcol = 4 is True
x = matrix(rnorm(n * pcol, 2, 1), nc= pcol)

# vx = (x[,1]^2 + exp(x[,3]) / 1.65) * 80
vx = exp(-x[,1] + x[,3])
# vx = rep(1, n) # To be removed

e = rnorm(n, 0, sqrt(vx))

y = 2 + 1 * x[,1] + 2 * x[,2] + 3 * x[,3] + e; print("OR Model 1") # Model 1
# y = 2 + x[,1]^2 + + exp(x[,2]) + e; print("OR Model 2") # Model 2

# sd(e); sd(y - e)

pi = 1 / (1 + exp(-(0.1 * x[,1] - 0.1 * x[,2] - 0.2 * x[,3]))); print("RM Model 1") # Model 1
# pi = 1 / (1 + exp(-(-0.1 * x[,1] + 0.1 * e))); print("RM Model 2") # Model 2
# pi = runif(n, 0, 1) # To be removed
# pi = rep(0.3, n) # To be removed

theta = sum(y)
mean(pi)


final_res <- foreach(
  simnum = 1:SIMNUM, 
  .packages = c("nleqslv", "CVXR", "caret"),
  .errorhandling="pass") %dopar% {
    theta_mat = NULL
    var_mat = NULL; var_mat2 = NULL
    CR_mat = NULL; CR_mat2 = NULL
    alpha_vec = NULL
    for(cnt in 1:nrow(modelcases)){
      # RM = modelcases[cnt, 1]
      OM = modelcases[cnt, 1]
      # VM = modelcases[cnt, 3]
      
      # Simulation setup of
      
      
      theta_res = NULL
      var_res = NULL; var_res2 = NULL
      CR_res = NULL
      
      delta = rbinom(n, 1, pi)
      Index_S = (delta == 1)
      
      y_S = y[Index_S]
      x_S = x[Index_S,]
      
      del = NA
      # type = "EL"
      
      data = data.frame(y, delta, vx = vx, x = x, pi = pi)
      data_S = data[Index_S,]
      
      theta_Hajek = sum(y_S / pi[Index_S]) / sum(1 / pi[Index_S]) * n
      
      theta_res = c(theta_res, Hajek = theta_Hajek)  # Hajek
      
      var1 = sum((y_S - theta_Hajek / n)^2 / pi[Index_S] * (1 / pi[Index_S] - 1))
      var_res = c(var_res, Hajek = var1)
      var_res2 = c(var_res2, Hajek = var1)
      
      # Omodel = lm(reformulate(paste0("x.", 1:pcol), response = "y"), data = data_S)
      # yhat = predict.lm(Omodel, data, type = "response")
      # theta_res = c(theta_res, GREG0 = sum(yhat) + sum((y_S - yhat[Index_S]) / pi[Index_S])) # GREG0
      # 
      # Omodel = lm(reformulate(paste0("x.", 1:pcol), response = "y"), weights = 1 / pi, data = data_S)
      # yhat = predict.lm(Omodel, data, type = "response")
      # theta_res = c(theta_res, GREG1 = sum(yhat) + sum((y_S - yhat[Index_S]) / pi[Index_S])) # GREG1
      # 
      # Omodel = lm(reformulate(paste0("x.", 1:pcol), response = "y"), weights = 1 / vx, data = data_S)
      # yhat = predict.lm(Omodel, data, type = "response")
      # theta_res = c(theta_res, GREG2 = sum(yhat) + sum((y_S - yhat[Index_S]) / pi[Index_S])) # GREG2
      
      
      # solve(t(cbind(1, x_S)) %*% (cbind(1, x_S) / pi[Index_S]), t(cbind(1, x_S)) %*% (y_S / pi[Index_S]))
      # c(cbind(1, x) %*% solve(t(cbind(1, x_S)) %*% (cbind(1, x_S) / vx[Index_S]), t(cbind(1, x_S)) %*% (y_S / vx[Index_S]))) - yhat
      
      # var_res = c(var_res, AIPW = var(eta) / n)
      
      # Z_S = cbind(1, x_S); Zbar = c(1, colMeans(x));  Z_St = t(Z_S); d_S = rep(1, sum(delta)); 
      
      prime1 = function(x, type, del){ # Inverse of gprime(d)
        switch(type,
               SL = rep(1, length(x)),
               EL = x^2,
               ET = x,
               CE = x * (x - 1),
               HD = 2 * x^(1.5),
               Hb = ifelse(abs(x) < del, 1, NA), # ???
               PH = (1 + (x / del)^2)^(1.5),
               GE = x^(1 - del))
      }
      
      G = function(x, type, del){
        switch(type,
               SL = x^2/2,
               EL = -log(x),
               ET = x * (log(x) - 1),
               CE = (x-1) * log(x-1) - x * log(x),
               HD = -2 * sqrt(x),
               PH = del^2 * sqrt(1 + (x / del)^2))
      }
      
      f = function(lambda, d_S, v_S, Z_S, Zbar, type, del, n, ..., returnw = F){
        # print(lambda)
        # w_S = d_S * ginv(drop(Z_S %*% lambda), type = type)
        if(type == "HD" & any(Z_S %*% lambda / v_S >= 0)) return(rep(Inf, length(lambda)))
        if(!returnw & type == "GE" & any(drop(Z_S %*% lambda / v_S) * del <= 0)) return(rep(Inf, length(lambda)))
        if(type == "PH" & any(abs(Z_S %*% lambda / v_S) >= del)) return(rep(Inf, length(lambda)))
        
        if(type == "SL"){
          w_S = d_S * drop(Z_S %*% lambda / v_S)
        }else if(type == "EL"){
          w_S = -d_S / drop(Z_S %*% lambda / v_S)
        }else if(type == "ET"){
          w_S = d_S * exp(drop(Z_S %*% lambda / v_S))
        }else if(type == "CE"){
          w_S = d_S / (1 - exp(drop(Z_S %*% lambda / v_S)))
        }else if(type == "HD"){
          w_S = d_S / drop(Z_S %*% lambda / v_S)^2
        }else if(type == "PH"){
          w_S = d_S / sqrt(1 / drop(Z_S %*% lambda / v_S)^2 - 1 / del^2)
        }else if(type == "GE"){
          w_S = d_S * (drop(Z_S %*% lambda / v_S) * del)^(1 / del)
        }
        # print(w_S)
        if(!returnw & type != "SL" & any(w_S <= 0)) return(rep(Inf, length(lambda)))
        if(type == "CE" & any(w_S <= 1)) return(rep(Inf, length(lambda)))
        if(returnw == T){
          return(w_S)
        }else{
          return(colSums(Z_S * w_S) / n - Zbar)
        }
      }
      
      h = function(lambda, d_S, v_S, Z_S, Z_St, Zbar, type, del, n){
        # return(Z_St %*% (Z_S * d_S * ginvprime(drop(Z_S %*% lambda), type = type) / n))
        if(type == "SL"){
          return(Z_St %*% (Z_S / v_S * d_S) / n)
        }else if(type == "EL"){
          w_S = -d_S / drop(Z_S %*% lambda / v_S)
          return(Z_St %*% (Z_S / v_S * (w_S^2 / d_S)) / n)
        }else if(type == "ET"){
          w_S = d_S * exp(drop(Z_S %*% lambda / v_S))
          return(Z_St %*% (Z_S / v_S * w_S) / n)
        }else if(type == "CE"){
          p_Stmp = 1 / (1 - exp(drop(Z_S %*% lambda / v_S)))
          return(Z_St %*% (Z_S / v_S * d_S * p_Stmp * (p_Stmp - 1)) / n)
        }else if(type == "HD"){
          return(Z_St %*% (Z_S / v_S * (-2 * d_S / drop(Z_S %*% lambda / v_S)^3)) / n)
        }else if(type == "PH"){
          return(Z_St %*% (Z_S / v_S * (d_S * (1 - (drop(Z_S %*% lambda / v_S) / del)^2)^(-1.5))) / n)
        }else if(type == "GE"){
          return(Z_St %*% (Z_S / v_S * (d_S * drop(Z_S %*% lambda  / v_S * del)^(1 / del - 1)) / n))
        }
      }
      
      targetftn0 = function(W, d_S, Z_S, Z_St, init, Zbar, type, del,..., returnw = F){
        if(type == "ET" | type == "SL" | type == "PH"){
          if(W < 0) return(.Machine$double.xmax)
        }else if(type == "EL" | type == "CE" | type == "HD"){
          if(W > 0) return(.Machine$double.xmax)
        }
        alphaHT = Zbar[length(Zbar)]
        Zbar[length(Zbar)] <- W
        nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, Z_S = Z_S, v_S = v_S,
                              Z_St = Z_St, Zbar = Zbar, type = type, del = del, n = n,
                              method = "Newton", control = list(maxit = 1e5, allowSingular = T))
        # if(nleqslv_res$termcd != 1 & nleqslv_res$termcd != 2){
        if(nleqslv_res$termcd != 1){
          if(max(abs(f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, n = n))) > 1e-5)
            return(.Machine$double.xmax)
        }
        w_S = f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, n = n, returnw = T)
        # if(any(is.infinite(w_S))) return(.Machine$double.xmax)
        
        
        
        if(returnw == F){
          # return(sum(G(w_S, type = type, del = del)) - n * W)
          return(sum(G(w_S, type = type, del = del)) - n * (alphaHT + 1) * log(abs(W + 1)))
        }else{
          return(w_S)
        }
      }
      
      # d_S = rep(1, sum(delta)); u_vec = 1 / pi # Reg
      # # d_S = vx[Index_S]; u_vec = -pi * vx # EL2
      # u_vec_S = u_vec[Index_S]; Uhat = mean(u_vec); 
      # Z_S = cbind(1, x_S, u_vec_S); Zbar = c(1, colMeans(x), Uhat); Z_St = t(Z_S)
      # 
      # init = rep(0, length(Zbar)); init[length(init)] = 1
      # 
      # nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, Z_S = Z_S, 
      #                       Z_St = Z_St, Zbar = Zbar, type = "SL", del = del,
      #                       method = "Newton", control = list(maxit = 1e5, allowSingular = T),
      #                       xscalm = "auto")
      # if(nleqslv_res$termcd != 1){
      #   if(max(abs(f(nleqslv_res$x, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = "SL", del = del))) > 1e-5)
      #     w_S = NA
      # }else{
      #   w_S = f(nleqslv_res$x, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = "SL", del = del, returnw = T)             
      # }
      # 
      # theta_res = c(theta_res, Reg = sum(y_S * w_S) / n) # Reg

      Omodel_d = lm(reformulate(paste0("x.", 1:pcol), response = "y"), weights = 1 / pi, data = data_S)
      # Omodel_d = lm(reformulate(paste0("x.", 1:pcol), response = "y"), data = data_S)      

      yhat_d = predict.lm(Omodel_d, data, type = "response")
      data = cbind(data, e2 = (yhat_d - y)^2); data_S = data[Index_S,]
      
      # glm(e[Index_S]^2 ~ x[Index_S,], family = gaussian(link = "log"))
      
      Vmodel_d = glm(reformulate(paste0("x.", 1:pcol), response = "e2"), family = gaussian(link = "log"),
                     weights = 1 / pi, data = data_S)
      # Vmodel_d = glm(reformulate(paste0("x.", 1:pcol), response = "e2"), family = gaussian(link = "log"),
      #                data = data_S)
      vhat = predict.glm(Vmodel_d, data, type = "response")
      
      for(type in c("EL", "ET", "SL")){
        for(DScnt in 1:2){
          if(DScnt == 1){
            v_S = rep(1, sum(delta)) # DS1; v_i = 1
          }else if(DScnt == 2){
            v_S = vhat[Index_S] # DS2; v_i = v_i
          }
          
          d_S = 1 / pi[Index_S]; 
          Z_S = cbind(1, x_S); Zbar = c(1, colMeans(x)); Z_St = t(Z_S)
          
          # init = rep(0, length(Zbar)); init[1] = -1
          # 
          # nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, v_S = v_S, Z_S = Z_S, 
          #                       Z_St = Z_St, Zbar = Zbar, type = "EL", del = del, n = n,
          #                       method = "Newton", control = list(maxit = 1e5, allowSingular = T),
          #                       xscalm = "auto")
          # if(nleqslv_res$termcd != 1){
          #   if(max(abs(f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = "EL", del = del, n = n))) > 1e-5)
          #     w_S = NA
          # }else{
          #   w_S = f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = "EL", del = del, n = n, returnw = T)             
          # }
          
          # lambda = CVXR::Variable(length(Zbar))
          # prob <- CVXR::Maximize(sum(d_S * v_S * log(-(Z_S %*% lambda) / v_S)) + sum(Zbar * lambda) * n)
          # res <- CVXR::solve(CVXR::Problem(prob))
          # lambdahat = drop(res$getValue(lambda))
          # w_S = f(lambdahat, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = "EL", del = del, n = n, returnw = T)
          
          w = CVXR::Variable(sum(Index_S))
          constraints <- list(Z_St %*% w / n == Zbar)
          if(type == "EL"){
            prob <- CVXR::Minimize(sum(v_S * (-d_S * log(w) + w) ))
          }else if(type == "ET"){
            prob <- CVXR::Minimize(sum(v_S * (-entr(w) - w * log(d_S)  - w)))
          }else if(type == "SL"){
            prob <- CVXR::Minimize(sum(v_S * ((w - d_S)^2 / d_S)))
          }
          res <- CVXR::solve(CVXR::Problem(prob, constraints), solver = "ECOS_BB")
          w_S = drop(res$getValue(w))
          
          # Z_St %*% w_S - Zbar * n
          
          gammahat = solve(Z_St %*% (Z_S / d_S / v_S),
                           Z_St %*% (y_S / d_S / v_S))
          var1 = sum((y_S - c(Z_S %*% gammahat))^2 / pi[Index_S] * (1 / pi[Index_S] - 1))
          var2 = sum((y_S - c(Z_S %*% gammahat))^2 * w_S^2 * (1 -  pi[Index_S]))
          
          theta_res = c(theta_res, setNames(sum(y_S * w_S), paste(type, "DS", DScnt, sep = ""))) # DS
          
          var_res = c(var_res, setNames(var1, paste(type, "DS", DScnt, sep = "")))
          var_res2 = c(var_res2, setNames(var2, paste(type, "DS", DScnt, sep = "")))
        }
        
        for(GECcnt in 1:2){
          if(type == "EL"){
            u_vec = -pi
          }else if(type == "ET"){
            u_vec = -log(pi)
          }else if(type == "SL"){
            u_vec = 1 / pi
          }
          if(GECcnt == 1){
            v_S = rep(1, sum(delta)); # GEC1; c_i = 1
          }else if(GECcnt == 2){
            v_S = vhat[Index_S]; u_vec = u_vec * vhat # GEC4; c_i = \hat v_i          
          }else if(GECcnt == 3){
            # v_S = vx[Index_S] / pi[Index_S]^2; u_vec = -vx / pi # GEC3; c_i = v_i / g'(d_i)
          }else if(GECcnt == 4){
            v_S = vx[Index_S]; u_vec = u_vec * vx # GEC2; c_i = v_i
          }else if(GECcnt == 5){
            # v_S = vhat[Index_S] / pi[Index_S]^2; u_vec = -vhat / pi # GEC5; c_i = \hat v_i / g'(d_i)
          }
          
          d_S = rep(1, sum(delta)); 
          u_vec_S = u_vec[Index_S]; Uhat = mean(u_vec); 
          Z_S = cbind(1, x_S, u_vec_S); Zbar = c(1, colMeans(x), Uhat); Z_St = t(Z_S)
          
          init = rep(0, length(Zbar)); init[length(init)] = 1
          
          nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, v_S = v_S, Z_S = Z_S, 
                                Z_St = Z_St, Zbar = Zbar, type = type, del = del, n = n,
                                method = "Newton", control = list(maxit = 1e5, allowSingular = T),
                                xscalm = "auto")
          if(nleqslv_res$termcd != 1){
            if(max(abs(f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, n = n))) > 1e-5)
              w_S = NA
          }else{
            w_S = f(nleqslv_res$x, d_S = d_S, v_S = v_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, n = n, returnw = T)             
          }
          
          gammahat = solve(Z_St %*% (Z_S * prime1(1 / pi[Index_S], type = type, del = del) / d_S / v_S),
                           Z_St %*% (y_S * prime1(1 / pi[Index_S], type = type, del = del) / d_S / v_S))
          var1 = sum((y_S - c(Z_S %*% gammahat))^2 / pi[Index_S] * (1 / pi[Index_S] - 1))
          var2 = sum((y_S - c(Z_S %*% gammahat))^2 * w_S^2 * (1 -  pi[Index_S]))
          
          theta_res = c(theta_res, setNames(sum(y_S * w_S), paste(type, "GEC", GECcnt, sep = ""))) # GEC
          var_res = c(var_res, setNames(var1, paste(type, "GEC", GECcnt, sep = "")))
          var_res2 = c(var_res2, setNames(var2, paste(type, "GEC", GECcnt, sep = "")))
          # var_res = c(var_res, EL0 = var(eta) / n)
        }
      }
      

      
      

        

            
      # glm(vx ~ x.1 + x.2 + x.3 + x.4, family = gaussian(link = "log"),
      #     weights = 1 / pi[Index_S], data = data_S)
      
      
      # var_res = c(var_res, EL = var(eta) / n)
      
      # nlmres= nlm(targetftn0, Zbar[length(Zbar)], d_S = d_S, Z_S = Z_S, Z_St = Z_St, 
      #             init = init, Zbar = Zbar, type = type, del = del)
      # # if(nlmres$code != 1 & nlmres$code != 2 & nlmres$code != 3) print(nlmres$estimate)
      # if(nlmres$code != 1 & nlmres$code != 2 & nlmres$code != 3) stop(nlmres$code)
      # W = nlmres$estimate
      # if(nlmres$minimum >= .Machine$double.xmax){
      #   # stop(targetftn(Zbar[length(Zbar)], d_S = d_S, Z_S = Z_S, Z_St = Z_St, init = init, Zbar = Zbar, type = type))
      #   w_S = NA
      # }else{
      #   w_S = targetftn0(W, d_S, Z_S, Z_St, init, Zbar, type, returnw = T, del = del)
      #   
      #   # nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, Z_S = Z_S, 
      #   #                       Z_St = Z_St, Zbar = c(1, colMeans(x), W), type = type, del = del,
      #   #                       method = "Newton", control = list(maxit = 1e5, allowSingular = T),
      #   #                       xscalm = "auto")
      #   # w_S = f(nleqslv_res$x, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, returnw = T)   
      # }
      # theta_res = c(theta_res, Qin = sum(y_S * w_S) / n) # Qin
      
      # var_res = c(var_res, GEC = var(eta) / n)
      
      CR_res = ifelse(abs(theta_res - theta) < qnorm(0.975) * sqrt(var_res), 1, 0)
      CR_res2 = ifelse(abs(theta_res - theta) < qnorm(0.975) * sqrt(var_res2), 1, 0)
      
      # alpha_vec = c(alpha_vec, nlmres$par)
      theta_mat = cbind(theta_mat, theta_res)
      var_mat = cbind(var_mat, var_res); var_mat2 = cbind(var_mat2, var_res2)
      CR_mat = cbind(CR_mat, CR_res); CR_mat2 = cbind(CR_mat2, CR_res2)
    }
    # list(theta_mat, alpha = alpha_vec, var_mat, CR_mat)
    list(theta_mat, var_mat, var_mat2, CR_mat, CR_mat2)
  }
# final_res
paste("# of failure:", sum(!sapply(lapply(final_res, function(x) x[[1]]), function(x) is.numeric(unlist(x))))); 
final_res0 = lapply(final_res, function(x) x[[1]])
print(final_res0[!sapply(final_res0, function(x) is.numeric(unlist(x)))])

final_res = final_res[sapply(final_res, length) == max(sapply(final_res, length))]

if(!interactive()) save.image(paste(timenow0, ".RData", sep = ""))

final_res1 = lapply(final_res, function(x) x[[1]])

stopCluster(cl)
timenow2 = Sys.time()
print(timenow2 - timenow1)
# if(sum(!sapply(final_res1, function(x) is.numeric(unlist(x)))) != 0) stop(paste(final_res1))

final_res1 = final_res1[sapply(final_res1, function(x) is.numeric(unlist(x)))]

# res1 <- Reduce(`+`, final_res1) / length(final_res1)
# colnames(res1) = rnames

tmpres = do.call(rbind, lapply(final_res1, c))

BIAS = colMeans(tmpres - theta, na.rm = TRUE)
SE = sqrt(colMeans(apply(tmpres, 2, function(x) (x - mean(x, na.rm = TRUE))^2), na.rm = TRUE))
RMSE = sqrt(colMeans((tmpres - theta)^2, na.rm = TRUE))

res <- cbind(BIAS, SE, RMSE)
rownames(res) = row.names(final_res1[[1]])

colSums(is.na(tmpres))

# xtable::xtable(res, digits = 4)

res

final_res2 = lapply(final_res, function(x) x[[2]])
final_res2 = final_res2[sapply(final_res2, function(x) is.numeric(unlist(x)))]
tmpres2 = do.call(rbind, lapply(final_res2, c))

final_res4 = lapply(final_res, function(x) x[[4]])
final_res4 = final_res4[sapply(final_res4, function(x) is.numeric(unlist(x)))]
tmpres4 = do.call(rbind, lapply(final_res4, c))
cbind()

final_res3 = lapply(final_res, function(x) x[[3]])
final_res3 = final_res3[sapply(final_res3, function(x) is.numeric(unlist(x)))]
tmpres3 = do.call(rbind, lapply(final_res3, c))

final_res5 = lapply(final_res, function(x) x[[5]])
final_res5 = final_res5[sapply(final_res5, function(x) is.numeric(unlist(x)))]
tmpres5 = do.call(rbind, lapply(final_res5, c))
res2 = cbind(RB1 = (colMeans(tmpres2, na.rm = TRUE) - SE^2) / SE^2, CR1 = colMeans(tmpres4, na.rm = TRUE),
      RB2 = (colMeans(tmpres3, na.rm = TRUE) - SE^2) / SE^2, CR2 = colMeans(tmpres5, na.rm = TRUE))
rownames(res2) = row.names(final_res1[[1]])

xtable::xtable(cbind(res, res2), digits = 4)
# res2
