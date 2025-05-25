# Model-assisted calibration using
# Generalized entropy calibration in survey sampling

# Simulation setup

if (!interactive()) {
  args <- as.numeric(commandArgs(trailingOnly = TRUE))
} else{
  args <- c(5)
}

timenow1 = Sys.time()
timenow0 = gsub(' ', '_', gsub('[-:]', '', timenow1))
timenow = paste(timenow0, ".txt", sep = "")

# install.packages( "MatrixModels", type="win.binary" )

library(nleqslv)
library(PracTools)
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
# library(caret)
library(tidyverse)
library(xtable)
library(GECal)
# library(mice)
# library(kableExtra)

set.seed(11)
SIMNUM = args[1]

if (!interactive()) {
  dir.create(timenow0)
  setwd(timenow0)
  
  sink(timenow, append = TRUE)
}
# setup parallel backend to use many processors
cores = min(detectCores() - 3, 101)
print(paste("cores =", cores))
# cl <- makeCluster(cores, outfile = timenow) #not to overload your computer
cl <- makeCluster(cores)
registerDoParallel(cl)
# summary(smho98)
# summary(smho.N874)

# nhis = read.csv("C:/Users/ghkfk/Box/Data/NHIS2021/NHIS_2021.CSV")
nhis = read.csv("C:/Users/User/Box/Data/NHIS2021/NHIS_2021.CSV")
# nhis = read.csv("/Users/yhkwon/Library/CloudStorage/Box-Box/Data/NHIS2021/NHIS_2021.CSV")

names(nhis)[names(nhis) == "Height"] <- "HT"
names(nhis)[names(nhis) == "Weight"] <- "WT"
names(nhis)[names(nhis) == "ProvCode"] <- "REGION1"
names(nhis)[names(nhis) == "Gender"] <- "SEX"
nhis$REGION1 <- as.factor(nhis$REGION1)
nhis$SEX <- as.factor(nhis$SEX)
nhis$AgeGroup <- as.factor(nhis$AgeGroup)

nhis = nhis[sample(1:nrow(nhis), 1e5), ]
N <- nrow(nhis)

nhis[c("LDL", "HDL", "CholTotal", "Triglyc", 
       "DentalCaries", "Calculus", "UrineProt")] <- NULL
# apply(nhis, 2, function(x) sum(is.na(x)))
# imp <- mice(nhis, 1, method = "sample", maxit = 1)
# nhis <- complete(imp)

# Custom function to impute missing values by sampling observed values
impute_by_sampling <- function(data) {
  data[] <- lapply(data, function(x) {
    if (any(is.na(x))) {
      x[is.na(x)] <- sample(x[!is.na(x)], sum(is.na(x)), replace = TRUE)
    }
    x
  })
  return(data)
}

# Apply the function to your dataset
nhis <- impute_by_sampling(nhis)

smp.size <- 1000

tab1 = table(nhis$AgeGroup, nhis$REGION1, nhis$SEX)
tab1 = ifelse(tab1 > 15, 5, round(tab1 / 3))

# tab1 = round(table(nhis$AgeGroup, nhis$REGION1, nhis$SEX) / N * smp.size)
# tab1 = ifelse(tab1 < 1, 1, tab1)
# tab1 = ifelse(table(nhis$AgeGroup, nhis$REGION1, nhis$SEX) == 0, 0, tab1)

n = sum(tab1)

nhis$Smoking <- ifelse(nhis$Smoking == 3, 1, 0)

# summary(nhis)

# nhis$Hemo <- nhis$OralExam  # Variable of interest y is Smoking

theta = mean(nhis$Hemo)

final_res <- foreach(
  simnum = 1:SIMNUM,
  .packages = c("nleqslv", "PracTools", "GECal"),
  .errorhandling = "pass"
) %dopar% {
  
  index = c()
  pi_S = c()
  pi = rep(0, N)
  for(AgeGroup in unique(nhis$AgeGroup)){
    for(REGION1 in unique(nhis$REGION1)){
      for(SEX in unique(nhis$SEX)){
        idx_tmp = which(nhis$AgeGroup == AgeGroup &
                          nhis$REGION1 == REGION1 & 
                          nhis$SEX == SEX)
        n_h = tab1[AgeGroup,REGION1,SEX]
        index = c(index, sample(idx_tmp, n_h, replace = FALSE))
        pi_S = c(pi_S, rep(n_h / length(idx_tmp), n_h))
        pi[idx_tmp] = n_h / length(idx_tmp)
      }
    }
  }
  # summary(pi)
  
  delta = as.integer(1:N %in% index)
  Rmodel = glm(delta ~ AgeGroup + REGION1 + SEX, family = binomial,
               data = nhis)
  
  # Rmodel = glm(delta ~ AgeGroup * REGION1 * SEX, family = binomial,
  #              data = nhis)
  
  # Rmodel = glm(delta ~ AgeGroup + REGION1 + SEX + HT + WT + Waist + Alcohol + SysBP + DiaBP+ FBS + Creatinine, family = binomial,
  #              data = nhis)
  
  # Rmodel = glm(delta ~ AgeGroup + SEX, family = binomial,
  #              data = nhis)
  pihat0 = predict.glm(Rmodel, nhis, type = "response")
  
  
  findphi2 = function(phi, x0, Z, delta, ..., returnw = F){
    pi_phi = drop(1 / (1 + exp(-x0 %*% phi)))
    w_phi = ifelse(delta, (1 - pi_phi) / pi_phi, -1)
    if(returnw){
      return(drop(1 / pi_phi))
    }else{
      return(drop(w_phi %*% Z))
    }
  }
  jacphi2 = function(phi, x0, Z, delta, ..., returnw = F){
    pi_phi = drop(1 / (1 + exp(-x0 %*% phi)))
    logit_phi = pi_phi / (1 - pi_phi)
    return(-t(Z) %*% (x0 * ifelse(delta, 1 / logit_phi, 0) ))
    # return(-t(x0) %*% (Z * ifelse(delta, 1 / logit_phi, logit_phi) ))
  }
  
  xR = model.matrix(Rmodel)
  nleqslv_res = nleqslv(Rmodel$coefficients, findphi2, jac = jacphi2, x0 = xR, Z = xR, delta = delta,
                         control = list(maxit = 1e5, allowSingular = T), xscalm = "auto",
                         method = "Newton")
  
  if(nleqslv_res$termcd != 1 & max(abs(findphi2(nleqslv_res$x, x0 = xR, Z = xR, delta = delta))) > 1e-5){
    w_S = NA
  }else{
    w_S = findphi2(nleqslv_res$x, x0 = xR, Z = xR, delta = delta, returnw = T)
  }
  # drop(t(xR[Index_S,]) %*% w_S2[Index_S]); colSums(xR)
  w_S2 = w_S
  
  # summary(lm(Smoking ~ ., data = nhis))

  # Omodel = lm(Hemo ~ AgeGroup + SEX + HT + WT, data = nhis)
    
  Omodel = lm(Hemo ~ AgeGroup + SEX + HT + WT + Waist + Alcohol + SysBP + DiaBP + FBS + Creatinine, data = nhis)
  yhat = predict.lm(Omodel, nhis, type = "response")
  
  nhis.samp <- nhis[index, ]
  # d_S       <- 1 / pi_S
  y_S <- nhis.samp$Hemo
  
  # omega2_S = matrix(0, nrow = length(index), ncol = length(index))
  # for (cnt in 1:5) {
  #   idx_strat2 = which(smho98.samp$stratum5 == cnt)
  #   # omega2_S[idx_strat2, idx_strat2] <-
  #   diag(omega2_S)[idx_strat2] <-
  #     1 / pi_S[idx_strat2] * (1 / pi_S[idx_strat2] - 1)
  # }
  
  theta_res = NULL
  se_res = NULL
  
  for(pimethod in 0:3){
    if(pimethod == 0){
      pihat = pi
    }else if(pimethod == 1){
      pihat = rep(n / N, N)
    }else if(pimethod == 2){
      pihat = pihat0
    }else if(pimethod == 3){
      pihat = 1 / w_S2
    }
    pihat = ifelse(pihat > 0.65, 0.65, pihat) # To make CE convergent
    d_S <- 1 / pihat[index]
  
  # #HT estimator
  # const = numeric(0)
  # calibration <- GEcalib(
  #   ~ 0,
  #   dweight = d_S,
  #   data = nhis.samp,
  #   const = const,
  #   entropy = 1,
  #   method = "DS"
  # )
  # res_est = estimate(y_S ~ 1, calibration = calibration)$estimate
  # 
  # theta_res = c(theta_res, setNames(res_est[1] / N, "IPW"))
  # se_res = c(se_res, setNames(res_est[2] / N, "IPW"))
  
  # Hajek estimator
  # const = N
  # calibration <- GEcalib(
  #   ~ 1,
  #   dweight = d_S,
  #   data = nhis.samp,
  #   const = const,
  #   entropy = 1,
  #   method = "DS"
  # )
  # res_est = estimate(y_S ~ 1, calibration = calibration)$estimate
  # 
  # theta_res = c(theta_res, setNames(res_est[1] / N, "SIPW"))
  # se_res = c(se_res, setNames(res_est[2] / N, "SIPW"))

  # SIPW estimator
  theta_res = c(theta_res, IPW = sum((y_S) / pihat[index]) / sum(1 / pihat[index])) # AIPW  
  hhat = pihat * model.matrix(~1, data = nhis)
  kappa = solve(t(hhat[index,]) %*% (hhat[index,] * (1 / pihat[index] - 1) / pihat[index]),
                t(hhat[index,]) %*% ((y_S - yhat[index]) * (1 / pihat[index] - 1) / pihat[index]))
  eta = yhat + drop(hhat %*% kappa)
  eta[index] = eta[index] + (y_S - yhat[index] - drop(hhat %*% kappa)[index]) / pihat[index]
  
  se_res = c(se_res, IPW = sqrt(var(eta) / N))

  # AIPW estimator
  # theta_res = c(theta_res, AIPW = (sum(yhat) + sum((y_S - yhat[index]) / pihat[index])) / N) # AIPW
  # 
  # if(pimethod == 1){
  #   hhat = pihat * model.matrix(Rmodel)
  # }else if(pimethod == 2){
  #   hhat = pihat / (1 - pihat) * model.matrix(Rmodel)
  # }
  # # hhat = pihat * model.matrix(Rmodel)
  # kappa = solve(t(hhat[index,]) %*% (hhat[index,] * (1 / pihat[index] - 1) / pihat[index]),
  #               t(hhat[index,]) %*% ((y_S - yhat[index]) * (1 / pihat[index] - 1) / pihat[index]))
  # eta = yhat + drop(hhat %*% kappa)
  # eta[index] = eta[index] + (y_S - yhat[index] - drop(hhat %*% kappa)[index]) / pihat[index]
  # 
  # se_res = c(se_res, AIPW = sqrt(var(eta) / N))

  vectmp <- as.character(formula(Omodel))[3]
  # vectmp <- c("Y_IP")
  
  fortmp <- formula(paste("~", paste(vectmp, collapse = "+")))
  fortmp2 <- formula(paste("~", paste(c(vectmp, "g(d_S)"), collapse = "+")))
  
  # Omodel_d = lm(reformulate(vectmp, response = "y_S"), weights = 1 / pi_S, data = smho98.samp)
  # # Omodel_d = lm(reformulate(paste0("x.", 1:pcol), response = "y"), data = data_S)      
  # 
  # yhat_d = predict.lm(Omodel_d, smho98.samp, type = "response")
  # e2 = (yhat_d - y_S)^2
  # Vmodel_d = glm(reformulate(vectmp, response = "e2"),
  #                weights = 1 / pi_S, data = smho98.samp)
  # vhat = predict.glm(Vmodel_d, smho98, type = "response")
  
  for (entropy in list(-1, -1/2, 1, "CE")) {
    # const = colSums(model.matrix(fortmp, nhis))
    # 
    # calibration <- GEcalib(
    #   fortmp,
    #   dweight = d_S,
    #   data = nhis.samp,
    #   const = const,
    #   entropy = entropy,
    #   method = "DS"
    # )
    # 
    # res_est = estimate(y_S ~ 1, calibration = calibration)$estimate
    # 
    # theta_res = c(theta_res, setNames(res_est[1] / N, paste("DS", entropy, sep = "_"))) # DS
    # se_res = c(se_res, setNames(res_est[2] / N, paste("DS", entropy, sep = "_")))
    
    if(pimethod == 1){
      const = colSums(model.matrix(fortmp, nhis))
      calibration <- GEcalib(
        fortmp,
        dweight = d_S,
        data = nhis.samp,
        const = const,
        entropy = entropy,
        method = "GEC0"
      )
    }else{
      const = colSums(cbind(model.matrix(fortmp, nhis), 
                            g(1 / pihat, entropy = entropy)))
      calibration <- GEcalib(
        fortmp2,
        dweight = d_S,
        data = nhis.samp,
        const = const,
        entropy = entropy,
        method = "GEC"
      )
    }
    
    res_est = estimate(y_S ~ 1, calibration = calibration)$estimate
    
    theta_res = c(theta_res, setNames(res_est[1] / N, paste("GEC", entropy, sep = "_"))) # DS
    se_res = c(se_res, setNames(res_est[2] / N, paste("GEC", entropy, sep = "_")))
  }
  }
  
  CR_res = ifelse(abs(theta_res - theta) < qnorm(0.975) * se_res, 1, 0)
  list(theta_res, se_res, CR_res)
}

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
rownames(res) = names(final_res1[[1]])
res

xtable(res * 1e2, digits = 3)

colSums(is.na(tmpres))


final_res2 = lapply(final_res, function(x) x[[2]])
final_res2 = final_res2[sapply(final_res2, function(x) is.numeric(unlist(x)))]
tmpres2 = do.call(rbind, lapply(final_res2, c))

final_res3 = lapply(final_res, function(x) x[[3]])
final_res3 = final_res3[sapply(final_res3, function(x) is.numeric(unlist(x)))]
tmpres3 = do.call(rbind, lapply(final_res3, c))

res2 = cbind(RB = (colMeans(tmpres2^2, na.rm = TRUE) - SE^2) / SE^2, 
             CR = colMeans(tmpres3, na.rm = TRUE))

# xtable(res2 * 1e2)

# hist(g(1 / pi, entropy = -1))
# hist(g(1 / pi, entropy = 0))
# hist(1 / pi)
# summary(1 / pi)

# summary(lm(EXPTOTAL ~ . + I(1/pi), data = smho98))
# summary(lm(smho98$EXPTOTAL ~ I(1/pi)))
# summary(lm(EXPTOTAL ~ I(1/pi) + Y_IP, data = smho98))
