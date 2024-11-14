# Model-assisted calibration using
# Generalized entropy calibration in survey sampling

# Simulation setup

if (!interactive()) {
  args <- as.numeric(commandArgs(trailingOnly = TRUE))
} else{
  args <- c(50)
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


names(nhis)[names(nhis) == "Height"] <- "HT"
names(nhis)[names(nhis) == "Weight"] <- "WT"
names(nhis)[names(nhis) == "ProvCode"] <- "REGION1"
names(nhis)[names(nhis) == "Gender"] <- "SEX"
nhis$REGION1 <- as.factor(nhis$REGION1)
nhis$SEX <- as.factor(nhis$SEX)
nhis$AgeGroup <- as.factor(nhis$AgeGroup)
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

tab1 = round(table(nhis$AgeGroup, nhis$REGION1, nhis$SEX) / N * smp.size)
tab1 = ifelse(tab1 < 1, 1, tab1)
n = sum(tab1)

theta = mean(nhis$Hemo)

final_res <- foreach(
  simnum = 1:SIMNUM,
  .packages = c("nleqslv", "PracTools", "GECal"),
  .errorhandling = "pass"
) %dopar% {
  
  index = c()
  pi_S = c()
  for(AgeGroup in unique(nhis$AgeGroup)){
    for(REGION1 in unique(nhis$REGION1)){
      for(SEX in unique(nhis$SEX)){
        idx_tmp = which(nhis$AgeGroup == AgeGroup &
                          nhis$REGION1 == REGION1 & 
                          nhis$SEX == SEX)
        n_h = tab1[AgeGroup,REGION1,SEX]
        index = c(index, sample(idx_tmp, n_h, replace = FALSE))
        pi_S = c(pi_S, rep(n_h / length(idx_tmp), n_h))
      }
    }
  }
  
  delta = as.integer(1:N %in% index)
  Rmodel = glm(delta ~ AgeGroup + REGION1 + SEX, family = binomial,
               data = nhis)
  pihat = predict.glm(Rmodel, nhis, type = "response")
  
  Omodel = lm(Hemo ~ AgeGroup + REGION1 + SEX + HT + WT, data = nhis)
  yhat = predict.lm(Omodel, nhis, type = "response")
  
  nhis.samp <- nhis[index, ]
  # d_S       <- 1 / pi_S
  d_S       <- 1 / pihat[index]
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
  
  #HT estimator
  const = numeric(0)
  calibration <- GEcalib(
    ~ 0,
    dweight = d_S,
    data = nhis.samp,
    const = const,
    entropy = 1,
    method = "DS"
  )
  res_est = estimate(y_S ~ 1, calibration = calibration)$estimate

  theta_res = c(theta_res, setNames(res_est[1] / N, "IPW"))
  se_res = c(se_res, setNames(res_est[2] / N, "IPW"))
  
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
  theta_res = c(theta_res, SIPW = sum((y_S) / pihat[index]) / sum(1 / pihat[index])) # AIPW  
  hhat = pihat * model.matrix(~1, data = nhis)
  kappa = solve(t(hhat[index,]) %*% (hhat[index,] * (1 / pihat[index] - 1) / pihat[index]),
                t(hhat[index,]) %*% ((y_S - yhat[index]) * (1 / pihat[index] - 1) / pihat[index]))
  eta = yhat + drop(hhat %*% kappa)
  eta[index] = eta[index] + (y_S - yhat[index] - drop(hhat %*% kappa)[index]) / pihat[index]
  
  se_res = c(se_res, SIPW = sqrt(var(eta) / N))

  # AIPW estimator
  theta_res = c(theta_res, AIPW = (sum(yhat) + sum((y_S - yhat[index]) / pihat[index])) / N) # AIPW
  
  hhat = pihat * model.matrix(Omodel)[, -1]
  kappa = solve(t(hhat[index,]) %*% (hhat[index,] * (1 / pihat[index] - 1) / pihat[index]),
                t(hhat[index,]) %*% ((y_S - yhat[index]) * (1 / pihat[index] - 1) / pihat[index]))
  eta = yhat + drop(hhat %*% kappa)
  eta[index] = eta[index] + (y_S - yhat[index] - drop(hhat %*% kappa)[index]) / pihat[index]
  
  se_res = c(se_res, AIPW = sqrt(var(eta) / N))

  vectmp <- c("AgeGroup","REGION1", "SEX", "HT", "WT")
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
    
    res_est = estimate(y_S ~ 1, calibration = calibration)$estimate
    
    theta_res = c(theta_res, setNames(res_est[1] / N, paste("GEC", entropy, sep = "_"))) # DS
    se_res = c(se_res, setNames(res_est[2] / N, paste("GEC", entropy, sep = "_")))
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

xtable(res * 1e2)

colSums(is.na(tmpres))


final_res2 = lapply(final_res, function(x) x[[2]])
final_res2 = final_res2[sapply(final_res2, function(x) is.numeric(unlist(x)))]
tmpres2 = do.call(rbind, lapply(final_res2, c))

final_res3 = lapply(final_res, function(x) x[[3]])
final_res3 = final_res3[sapply(final_res3, function(x) is.numeric(unlist(x)))]
tmpres3 = do.call(rbind, lapply(final_res3, c))

res2 = cbind(RB = (colMeans(tmpres2^2, na.rm = TRUE) - SE^2) / SE^2, 
             CR = colMeans(tmpres3, na.rm = TRUE))

# hist(g(1 / pi, entropy = -1))
# hist(g(1 / pi, entropy = 0))
# hist(1 / pi)
# summary(1 / pi)

# summary(lm(EXPTOTAL ~ . + I(1/pi), data = smho98))
# summary(lm(smho98$EXPTOTAL ~ I(1/pi)))
# summary(lm(EXPTOTAL ~ I(1/pi) + Y_IP, data = smho98))
