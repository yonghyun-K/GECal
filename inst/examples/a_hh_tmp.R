# Load the survey package
library(survey)

setwd("C:/Users/User/Box/Data/KHP2.2_sas")

a_ind = read.csv("a_ind.csv")
a_ind$EDU_STAT[is.na(a_ind$EDU_STAT)] <- 5

a_ind = a_ind[,sapply(a_ind, function(x) sum(is.na(x)) <= 4994)]

a_ind = a_ind[complete.cases(a_ind),]

# Specify the columns to convert to factors
cols_to_factor <- c("REGION1", "PRE_RES", "DEATH_I_YN", "SEX", "MARR", "EDU", 
                    "EDU_STAT", "LIVE_T_YN", "HEALTH_INS", "DISA_YN")

# Convert specified columns to factors
a_ind[cols_to_factor] <- lapply(a_ind[cols_to_factor], as.factor)


which(names(a_ind) == "P1")
which(names(a_ind) == "HS8_YN")
which(names(a_ind) == "CD1_HTN")

a_ind = a_ind[,c(1:10, 23, 56:85)]

summary(a_ind)
summary(a_ind[,155:253])
mice::mice



a_hh = read.csv("a_hh.csv")

library(dplyr)
a_ind <- left_join(a_ind, a_hh %>% select(HHID, REGION2), by = join_by(HHID))

# subset(a_ind, select = -I_INC1)

# names(a_ind)
# a_ind[, !grepl("^CARE", names(a_ind))]

head(a_ind$CD1_HTN)

names(a_ind[, !grepl("^CARE|CD", names(a_ind))])
a_ind[-I_INC1]

library(missForest)
a_ind <- missForest(a_ind)$ximp
a_ind_imp <- a_ind[, !names(a_ind) %in% c("N_LT_REASON", "DISA_TY")][,1:18]
a_ind_imp <- cbind(a_ind_imp, CD1_HTN = a_ind$CD1_HTN)
a_ind_imp <- a_ind_imp[!is.na(a_ind_imp$CD1_HTN),]
a_ind_imp <- a_ind_imp[!is.na(a_ind_imp$SEX),]
summary(a_ind_imp)
# a_ind_imp <- missForest(a_ind_imp)$ximp


summary(lm(CD1_HTN ~ ., data = a_ind))



library(glmnet)

lasso_model <- cv.glmnet(as.matrix(subset(a_ind, select = -CD1_HTN)), 
                         a_ind$CD1_HTN)

lasso_model <- cv.glmnet(as.matrix(a_ind[, !grepl("^CD", names(a_ind))]), 
                         a_ind$CD1_HTN)

lasso_model <- glmnet(as.matrix(a_ind[, !grepl("^CD", names(a_ind))]), 
                         a_ind$CD1_HTN, lambda = 0.02)


lasso_model <- glmnet(as.matrix(subset(a_ind, select = -CD1_HTN)), 
        a_ind$CD1_HTN, lambda = 0.01)

selected_variables <- coef(lasso_model, s = "lambda.min")
selected_variables
