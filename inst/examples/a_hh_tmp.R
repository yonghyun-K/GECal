# Load the survey package
library(survey)
library(dplyr)
library(GECal)

setwd("C:/Users/ghkfk/Box/Data/KHP2.2_sas")

a_ind = read.csv("a_ind.csv")
a_hh = read.csv("a_hh.csv")
a_ind$EDU_STAT[is.na(a_ind$EDU_STAT)] <- 5

a_ind <- left_join(a_ind, a_hh %>% select(HHID, REGION2), by = join_by(HHID))

a_ind = a_ind[,sapply(a_ind, function(x) sum(is.na(x)) <= 4994)]

a_ind = a_ind[complete.cases(a_ind),]

which(names(a_ind) == "REGION1")
which(names(a_ind) == "REGION2")
which(names(a_ind) == "LIVE_T_YN")

which(names(a_ind) == "P1")
which(names(a_ind) == "HS8_YN")
which(names(a_ind) == "CD1_HTN")

a_ind = a_ind[,c(1:16, 23, 56:85, 96)]

# Specify the columns to convert to factors
cols_to_factor <- c("REGION1", "PRE_RES", "DEATH_I_YN", "SEX", "MARR", "EDU", 
                    "EDU_STAT", "LIVE_T_YN")

# Convert specified columns to factors
a_ind[cols_to_factor] <- lapply(a_ind[cols_to_factor], as.factor)

which(names(a_ind) == "MA1")
which(names(a_ind) == "HS6_YN")

a_ind[28:41] <- lapply(a_ind[28:41], as.factor)

# summary(a_ind)


# library(missForest)
# a_ind <- missForest(a_ind)$ximp
# a_ind_imp <- a_ind[, !names(a_ind) %in% c("N_LT_REASON", "DISA_TY")][,1:18]
# a_ind_imp <- cbind(a_ind_imp, CD1_HTN = a_ind$CD1_HTN)
# a_ind_imp <- a_ind_imp[!is.na(a_ind_imp$CD1_HTN),]
# a_ind_imp <- a_ind_imp[!is.na(a_ind_imp$SEX),]
# summary(a_ind_imp)
# # a_ind_imp <- missForest(a_ind_imp)$ximp

# Variable selection
library(glmnet)

lasso_model <- cv.glmnet(as.matrix(subset(a_ind, select = -CD1_HTN)), 
                         a_ind$CD1_HTN)

lasso_model <- cv.glmnet(as.matrix(a_ind[, !grepl("^CD", names(a_ind))]), 
                         a_ind$CD1_HTN)

lasso_model <- glmnet(as.matrix(a_ind[, !grepl("^CD", names(a_ind))]), 
                         a_ind$CD1_HTN, lambda = 0.02)

lasso_model <- glmnet(as.matrix(a_ind[, !grepl("^CD", names(a_ind))]), 
                      a_ind$CD1_HTN, lambda = 0.01)


lasso_model <- glmnet(as.matrix(subset(a_ind, select = -CD1_HTN)), 
        a_ind$CD1_HTN, lambda = 0.01)

lasso_model <- glmnet(as.matrix(a_ind[, !grepl("^CD", names(a_ind))]), 
                      a_ind$CD1_HTN, lambda = 0.02)
(selected_variables <- coef(lasso_model, s = "lambda.min"))
# Convert sparse matrix to a regular vector
selected_variables <- as.matrix(selected_variables)

# Identify variables with non-zero coefficients
selected_variables <- rownames(selected_variables)[selected_variables[, 1] != 0]
selected_variables <- selected_variables[-1] # Remove the intercept

# Subset data for selected variables
selected_data <- a_ind[, selected_variables]

# Perform regression using lm
lm_model <- lm(a_ind$CD1_HTN ~ ., data = selected_data)

summary(lm_model)

# plot(lm_model)


# Define the survey design object
survey_design <- svydesign(
  id = ~1,  # Assuming no clustering
  # cluster = ~PIDWON,
  strata = ~REGION1 + REGION2,  # Strata variables
  weights = ~I_WGC,  # Weight variable
  data = a_ind
)


# svymean(~I_INC1, design = survey_design, vartype = c("se"))
# svytotal(~I_INC1, design = survey_design, vartype = c("se"))

# Compute the mean for I_INC1
mean_income <- svymean(~CD1_HTN, design = survey_design, vartype = c("se"))

# Display the results
print(mean_income)


formulatmp = ~BIRTH_Y + HB11_1 + WT + 
  HS_EQ1 + HS_EQ2 + HS_EQ3 + HS_SRH
formulatmp2 = update(formulatmp, ~ . + g(I_WGC))
const = a_ind$I_WGC %*% model.matrix(formulatmp, data = a_ind)
const2 = c(const, -nrow(a_ind))

calibration = GECal::GEcalib(formulatmp, dweight = I_WGC, data = a_ind, const = const,
                             entropy = "SL", method = "GEC0")

calibration = GECal::GEcalib(formulatmp2, dweight = I_WGC, data = a_ind, const = const2,
                             entropy = "EL", method = "GEC")

plot(calibration$w, a_ind$I_WGC)

