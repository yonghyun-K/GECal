# Load the survey package
library(survey)

setwd("C:/Users/User/Box/Data/KHP2.2_sas")

a_hh = read.csv("a_hh.csv")

a_hh$REGION1 <- as.factor(a_hh$REGION1)
a_hh$REGION2 <- as.factor(a_hh$REGION2)

a_hh <- a_hh[!is.na(a_hh$REGION1), ] # Maybe it has to be removed

# table(a_hh$REGION1, a_hh$REGION2, useNA ="always")

# Assuming your data frame is named `a_hh`
# Define the survey design object
survey_design <- svydesign(
  id = ~1,  # Assuming no clustering (if there's a cluster variable, replace ~1 with ~cluster_variable)
  strata = ~REGION1 + REGION2,  # Specify the strata variables
  weights = ~H_WGC,  # Specify the weight variable
  data = a_hh  # Specify the dataset
)

# Calculate the mean for the variable OTC_MED
survey_means <- svymean(~OTC_MED, design = survey_design)

# Display the results
print(survey_means)

survey_freq <- svytable(~HEXP1, design = survey_design)

survey_freq

# Compute the weighted totals and their standard errors for HEXP1
weighted_totals <- svytotal(~as.factor(HEXP1), design = survey_design, vartype = c("se"))

a_ind = read.csv("a_ind.csv")

library(dplyr)
a_ind <- left_join(a_ind, a_hh %>% select(HHID, REGION2), by = join_by(HHID))

a_ind <- subset(a_ind, I_INC1 != 0) # To be removed
a_ind <- a_ind[!is.na(a_ind$REGION1), ] # To be removed

tmptable = tapply(a_ind$I_WGC, a_ind %>% select(REGION1), sum) / tapply(a_ind$I_WGC, a_ind %>% select(REGION1), length)
# table(a_ind %>% select(REGION1, REGION2), useNA = "always")

# Define the survey design object
survey_design <- svydesign(
  id = ~1,  # Assuming no clustering
  # strata = ~REGION1,
  # strata = ~REGION1 + REGION2,  # Strata variables
  weights = ~I_WGC,  # Weight variable
  data = a_ind
  # data = subset(a_ind, I_INC1 != 0)  # Filtered dataset
)

(weighted_mean <- sum(a_ind$I_INC1 * a_ind$I_WGC) / sum(a_ind$I_WGC))

sqrt(var(a_ind$I_INC1) / nrow(a_ind) * sum(a_ind$I_WGC)^2)

# Step 2: Add a column for strata
a_ind$strata <- interaction(a_ind$REGION1)

# Step 3: Split data by strata
library(dplyr)

# Step 3: Split data by strata and calculate variance
variance <- a_ind %>%
  group_by(strata) %>%
  reframe(
    n_h = n(),
    N_h = sum(I_WGC),
    stratum_var = (N_h / n_h)^2 * n_h * var(I_INC1)
    # w_delta = I_WGC * (I_INC1 - weighted_mean),
    # mean_w_delta = mean(I_WGC * (I_INC1 - weighted_mean)),
    # stratum_var = sum((w_delta - mean_w_delta)^2) / (n_h - 1)
  ) %>%
  summarise(total_variance = sum(stratum_var)) %>%
  pull(total_variance)

sum(a_ind$I_WGC)
summary(a_ind$I_WGC)

# Step 4: Calculate the standard error
(standard_error <- sqrt(variance))
log(standard_error)

svymean(~I_INC1, design = survey_design, vartype = c("se"))

tmp <- svytotal(~I_INC1, design = survey_design, vartype = c("se"))
log(sqrt(attr(tmp,"var")))

svytotal(~I_INC1, design = survey_design, vartype = c("se"))


svytotal(~intercept, design = survey_design, vartype = c("se"))
svymean(~I_INC1, design = survey_design, vartype = c("se"))[1] * sum(subset(a_ind, I_INC1 != 0) %>% select(I_WGC))
 
svyvar(~I_INC1, survey_design)

calibration <- GECal::GEcalib(~ 0, dweight = I_WGC, data = subset(a_ind, I_INC1 != 0),
                              const = numeric(0),
                              entropy = "SL", method = "DS")
GECal::estimate(I_INC1 ~ 1, calibration = calibration, data = subset(a_ind, I_INC1 != 0))$estimate # HT estimator

calibration <- GECal::GEcalib(~ 1, dweight = I_WGC, data = subset(a_ind, I_INC1 != 0),
                              const = c(sum(subset(a_ind, I_INC1 != 0) %>% select(I_WGC))),
                              entropy = "SL", method = "DS")
GECal::estimate(I_INC1 ~ 1, calibration = calibration, data = subset(a_ind, I_INC1 != 0))$estimate # HT estimator


sqrt(var(unlist(subset(a_ind, I_INC1 != 0) %>% select(I_INC1))) * 
       (1 - nrow(subset(a_ind, I_INC1 != 0)) / sum(subset(a_ind, I_INC1 != 0) %>% select(I_WGC))) / nrow(subset(a_ind, I_INC1 != 0))) *
  sum(subset(a_ind, I_INC1 != 0) %>% select(I_WGC))

# Compute the mean for I_INC1
mean_income <- svymean(~I_INC1, design = survey_design, vartype = c("se"))

# Display the results
print(mean_income)

# table(a_ind$HS1, useNA = "always")

table(a_ind$REGION1, a_ind$REGION2)

sum(is.na(a_ind$I_WGC))
a_ind$I_WGC[a_ind$REGION1[!is.na(a_ind$REGION1)] == 11]

survey_design <- svydesign(
  id = ~1,  # Assuming no clustering
  strata = ~REGION1 + REGION2,  # Strata variables
  weights = ~I_WGC,  # Weight variable
  data = a_ind[!is.na(a_ind$REGION1) & !is.na(a_ind$HS1), ]
  # data = a_ind[!is.na(a_ind$REGION1), ]
)


# Compute weighted frequencies for HS1
weighted_frequencies <- svytotal(~as.factor(HS1), design = survey_design,
                                 na.rm=T)
weighted_frequencies
