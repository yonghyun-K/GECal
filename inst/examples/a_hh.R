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

# a_ind <- subset(a_ind, I_INC1 != 0) # To be removed
# a_ind <- a_ind[!is.na(a_ind$REGION1), ] # To be removed

# tmptable = tapply(a_ind$I_WGC, a_ind %>% select(REGION1), sum) / tapply(a_ind$I_WGC, a_ind %>% select(REGION1), length)
# table(a_ind %>% select(REGION1, REGION2), useNA = "always")

# Define the survey design object
survey_design <- svydesign(
  id = ~1,  # Assuming no clustering
  # strata = ~REGION1,
  # strata = ~REGION1 + REGION2,  # Strata variables
  weights = ~I_WGC,  # Weight variable
  # data = a_ind
  data = subset(a_ind, I_INC1 != 0)  # Filtered dataset
)

# svymean(~I_INC1, design = survey_design, vartype = c("se"))
# svytotal(~I_INC1, design = survey_design, vartype = c("se"))

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
