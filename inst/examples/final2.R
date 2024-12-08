setwd("C:/Users/User/Box/Data")

library(survey)
library(dplyr)

c_ind = read.csv("KHP2.2_sas/c_ind.txt")
c_hh = read.csv("KHP2.2_sas/c_hh.txt")
c_ind$REGION1 <- NULL
c_ind <- left_join(c_ind, c_hh %>% select(HHID, REGION2, REGION1), by = join_by(HHID))
c_ind$REGION1 <- as.factor(c_ind$REGION1)
c_ind$REGION2 <- as.factor(c_ind$REGION2)
c_ind$AGE = 2021 - c_ind$BIRTH_Y

c_ind_new <- c_ind %>% subset(!is.na(WT)) %>% subset(AGE > 18) 
sort(sapply(c_ind_new, function(x) sum(is.na(x))))

c_ind_new$REGION_SEX <- interaction(
  as.factor(c_ind_new$REGION1),
  as.factor(c_ind_new$SEX),
  drop = TRUE
)

# HN22 = readxl::read_xlsx("KNHANES 2022/HN22_ALL.xlsx")
HN22_ALL = read.csv("KNHANES 2022/hn22_all.txt")

HN22_ALL$DI1_ag_group <- 
  cut(HN22_ALL$DI1_ag, breaks = c(-Inf, 19, 29, 49, 64, 887, 888, 999),
      labels = c(1:5, 888, 999),
      right = T)

HN22_ALL$region_sex <- interaction(
  as.factor(HN22_ALL$region),
  as.factor(HN22_ALL$sex),
  drop = TRUE
)

# Create the survey design object
survey_design <- svydesign(
  id = ~psu,              # Specify the cluster variable
  strata = ~kstrata,      # Specify the strata variable
  weights = ~wt_itvex,    # Specify the weight variable
  data = HN22_ALL,        # Specify the dataset
  nest = TRUE             # Specify nesting (if applicable)
)

sex_sum <- unclass(svytotal(~as.factor(sex),
    design = subset(survey_design, age > 18), na.rm=T))
attributes(sex_sum) <- NULL

region_sum <- unclass(svytotal(~as.factor(region),
    design = subset(survey_design, age > 18), na.rm=T))
attributes(region_sum) <- NULL

region_sex_sum <- svytotal(~ region_sex,
    design = subset(survey_design, age > 18), na.rm = TRUE)
attributes(region_sex_sum) <- NULL

age_sum <- unclass(svytotal(~age,
    design = subset(survey_design, age > 18), na.rm=T))
attributes(age_sum) <- NULL

HE_ht_sum <- unclass(svytotal(~HE_ht,
    design = subset(survey_design, age > 18), na.rm=T))
attributes(HE_ht_sum) <- NULL

HE_wt_sum <- unclass(svytotal(~HE_wt,
   design = subset(survey_design, age > 18), na.rm=T))
attributes(HE_wt_sum) <- NULL

c_ind_new$SEX <- factor(c_ind_new$SEX)
c_ind_new$REGION1 <- factor(c_ind_new$REGION1)
c_ind_new$REGION_SEX <- factor(c_ind_new$REGION_SEX)
c_ind_new$CD1_HTN <- factor(c_ind_new$CD1_HTN)
c_ind_new$CD2_HTN <- factor(c_ind_new$CD2_HTN)

# calibration = GECal::GEcalib(~ 0 + REGION_SEX + WT + HT, method = "GEC0", entropy = "ET",
#                              const = c(region_sex_sum, HE_wt_sum, HE_ht_sum),
#                              dweight = rep(1), data = c_ind_new)
# 
# calibration = GECal::GEcalib(~ 0 + REGION_SEX + AGE, method = "GEC0", entropy = "ET",
#                              const = c(region_sex_sum, age_sum),
#                              dweight = rep(1), data = c_ind_new)

calibration_ET = GECal::GEcalib(~ 0 + REGION_SEX + AGE + WT, method = "GEC0", entropy = "ET",
                             const = c(region_sex_sum, age_sum, HE_wt_sum),
                             dweight = rep(1), data = c_ind_new)

calibration_HD = GECal::GEcalib(~ 0 + REGION_SEX + AGE + WT, method = "GEC0", entropy = "HD",
                             const = c(region_sex_sum, age_sum, HE_wt_sum),
                             dweight = rep(1), data = c_ind_new)

calibration_CE = GECal::GEcalib(~ 0 + REGION_SEX + AGE + WT, method = "GEC0", entropy = "CE",
                             const = c(region_sex_sum, age_sum, HE_wt_sum),
                             dweight = rep(10), data = c_ind_new)

calibration_SL = GECal::GEcalib(~ 0 + REGION_SEX + AGE + WT, method = "GEC0", entropy = "SL",
                                const = c(region_sex_sum, age_sum, HE_wt_sum),
                                dweight = rep(10), data = c_ind_new)

# summary(calibration$w)
# benchmark ############

survey_design_c_ind <- svydesign(
  id = ~1,  # Assuming no clustering
  strata = ~REGION1 + REGION2,  # Strata variables
  weights = ~I_WGC,  # Weight variable
  # data = c_ind[!is.na(c_ind$REGION1) & !is.na(c_ind$HS1), ]
  data = c_ind[!is.na(c_ind$REGION1), ]
  # data = c_ind
)

calibrations <- list(
  calibration_ET, 
  calibration_HD, 
  calibration_CE, 
  calibration_SL
)

calibration_names <- c("survey_design_ET", "survey_design_HD", "survey_design_CE", "survey_design_SL")
res_names <- c("res_ET", "res_HD", "res_CE", "res_SL")

for (i in seq_along(calibrations)) {
  survey_design0 = svydesign(
    id = ~1,  # Assuming no clustering
    strata = ~REGION1 + REGION2,  # Strata variables
    weights = ~calibrations[[i]]$w,  # Weight variable
    data = c_ind_new
  )
  
  res_tmp = unclass(svytotal(~factor(CD2_HTN), design = survey_design0, na.rm=T))
  res_tmp = data.frame(cbind(total = res_tmp, SE = sqrt(diag(attr(res_tmp, "var"))))[-1,])
  
  assign(
    res_names[i],
    res_tmp
  )
  
}

svytotal(~factor(CD1_HTN), design = survey_design_c_ind, na.rm=T)

svytotal(~factor(CD1_HTN), design = survey_design0, na.rm=T)

svytotal(~as.factor(DI1_dg), design = survey_design, na.rm=T)

unlist(svytotal(~factor(CD2_HTN), design = survey_design_c_ind, na.rm=T))

res1 = unclass(svytotal(~factor(CD2_HTN), design = survey_design_c_ind, na.rm=T))
res1 = data.frame(cbind(total = res1, SE = sqrt(diag(attr(res1, "var"))))[-1,])

res2 = unclass(svytotal(~as.factor(DI1_ag_group), design = survey_design, na.rm=T))
res2 = data.frame(cbind(total = res2, SE = sqrt(diag(attr(res2, "var"))))[c(-6, -7),])

res3 = unclass(svytotal(~factor(CD2_HTN), design = survey_design0, na.rm=T))
res3 = data.frame(cbind(total = res3, SE = sqrt(diag(attr(res3, "var"))))[-1,])


row.names(res1) <- row.names(res2) <- row.names(res3) <- NULL

# xtable::xtable(cbind(res1, res2, res3) / 1e3, digits = 0)

xtable::xtable(cbind(res1, res2, res_ET, res_HD, res_SL, res_CE) / 1e3, digits = 0)

# # Plotting
# par(mfrow = c(5, 1), mar = c(4, 4, 2, 2)) # Set 5 rows, 1 column layout
# colors <- c("red", "blue", "green")
# labels <- c("Naive", "HT", "GEC")
# 
# for (i in 1:5) {
#   # Extract data for the row
#   x_vals <- c(res1$total[i], res2$total[i], res3$total[i])
#   y_vals <- c(1, 2, 3) # Position for each res
#   lower <- x_vals - c(res1$SE[i], res2$SE[i], res3$SE[i]) * 2
#   upper <- x_vals + c(res1$SE[i], res2$SE[i], res3$SE[i]) * 2
#   
#   # Plot
#   plot(NULL, xlim = range(c(lower, upper)), ylim = c(0.5, 3.5),
#        xlab = "Total", ylab = "", yaxt = "n", main = paste("Row", i, "CIs"))
#   axis(2, at = 1:3, labels = labels)
#   for (j in 1:3) {
#     points(x_vals[j], y_vals[j], pch = 16, col = colors[j])
#     lines(c(lower[j], upper[j]), c(y_vals[j], y_vals[j]), col = colors[j], lwd = 2)
#   }
# }


library(ggplot2)

resnames<- c("Naive", "HT", "KL", "HD", "SL", "SKL")

Rows <- c("19 years old and below", "20–29 years old", "30–49 years old",
          "50–64 years old", "65 years old and above")


# Create a combined data frame
data <- bind_rows(
  res1 %>% mutate(Label = resnames[1], Row = Rows),
  res2 %>% mutate(Label = resnames[2], Row = Rows),
  res_ET %>% mutate(Label = resnames[3], Row = Rows),
  res_HD %>% mutate(Label = resnames[4], Row = Rows),
  res_SL %>% mutate(Label = resnames[5], Row = Rows),
  res_CE %>% mutate(Label = resnames[6], Row = Rows)
)

# Calculate lower and upper bounds for CIs
data <- data %>%
  mutate(Lower = total - 2 * SE,
         Upper = total + 2 * SE)
# Reorder the factor levels for Label
data <- data %>%
  mutate(Label = factor(Label, levels = rev(resnames)))  # Reverse order for y-axis

# Plot
ggplot(data, aes(x = total, y = Label, color = Label)) +
  geom_point(size = 2.5) +  # Smaller points for a clean look
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, size = 0.7) + # Thinner error bars
  facet_wrap(~ Row, ncol = 1, scales = "free_x", strip.position = "top") + # Allow free x-scales
  # scale_color_manual(values = c("black", "darkgray", "lightgray")) + # Grayscale color scheme
  scale_color_manual(values = c("black", "gray20", "gray40", "gray60", "gray70", "gray80")) + # Grayscale color scheme
  labs(x = "Estimated Total", y = "", title = "") + # Simplified axis titles
  theme_minimal(base_size = 12) + # Smaller base size for journal style
  theme(
    strip.text = element_text(size = 10, face = "bold", hjust = 0.5), # Centered and bold strip text
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold"),
    panel.grid.major = element_line(color = "lightgray", linetype = "dotted"), # Subtle grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(color = "black", size = 0.5), # Add a border for clarity
    legend.position = "none" # No legend for a cleaner look
  )


