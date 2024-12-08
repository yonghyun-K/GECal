setwd("C:/Users/User/Box/Data")

library(survey)
library(dplyr)

c_ind = read.csv("KHP2.2_sas/c_ind.txt")
c_hh = read.csv("KHP2.2_sas/c_hh.txt")
c_ind$REGION1 <- NULL
c_ind <- left_join(c_ind, c_hh %>% select(HHID, REGION2, REGION1), by = join_by(HHID))
c_ind$REGION1 <- as.factor(c_ind$REGION1)
c_ind$REGION2 <- as.factor(c_ind$REGION2)
c_ind$AGE = 2022 - c_ind$BIRTH_Y


sort(sapply(c_ind, function(x) sum(is.na(x))))

# c_ind <- left_join(c_ind, c_hh %>% select(HHID, REGION2), by = join_by(HHID))

survey_design0 <- svydesign(
  id = ~1,  # Assuming no clustering
  strata = ~REGION1 + REGION2,  # Strata variables
  weights = ~I_WGC,  # Weight variable
  # data = c_ind[!is.na(c_ind$REGION1) & !is.na(c_ind$HS1), ]
  data = c_ind[!is.na(c_ind$REGION1), ]
  # data = c_ind
)

# 성별
svytotal(~as.factor(SEX), design = survey_design0, na.rm=T)
sex_sum <- unclass(svytotal(~as.factor(SEX), design = survey_design0, na.rm=T))
# sex_sum <- unclass(svytotal(~as.factor(SEX), 
#     design = subset(survey_design0, !is.na(CD1_HTN)), na.rm=T))
attributes(sex_sum) <- NULL

# # 지역
# original_values <- levels(unique(c_ind[!is.na(c_ind$REGION1), "REGION1"]))
# new_values <- 1:17
# 
# mapping <- setNames(new_values, original_values)
# c_ind$REGION1 <- mapping[as.character(c_ind$REGION1)] # Map using the mapping vector

svytotal(~as.factor(REGION1), design = survey_design0, na.rm=T)
region_sum <- unclass(svytotal(~as.factor(REGION1), design = survey_design0, na.rm=T))
attributes(region_sum) <- NULL

# 성별, 지역
survey_design0$variables$region_sex <- interaction(
  as.factor(survey_design0$variables$REGION1),
  as.factor(survey_design0$variables$SEX),
  drop = TRUE
)
region_sex_sum <- svytotal(~ region_sex, design = survey_design0, na.rm = TRUE)
attributes(region_sex_sum) <- NULL

# 나이
svytotal(~as.factor(AGE), design = survey_design0, na.rm=T)
age_sum <- unclass(svytotal(~AGE, design = survey_design0, na.rm = TRUE))
# age_sum <- unclass(svytotal(~AGE, 
#     design = subset(survey_design0, !is.na(CD1_HTN)), na.rm=T))
attributes(age_sum) <- NULL

# 신장(cm)
svymean(~HT, design = survey_design0, na.rm = TRUE)
HE_ht_sum <- unclass(svytotal(~HT, design = survey_design0, na.rm = TRUE))
# HE_ht_sum <- unclass(svytotal(~HT, 
#     design = subset(survey_design0, !is.na(CD1_HTN)), na.rm=T))
attributes(HE_ht_sum) <- NULL

svymean(~HT, design = subset(survey_design0, AGE >= 20), na.rm = TRUE) # male height


svymean(~HT, design = subset(survey_design0, SEX == 1 &
                               AGE >= 20 & AGE < 70), na.rm = TRUE) # male height

svymean(~HT, design = subset(survey_design0, SEX == 2 &
                               AGE >= 20 & AGE < 70), na.rm = TRUE) # female height


# 체중
svymean(~WT, design = survey_design0, na.rm = TRUE)

HE_wt_sum <- unclass(svytotal(~WT, design = survey_design0, na.rm = TRUE))
# HE_wt_sum <- unclass(svytotal(~WT, 
#    design = subset(survey_design0, !is.na(CD1_HTN)), na.rm=T))
attributes(HE_wt_sum) <- NULL

# 만성질환 진단: 고혈압
svytotal(~as.factor(CD1_HTN), design = survey_design0, na.rm=T)
HE_HP_sum <- unclass(svytotal(~as.factor(CD1_HTN), design = survey_design0, na.rm=T))
attributes(HE_HP_sum) <- NULL

# 고혈압 환자 나이
table(c_ind[!is.na(c_ind$REGION1) & c_ind$CD1_HTN == 1,"AGE"]) 

table(c_ind$CD1_HTN, useNA = "ifany")

# 만성질환 진단 시기: 고혈압
svytotal(~as.factor(CD2_HTN), design = survey_design0, na.rm=T)
# svytotal(~as.factor(DI1_ag_group), design = survey_design, na.rm=T)
DI1_ag_group_sum <- unclass(svytotal(~as.factor(CD2_HTN), design = survey_design0, na.rm=T))
attributes(DI1_ag_group_sum) <- NULL
DI1_ag_group_sum <- DI1_ag_group_sum[-1]
# sum(DI1_ag_group_sum)
# sum(unclass(svytotal(~as.factor(DI1_ag_group), design = survey_design, na.rm=T))[c(-6, -7)])

# 만성질환 진단: 당뇨
svytotal(~as.factor(CD1_DM), design = survey_design0, na.rm=T)

# 만성질환 진단: 천식
svytotal(~as.factor(CD1_AST), design = survey_design0, na.rm=T)

#  만성질환 유무_ 무릎관절증
svytotal(~as.factor(CD1_OAK), design = survey_design0, na.rm=T)
table(c_ind$CD1_OAK, useNA = "ifany")

# 1 년간 자살 생각 여부
svytotal(~as.factor(HS4_YN), design = survey_design0, na.rm=T)

# 걷기-최근 일주일 하루 걸은 시간(분)
svymean(~P2_2, design = survey_design0, na.rm = TRUE)

# 연간 외래서비스 이용 수납금액(원)
svymean(~OUOOP_1, design = survey_design0, na.rm = TRUE)

which(names(c_ind) == "CD1_HTN")
which(names(c_ind) == "CD1_1")

# Create the regression formula dynamically
formula <- as.formula(paste("OUOOP_1 ~", paste(names(c_ind[,(19:25) * 2]), collapse = " + ")))

summary(lm(formula, data = c_ind, na.action = na.omit))


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

sort(sapply(c_ind, function(x) sum(is.na(x))))
table(HN22_ALL$DI1_dg, useNA = "ifany")
table(HN22_ALL$DI1_ag, useNA = "ifany")

sum(c_ind[is.na(c_ind$WT),"AGE"] >= 19, na.rm = T)

table(c_ind[is.na(c_ind$WT),"SEX"], useNA = "ifany")

table(c_ind[!is.na(c_ind$SEX),"AGE"])

table(c_ind[!is.na(c_ind$SEX),"WT"], useNA = "ifany")

table(c_ind[c_ind$CD1_HTN == 1,"AGE"])

table(c_ind$REGION1, useNA = "ifany")
table(c_ind$D1, useNA = "ifany")

table(c_ind$CD1_HTN, useNA = "ifany")
sum(table(c_ind$CD2_HTN, useNA = "ifany")[2:6])

table(c_ind$CD1_CLD, useNA = "ifany")

table(c_ind$CD, useNA = "ifany")
table(c_ind$REGION1, useNA = "ifany")

# Create the survey design object
survey_design <- svydesign(
  id = ~psu,              # Specify the cluster variable
  strata = ~kstrata,      # Specify the strata variable
  weights = ~wt_itvex,    # Specify the weight variable
  data = HN22_ALL,        # Specify the dataset
  nest = TRUE             # Specify nesting (if applicable)
)

# 지역
svytotal(~ as.factor(region), design = survey_design, na.rm = TRUE)

# 성별, 지역
svytotal(~ region_sex, design = survey_design, na.rm = TRUE)

# Compute survey means for the variable HE_BMI
# 체질량지수
svymean(~HE_BMI, design = survey_design, na.rm = TRUE)

# 성별
svytotal(~as.factor(sex), design = survey_design, na.rm=T)

# 나이
svytotal(~as.factor(age), design = survey_design, na.rm=T)

# 신장(cm)
svymean(~HE_ht, design = survey_design, na.rm = TRUE)

# 체중
svymean(~HE_wt, design = survey_design, na.rm = TRUE)

# Compute the mean height for males
svymean(~HE_ht, design = subset(survey_design, sex == 1 &
                                age >= 20 & age < 70), na.rm = TRUE) # male height

svymean(~HE_ht, design = subset(survey_design, sex == 2 &
                                  age >= 20 & age < 70), na.rm = TRUE) # female height

# Compute the frequency table for HE_HP
#고혈압 유병여부
svytotal(~as.factor(HE_HP), design = survey_design, na.rm=T)

#고혈압 의사진단 여부
svytotal(~as.factor(DI1_dg), design = survey_design, na.rm=T)

HN22_ALL_tmp <- HN22_ALL[,c("DI1_ag","DI1_dg", "age","sex", "HE_wt", "HE_ht", "DE1_dg",
                            "wt_itvex", "psu", "kstrata")]

HN22_ALL_tmp <- HN22_ALL[,c("DI1_ag","DI1_dg", "age","sex", "DE1_dg", "region", "region_sex",
                            "wt_itvex", "psu", "kstrata")]

HN22_ALL_tmp <- HN22_ALL[,c("DI1_ag","DI1_ag_group", "DI1_dg", "DI1_2",
                            "age", "sex", "DE1_dg", "region", "region_sex",
                            "wt_itvex", "psu", "kstrata")]

# sapply(HN22_ALL, function(x) sum(is.na(x) & !is.na(HN22_ALL$DI1_dg)) == 0)
table(HN22_ALL$DI1_dg, useNA = "ifany"); sum(complete.cases(HN22_ALL_tmp))
# sum(HN22_ALL$DI1_ag < 888 & HN22_ALL$DI1_dg == 1, na.rm = T)
if(!is.null(HN22_ALL_tmp$DI1_dg)){
  HN22_ALL_tmp <- HN22_ALL_tmp %>% subset(DI1_dg != 8 & DI1_dg != 9)
  HN22_ALL_tmp$DI1_dg <- as.factor(HN22_ALL_tmp$DI1_dg) 
}
# if(!is.null(HN22_ALL_tmp$DI1_ag_group)){
#   HN22_ALL_tmp <- HN22_ALL_tmp %>% subset(DI1_ag_group != 888 & DI1_ag_group != 999)
#   HN22_ALL_tmp$DI1_ag_group <- factor(HN22_ALL_tmp$DI1_ag_group) 
# }
if(!is.null(HN22_ALL_tmp$sex)) HN22_ALL_tmp$sex <- as.factor(HN22_ALL_tmp$sex)
if(!is.null(HN22_ALL_tmp$region)) HN22_ALL_tmp$region <- as.factor(HN22_ALL_tmp$region)
if(!is.null(HN22_ALL_tmp$region_sex)) HN22_ALL_tmp$region_sex <- as.factor(HN22_ALL_tmp$region_sex)
HN22_ALL_tmp <- HN22_ALL_tmp[complete.cases(HN22_ALL_tmp),]
nrow(HN22_ALL_tmp)
# model.matrix(~ 0 + DI1_dg + sex+ age + HE_wt + HE_ht, data = HN22_ALL_tmp)

# 고혈압 진단시기
# calibration = GECal::GEcalib(~ 0 + DI1_dg + sex + age, method = "DS", entropy = "SL",
#                const = c(c(HE_HP_sum[2], HE_HP_sum[1]), sex_sum[2], age_sum),
#                dweight = wt_itvex, data = HN22_ALL_tmp)

calibration = GECal::GEcalib(~ 0 + DI1_dg + sex, method = "DS", entropy = "SL",
                             const = c(c(HE_HP_sum[2], HE_HP_sum[1]), sex_sum[2]),
                             dweight = wt_itvex, data = HN22_ALL_tmp)

# colSums(model.matrix(~ 0 + sex + DI1_ag_group, data = HN22_ALL_tmp))

# calibration = GECal::GEcalib(~ 0 + DI1_dg, method = "DS", entropy = "SL",
#                              const = c(c(HE_HP_sum[2], HE_HP_sum[1])),
#                              dweight = wt_itvex, data = HN22_ALL_tmp)

# calibration = GECal::GEcalib(~ 0 + sex, method = "DS", entropy = "SL",
#                              const = sex_sum,
#                              dweight = wt_itvex, data = HN22_ALL_tmp)

summary(calibration$w)

cor(calibration$w, HN22_ALL_tmp$wt_itvex)

survey_design_subset <- svydesign(
  id = ~psu,              # Specify the cluster variable
  strata = ~kstrata,      # Specify the strata variable
  weights = ~wt_itvex,    # Specify the weight variable
  data = HN22_ALL_tmp,        # Specify the dataset
  nest = TRUE             # Specify nesting (if applicable)
)

survey_design_subset2 <- svydesign(
  id = ~psu,              # Specify the cluster variable
  strata = ~kstrata,      # Specify the strata variable
  weights = ~calibration$w,    # Specify the weight variable
  data = HN22_ALL_tmp,        # Specify the dataset
  nest = TRUE             # Specify nesting (if applicable)
)

nrow(HN22_ALL %>% subset(DI1_ag != 888 & DI1_ag != 999))

nrow(subset(survey_design, DI1_ag != 888 & DI1_ag != 999))

svymean(~DI1_ag, design = subset(survey_design_subset, DI1_ag != 888 & DI1_ag != 999), na.rm = TRUE)

svymean(~DI1_ag, design = subset(survey_design_subset2, DI1_ag != 888 & DI1_ag != 999), na.rm = TRUE)

svymean(~DI1_ag, design = subset(survey_design, DI1_ag != 888 & DI1_ag != 999), na.rm = TRUE)

# 고혈압 진단시기 그룹
svytotal(~as.factor(DI1_ag_group), design = survey_design, na.rm=T)

# 혈압조절제 복용 정도
svytotal(~as.factor(DI1_2), design = survey_design, na.rm=T)

#당뇨병 유병여부
svytotal(~as.factor(DE1_dg), design = survey_design, na.rm=T)

# 이상지질혈증
svytotal(~as.factor(DI2_dg), design = survey_design, na.rm=T)

# 천식
svytotal(~as.factor(DJ4_dg), design = survey_design, na.rm=T)

# 아토피피부염
svytotal(~as.factor(DL1_dg), design = survey_design, na.rm=T)

# 알레르기비염
svytotal(~as.factor(DJ8_dg), design = survey_design, na.rm=T)

# 부비동염
svytotal(~as.factor(DJ6_dg), design = survey_design, na.rm=T)

# 중이염
svytotal(~as.factor(DH4_dg), design = survey_design, na.rm=T)

# 주중 또는 일하는 날 하루 평균 수면 시간
svymean(~BP16_1, design = subset(survey_design, BP16_1 != 88 & BP16_1 != 99), na.rm=T)

# 주말 또는 일하지 않는 날 일하지 않는 전날 하루 평균 수면 시간
svymean(~BP16_2, design = subset(survey_design, BP16_2 != 88 & BP16_2 != 99), na.rm=T)

# 1 년간 2주이상 연속 우울감 여부
svytotal(~as.factor(BP5), design = subset(survey_design, BP5 != 8), na.rm=T)

# 1 년간 자살 생각 여부
svytotal(~as.factor(BP6_10), design = survey_design, na.rm=T)

# 1 년간 자살 계획 여부
svytotal(~as.factor(BP6_2), design = survey_design, na.rm=T)

# 걷기 지속 시간 분
svymean(~BE3_33, design = subset(survey_design, BE3_33 != 88 & BE3_33 != 99), na.rm=T)

# 고혈압 유병 여부
svytotal(~as.factor(DI1_pr), design = subset(survey_design, DI1_pr != 8 & DI1_pr != 9), na.rm=T)

# 고혈압 의사진단 여부
svytotal(~as.factor(DI1_dg), design = subset(survey_design, DI1_dg != 8 & DI1_dg != 9), na.rm=T)
table(HN22_ALL$DI1_dg, useNA = "ifany")

summary(HN22_ALL$HE_TG) # 중성지방
summary(HN22_ALL$HE_chol) # 총콜레스테롤

summary(HN22$wt_itvex) # 건강설문 검진조사 가중치


summary(HN22$체질량지수)



summary(HN22$)

summary(HN22$`체성분,영양 가중치`)
