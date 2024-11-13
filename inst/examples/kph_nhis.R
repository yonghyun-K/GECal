# nhis = read.csv("C:/Users/ghkfk/Box/Data/NHIS2023/NHIS_2023.CSV")
nhis = read.csv("C:/Users/ghkfk/Box/Data/NHIS2021/NHIS_2021.CSV")

names(nhis)[names(nhis) == "Height"] <- "HT"
names(nhis)[names(nhis) == "Weight"] <- "WT"
names(nhis)[names(nhis) == "ProvCode"] <- "REGION1"
names(nhis)[names(nhis) == "Gender"] <- "SEX"
nhis$REGION1 <- as.factor(nhis$REGION1)
nhis$SEX <- as.factor(nhis$SEX)
nhis$AgeGroup <- as.factor(nhis$AgeGroup)

kph = read.csv("C:/Users/ghkfk/Box/Data/KHP2.2_sas/2021_ind.csv")
kph = kph[!is.na(kph$SEX),]
kph = kph[!is.na(kph$HT),]

kph = kph[2021 - kph$BIRTH_Y >= 20,]

kph$REGION1 <- as.factor(kph$REGION1)
kph$SEX <- as.factor(kph$SEX)
# table(kph$SEX, useNA= "ifany")
# table(kph$REGION1, useNA= "ifany")

AgeGroup = cut(2021 - kph$BIRTH_Y, breaks = c((4:17) * 5, Inf), 
    labels = 5:18, right = FALSE)

kph$AgeGroup <- as.factor(AgeGroup)

table(kph$AgeGroup)

ncol(model.matrix(~ SEX + REGION1 + AgeGroup, data = kph))
head(2021 - kph$BIRTH_Y)
hist(2021 - kph$BIRTH_Y)

head(kph$CD1_HTN)

# sum((kph$CD1_HTN == 2), na.rm = T) / 
#   sum((!is.na(kph$CD1_HTN)))
# table(kph$CD1_HTN)

sum((kph$CD1_HTN == 1) * kph$I_WGC, na.rm = T) / 
  sum((!is.na(kph$CD1_HTN)) * kph$I_WGC)

form1 <- (~ SEX + REGION1 + AgeGroup + HT + WT)
form1 <- (~ HT + WT)

mean(kph$WT)
mean(nhis$WT)

sum(kph$HT * kph$I_WGC, na.rm = T)

total = colSums(model.matrix(form1, nhis)) * 50
total = colMeans(model.matrix(form1, nhis))

cal1 = GECal::GEcalib(form1, dweight = rep(1, nrow(kph)), 
               data = kph, const = total, method = "GEC0",
               entropy = "SL")

cal1 = GECal::GEcalib(form1, dweight = I_WGC, 
                      data = kph, const = total, method = "GEC0",
                      entropy = "EL")



plot(cal1$w,
kph$I_WGC)

cal1 = GECal::GEcalib(~ SEX, dweight = rep(1, nrow(kph)), 
                      data = kph, const = c(nrow(nhis), table(nhis$Gender)[-1]), 
                      method = "GEC0", entropy = "EL")

cal1 = GECal::GEcalib(~ SEX, dweight = rep(1, nrow(kph)), 
                      data = kph, const = c(nrow(nhis), table(nhis$Gender)[-1]), 
                      method = "DS", entropy = "SL")




View(nhis)
names(nhis)

head(nhis)

table(nhis$ProvCode)[1] * 50

length(table(nhis$ProvCode))

table(nhis$AgeGroup)





