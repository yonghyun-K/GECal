# install.packages("sampling") # Install "sampling" if it is unavailable.
library(GECal)

############
## Example 1
############
# matrix of sample calibration variables 
Xs=cbind(
  c(1,1,1,1,1,1,1,1,1,1),
  c(1,1,1,1,1,0,0,0,0,0),
  c(1,2,3,4,5,6,7,8,9,10)
)
# inclusion probabilities
piks=rep(0.2,times=10); d=1/piks
# vector of population totals
total=c(50,24,290)

# Calibration weights
sampling::calib(Xs, total, d = d, method="raking") * d
calibration <- GEcalib(~ 0 + Xs, dweight = d, const = total, method = "DS", entropy = "ET")
calibration$w

sampling::calib(Xs, total, d = d, method="linear") * d
calibration <- GEcalib(~ 0 + Xs, dweight = d, const = total, method = "DS", entropy = "SL")
calibration$w

############
## Example 2
############
# Example of g-weights (linear, raking, truncated, logit),
# with the data of Belgian municipalities as population.
# Firstly, a sample is selected by means of Poisson sampling.
# Secondly, the g-weights are calculated.
data(belgianmunicipalities, package = "sampling")
attach(belgianmunicipalities)
# matrix of calibration variables for the population
X=cbind(1, 
        Men03/mean(Men03),
        Women03/mean(Women03),
        Diffmen,
        Diffwom,
        TaxableIncome/mean(TaxableIncome),
        Totaltaxation/mean(Totaltaxation),
        medianincome/mean(medianincome),
        averageincome/mean(averageincome))
# selection of a sample with expectation size equal to 200
# by means of Poisson sampling
# the inclusion probabilities are proportional to the average income 
pik=sampling::inclusionprobabilities(averageincome,200)
N=length(pik)               # population size
s=sampling::UPpoisson(pik)            # sample
Xs=X[s==1,]                 # sample matrix of calibration variables
piks=pik[s==1]; d = 1 / piks# sample inclusion probabilities
n=length(piks)              # expected sample size
# vector of population totals of the calibration variables
total=c(t(rep(1,times=N))%*%X)  
# computation of the g-weights
# by means of different calibration methods

all.equal(sampling::calib(Xs, total, d = d, method="linear") * d,
          unname(GEcalib(~ 0 + Xs, dweight = d, const = c(total), 
                  method = "DS", entropy = "SL")$w))
GEcalib(~ 0 + Xs + g(d), dweight = d, const = c(total, sum(g(1 / pik, 1))), 
                       method = "GEC", entropy = "SL")$w

all.equal(sampling::calib(Xs, total, d = d, method="raking") * d,
          unname(GEcalib(~ 0 + Xs, dweight = d, const = c(total), 
                         method = "DS", entropy = "ET")$w))
GEcalib(~ 0 + Xs + g(d), dweight = d, const = c(total, sum(g(1 / pik, 0))), 
        method = "GEC", entropy = "ET")$w

# GEcalib(~ 0 + Xs, dweight = d, const = c(total), 
#         method = "DS", entropy = "EL")$w

# GEcalib(~ 0 + Xs, dweight = d, const = c(total), 
#         method = "DS", entropy = "HD")$w

GEcalib(~ 0 + Xs + g(d), dweight = d, const = c(total, sum(g(1 / pik, -1/2))), 
        method = "GEC", entropy = "HD")$w
