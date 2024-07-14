# install.packages("sampling") # Install "sampling" if it is unavailable.
library(sampling)
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
calib(Xs, total, d = d, method="raking") * d
GEcalib(Xs, total, d, entropy = "ET")
GEcalib(Xs, total, entropy = "ET")

calib(Xs, total, d = d, method="linear") * d
GEcalib(Xs, total, d, entropy = "SL")
GEcalib(Xs, total, entropy = "SL")

# GEcalib(Xs, total, d, entropy = "EL")

############
## Example 2
############
# Example of g-weights (linear, raking, truncated, logit),
# with the data of Belgian municipalities as population.
# Firstly, a sample is selected by means of Poisson sampling.
# Secondly, the g-weights are calculated.
data(belgianmunicipalities)
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
pik=inclusionprobabilities(averageincome,200)
N=length(pik)               # population size
s=UPpoisson(pik)            # sample
Xs=X[s==1,]                 # sample matrix of calibration variables
piks=pik[s==1]; d = 1 / piks# sample inclusion probabilities
n=length(piks)              # expected sample size
# vector of population totals of the calibration variables
total=c(t(rep(1,times=N))%*%X)  
# computation of the g-weights
# by means of different calibration methods

all.equal(calib(Xs, total, d = d, method="linear") * d,
          GEcalib(Xs, total, d = d, entropy = "SL"))
GEcalib(cbind(Xs, d), c(total, sum(1 / pik)), entropy = "SL")

all.equal(calib(Xs, total, d = d, method="raking") * d,
          GEcalib(Xs, total, d = d, entropy = "ET"))
GEcalib(cbind(Xs, log(d)), c(total, sum(log(1 / pik))), entropy = "ET")

all.equal(GEcalib(Xs, total, d = d, entropy = "EL"),
          GEcalib(cbind(Xs[,-ncol(Xs)], -Xs[,ncol(Xs)]), c(total[-ncol(Xs)], -total[ncol(Xs)]), d = d, 
                  entropy = "EL"))

GEcalib(cbind(Xs[,-ncol(Xs)], -piks), c(total[-ncol(Xs)], sum(-pik)), 
        entropy = "EL")

# GEcalib(cbind(Xs), c(total), entropy = "EL")


GEcalib(Xs, total, d = d, entropy = "HD")
GEcalib(cbind(Xs, -sqrt(1 / d)), c(total, sum(-sqrt(pik))), entropy = "HD")
