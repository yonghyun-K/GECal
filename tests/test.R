# library(sampling)

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
calib(Xs,d=d,total,method="raking") * d
GEcalib(Xs, d, total, entropy = "ET", DS = T)
GEcalib(Xs, d, total, entropy = "ET", DS = F)

calib(Xs,d=d,total,method="linear") * d
GEcalib(Xs, d, total, entropy = "SL", DS = T)
GEcalib(Xs, d, total, entropy = "SL", DS = F)

# GEcalib(Xs, d, total, entropy = "EL")

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

all.equal(calib(Xs,d=d,total,method="linear") * d,
          GEcalib(Xs, d, total, entropy = "SL", DS = T))
GEcalib(cbind(Xs, d), d, c(total, sum(1 / pik)), entropy = "SL", DS = F)

all.equal(calib(Xs,d=d,total,method="raking") * d,
          GEcalib(Xs, d, total, entropy = "ET", DS = T))
GEcalib(cbind(Xs, log(d)), d, c(total, sum(log(1 / pik))), entropy = "ET", DS = F)

all.equal(GEcalib(Xs, d, total, entropy = "EL", DS = T),
GEcalib(cbind(Xs[,-ncol(Xs)], -Xs[,ncol(Xs)]), d, c(total[-ncol(Xs)], -total[ncol(Xs)]), 
        entropy = "EL", DS = F))

GEcalib(Xs, d, total, entropy = "HD", DS = T)
GEcalib(cbind(Xs, -sqrt(1 / d)), d, c(total, sum(-sqrt(pik))), entropy = "HD", DS = F)


# GEcalib(cbind(Xs, log(1 - piks)), d, c(total, sum(log(1 - pik))), entropy = "CE", DS = F)

