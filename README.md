# GECal

Generalized Entropy Calibration produces calibration weights based on generalized entropy as the objective function for optimization. In GECcal, design weights play a role in the constraints to ensure design consistency, rather than being part of the objective function.

**Paper**: Kwon, Y., Kim, J., & Qiu, Y. (2024). [Debiased calibration estimation using generalized entropy in survey sampling.](https://arxiv.org/abs/2404.01076) Submitted.  
<!--Kwon, Y., & Kim, J. (2023). [Ensemble Fractional Imputation for Incomplete Categorical Data with a Graphical Model.](https://dmlr.ai/assets/accepted-papers/135/CameraReady/DMLR_paper.pdf) *In Workshop on Data-centric Machine Learning Research, International Conference on Machine Learning (ICML).* -->

## Installation
GECal can be installed using the package *devtools*. Please install *devtools* if it is unavailable.
``` r
# install.packages("devtools") # Install "devtools" if it is unavailable.
devtools::install_github("yonghyun-K/GECal", dependencies = T)

library(GECal)
```

## Example commands
``` r
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
```

## Externel Links
<!--
-->
- [CRAN Task View: CRAN Task View: Official Statistics & Survey Statistics]([https://cran.r-project.org/web/views/MissingData.html](https://cran.r-project.org/web/views/OfficialStatistics.html))

- [sampling](https://cran.r-project.org/web/packages/sampling/)

- [laeken](https://cran.r-project.org/web/packages/laeken/)
