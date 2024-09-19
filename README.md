# GECal

Generalized Entropy Calibration produces calibration weights based on generalized entropy as the objective function for optimization. In GECcal, design weights play a role in the constraints to ensure design consistency, rather than being part of the objective function.

**Paper**: Kwon, Y., Kim, J., & Qiu, Y. (2024). [Debiased calibration estimation using generalized entropy in survey sampling.](https://arxiv.org/abs/2404.01076) Submitted.  

## Installation
GECal is now available in CRAN. Use install.packages for installation.
``` r
install.packages("GECal")
```

GECal can instead be installed using the package *devtools*. Please install *devtools* if it is unavailable.
``` r
install.packages("devtools") # Install "devtools" if it is unavailable.
devtools::install_github("yonghyun-K/GECal", dependencies = T)
```

## Example commands
``` r
library(GECal)

set.seed(11)
N = 10000
x = data.frame(x1 = rnorm(N, 2, 1), x2= runif(N, 0, 4))
pi = pt((-x[,1] / 2 - x[,2] / 2), 3);
pi = ifelse(pi >.7, .7, pi)

delta = rbinom(N, 1, pi)
Index_S = (delta == 1)
pi_S = pi[Index_S]; d_S = 1 / pi_S
x_S = x[Index_S,,drop = FALSE]
# pimat = diag(d_S^2 - d_S) / N^2 # 1 / pi_i * (1 - 1 / pi_i)

e = rnorm(N, 0, 1)
y = x[,1] + x[,2] + e;
y_S = y[Index_S] # plot(x_S, y_S)

# Hajek estimator of population mean
calibration0 <- GECal::GEcalib(~ 1, dweight = d_S, data = x_S,
                               const = 1, 
                               entropy = "SL", method = "DS")
GECal::estimate(y_S ~ 1, calibration = calibration0)$estimate 
# sum(y_S * d_S) / sum(d_S)

# Hajek estimator of population total
calibration0 <- GECal::GEcalib(~ 1, dweight = d_S, data = x_S,
                               const = N,
                               entropy = "SL", method = "DS")
GECal::estimate(y_S ~ 1, calibration = calibration0)$estimate 
# sum(y_S * d_S) * N / sum(d_S)

# HT estimator of population total
calibration <- GECal::GEcalib(~ 0, dweight = d_S, data = x_S,
                              const = numeric(0),
                              entropy = "SL", method = "DS")
GECal::estimate(y_S ~ 1, calibration = calibration)$estimate

# DS estimator using ET(exponential tilting) divergence
calibration1 <- GECal::GEcalib(~ ., dweight = d_S, data = x_S,
                               const = colSums(cbind(1, x)),
                               entropy = "ET", method = "DS")
GECal::estimate(y_S ~ 1, calibration = calibration1)$estimate

# GEC0 estimator using ET entropy
calibration2 <- GECal::GEcalib(~ ., dweight = d_S, data = x_S,
                               const = colSums(cbind(1, x)),
                               entropy = "ET", method = "GEC0")
GECal::estimate(y_S ~ 1, calibration = calibration2)$estimate

# GEC estimator using ET entropy 
# when the population total for log(d_S) is known.
calibration3 <- GECal::GEcalib(~ . + g(d_S), dweight = d_S, data = x_S,
                               const = colSums(cbind(1, x, log(1 / pi))),
                               entropy = "ET", method = "GEC")
GECal::estimate(y_S ~ 1, calibration = calibration3)$estimate

# GEC estimator using ET entropy 
# when the population total for log(d_S) is unknown.
calibration4 <- GECal::GEcalib(~ . + g(d_S), dweight = d_S, data = x_S,
                               const = colSums(cbind(1, x, NA)),
                               entropy = "ET", method = "GEC")
GECal::estimate(y_S ~ 1, calibration = calibration4)$estimate

# GEC estimator using ET entropy using different K_alpha
# when the population total for log(d_S) is unknown.
calibration5 <- GECal::GEcalib(~ . + g(d_S), dweight = d_S, data = x_S,
                               const = colSums(cbind(1, x, NA)),
                               entropy = "ET", method = "GEC", K_alpha = "log")
GECal::estimate(y_S ~ 1, calibration = calibration5)$estimate

```

## Externel Links
<!--
-->
- [CRAN Task View: Official Statistics & Survey Statistics](https://CRAN.R-project.org/view=OfficialStatistics)

- [survey](https://CRAN.R-project.org/package=survey)

- [sampling](https://CRAN.R-project.org/package=sampling)

<!--
- [laeken](https://CRAN.R-project.org/package=laeken)
-->
