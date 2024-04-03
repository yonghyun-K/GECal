# GECal

Generalized Entropy Calibration produces calibration weights based on generalized entropy as the objective function for optimization. In GECcal, design weights play a role in the constraints to ensure design consistency, rather than being part of the objective function.

**Paper**: *   **Kwon, Y.**, Kim, J., & Qiu, Y. (2024). [Debiased calibration estimation using generalized entropy in survey sampling.](https://arxiv.org/abs/2404.01076) Submitted.  
<!--Kwon, Y., & Kim, J. (2023). [Ensemble Fractional Imputation for Incomplete Categorical Data with a Graphical Model.](https://dmlr.ai/assets/accepted-papers/135/CameraReady/DMLR_paper.pdf) *In Workshop on Data-centric Machine Learning Research, International Conference on Machine Learning (ICML).* -->

## Installation

``` r
library(devtools)
install_github("yonghyun-K/GECal", dependencies = T)

library(GECal)
```

## Example commands
<!--
``` r
library(EFI)

# Import data and generate missingness.
Y = as.data.frame.table(HairEyeColor, stringsAsFactors = TRUE)
Y = Y[rep(seq_len(nrow(Y)), Y$Freq), ]; Y$Freq <- NULL # OR # Y = tidyr::uncount(Y, Freq)
n = nrow(Y); p = ncol(Y); rownames(Y) <- 1:n
delta = matrix(rbinom(n * p, 1, 0.5), nr = n, nc = p); Y[delta == 0] = NA

# Ensemble Fractional Imputation.
cand.edges = apply(combn(p, 2), 2, list)
dp = doublep(Y, cand.edges, freq = F)
plot(dp)
EFI = efi(Y, dp, freq = F)

estimate(EFI, "(Hair == \"Black\") & (Eye == \"Brown\")")
estimate(EFI, "(Hair == \"Black\") & (Sex == \"Male\")")
```
-->
## Externel Links
<!--
- [CRAN Task View: Missing Data](https://cran.r-project.org/web/views/MissingData.html)

- [FHDI](https://github.com/cran/FHDI)

- [mice](https://github.com/amices/mice)
- https://stefvanbuuren.name/fimd/

- [missForest](https://github.com/stekhoven/missForest)

- [GAIN](https://github.com/jsyoon0823/GAIN)
-->
