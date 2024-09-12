# GECal:::GEcalib(~Sepal.Length + Sepal.Width + g(Petal.Width), dweight = Petal.Width,
#          data = iris, const = rep(1, 4),
#          method = "GEC",entropy = 1)
# 
# GECal:::GEcalib(~Sepal.Length + Sepal.Width + g(Petal.Width), dweight = Petal.Width,
#                  data = iris, const = c(1,1,1,2),
#                  method = "GEC",entropy = 1, G.scale = 2)
# 
# GECal:::GEcalib(~Sepal.Length + Sepal.Width + g(Petal.Width), dweight = Petal.Width,
#                  data = iris, const = c(1,1,1,2 + 2 * log(2) * nrow(iris)),
#                  method = "GEC",entropy = 1, weight.scale = 2)
# 
# GECal:::GEcalib(~Sepal.Length + Sepal.Width + g(2 * Petal.Width), dweight = 2 * Petal.Width,
#                  data = iris, const = c(1,1,1,2),
#                  method = "GEC",entropy = 1)
# 
# GECal:::GEcalib(~ + g(2 * Sepal.Length), dweight = Sepal.Length,
#                  data = iris, const = c(870, 2 * 1500),
#                  method = "GEC",entropy = 0)
# 
# 
# c(sum(iris$Sepal.Length), sum(g(iris$Sepal.Length, entropy = 0) * iris$Sepal.Length))
# 
# GECal:::GEcalib(~ + g(2 * Petal.Width), dweight = 2 * Petal.Width,
#                  data = iris, const = c(1,1,1,2),
#                  method = "GEC",entropy = 0)
# 
# 
# 
# GECal:::GEcalib(~Sepal.Length + Sepal.Width, dweight = Petal.Width,
#          data = iris, const = rep(1, 3),
#          method = "DS",entropy = 1)
# 
# 
# 
# GECal:::GEcalib(~iris$Sepal.Length + iris$Sepal.Width + g(iris$Petal.Width), dweight = iris$Petal.Width,
#          const = rep(1, 4), data = NULL,
#          method = "GEC",entropy = 1)
# 
# GECal:::GEcalib(~iris$Sepal.Length + iris$Sepal.Width, dweight = iris$Petal.Width,
#          const = rep(1, 3),
#          method = "DS",entropy = 1)
# 
# 
# GECal:::GEcalib(~iris$Sepal.Length + Sepal.Width, dweight = iris$Petal.Width,
#                  data = iris, const = rep(1, 3),
#                  method = "DS",entropy = 1)
# 
# 
# Xs=cbind(
#   c(1,1,1,1,1,1,1,1,1,1),
#   c(1,1,1,1,1,0,0,0,0,0),
#   c(1,3,5,7,9,6,7,8,9,10)
# )
# 
# y=rnorm(10)
# # inclusion probabilities
# piks=1:10/10; piks=piks/ sum(piks); d=1/piks
# # vector of population totals
# total=c(160,124,700); colSums(Xs * d)
# 
# # Calibration weights
# sampling::calib(Xs, total, d = d, method="raking") * d
# calibration <- GECal:::GEcalib(~ 0 + Xs, dweight = d, const = total, method = "DS", entropy = "ET")
# GECal::estimate(y ~ 1, calibration = calibration)
# 
# 
# GECal:::GEcalib(~ 0 + Xs, dweight = d, const = total, method = "GEC0", entropy = "ET")
# 
# GECal:::GEcalib(~ 0 + Xs + g(d), dweight = d, const = c(total, NA), method = "GEC", entropy = "ET",
#                  K_alpha = identity)
# 
# GECal:::GEcalib(~ 0 + Xs + g(d), dweight = d, const = c(total, NA), method = "GEC", entropy = "ET")
# 
# GECal:::GEcalib(~ 0 + g(d) + Xs, dweight = d, const = c(sum(g(d, entropy = 0) * d), total),
#                  method = "GEC", entropy = "ET")
# 
# GECal:::GEcalib(~ 0 + Xs + g(d), dweight = d, const = c(total, sum(g(d, entropy = 0) * d)),
#                  method = "GEC", entropy = "ET")
# 
# GECal:::GEcalib(~ 0 + Xs + g(d), dweight = d, const = c(total, sum(g(d, entropy = -1/2) * d) + 1),
#                  method = "GEC", entropy = -1/2)
# 
# # colSums(cbind(Xs, g(d, entropy = 0)) * GECal:::GEcalib(~ 0 + Xs + g(d), dweight = d, const = c(total, sum(g(d, entropy = 0) * d)),
# #                               method = "GEC", entropy = "ET"))
# 
# GECal:::GEcalib(~ g(d), dweight = d, const = c(150, sum(g(d, entropy = -1) * d)),
#                  method = "GEC", entropy = "EL")
# 
# GECal:::GEcalib(~ 0 + Xs, dweight = d, const = total, method = "DS", entropy = "EL")
# GECal:::GEcalib(~ 0 + Xs, dweight = d, const = total, method = "DS", entropy = "SL")
# 
# GECal:::GEcalib(~ 0 + Xs, dweight = d, const = total, method = "DS", entropy = -0.5)
# 
# sapply(seq(-10, 10, length = 100),
#        function(x) {
#          sum(GECal:::GEcalib(~ 0 + Xs, dweight = d, const = total, method = "DS", entropy = x));
#        })
# 
# sapply(c(seq(-8, 8, length = 100)),
#        function(x) {
#          sum(GECal:::GEcalib(~ g(d), dweight = d,
#                               const = c(150, sum(g(d, entropy = x) * d)), method = "GEC", entropy = x));
#        })
# 
# GECal:::GEcalib(~ g(d), dweight = d,
#                  const = c(150, sum(g(d, entropy = "CE") * d)), method = "GEC", entropy = "CE")
# 
# GECal:::GEcalib(~ g(2 * d), dweight = 2 * d,
#                  const = c(150, sum(g(d, entropy = "CE") * d) * 2), method = "GEC", entropy = "CE")
# 
# GECal:::GEcalib(~ g(d), dweight = d,
#                  const = c(150, sum(g(d, entropy = "PH", del = quantile(d, 0.5)) * d)),
#                  method = "GEC", entropy = "PH", del = quantile(d, 0.5))
# 
# GECal:::GEcalib(~ g(d), dweight = d,
#                  const = c(150, sum(g(d, entropy = "PH", del = quantile(d, 1)) * d)),
#                  method = "GEC", entropy = "PH", del = quantile(d, 1))
# 
# GECal:::GEcalib(~ g(d), dweight = d,
#                  const = c(150, NA), method = "GEC", entropy = "SL")
# 
# 
# GECal:::GEcalib(~ g(d), dweight = d, G.scale = c(rep(c(1,2), 5)),
#                  const = c(150, NA), method = "GEC", entropy = "CE")
# 
# 
# 
# attr(attr(model.frame(~ 0 + Sepal.Length + 1, iris), "terms"), "intercept")
# model.matrix(~ 0 + Sepal.Length + 1,iris)
# head(iris)
# 
# attributes(model.frame(~ 0 + Sepal.Length + 1, iris))
# 
# 
# library(survey)
# data(api)
# 
# dsrs<-svydesign(id=~1, weights=~pw, data=apisrs)
# 
# pop.totals<-c(`(Intercept)`=6194, stypeH=755, stypeM=1018)
# 
# # tmplm <- lm(api00~stype, data = apisrs)
# # head(cbind(tmplm$fitted.values, apisrs$api00))
# 
# (dsrsg<-calibrate(dsrs, ~stype, pop.totals))
# svytotal(~api00, dsrsg)
# # str(unclass(dsrs))
# 
# calibration <- GECal:::GEcalib(~ stype, dweight = pw, data = apisrs, 
#                                const = pop.totals, method = "DS", entropy = "SL")
# 
# n = nrow(apisrs); N = sum(apisrs$pw)
# 
# # pimat <- matrix(n * (n-1) / N / (N-1) - (n / N)^2, nrow = n, ncol = n)
# # diag(pimat) = n / N - (n / N)^2
# pimat <- matrix(0, nrow = n, ncol = n)
# diag(pimat) = n / N
# GECal::estimate(api00 ~ 1, data = apisrs, calibration = calibration,
#                 pimat = pimat * N^2)
# 
