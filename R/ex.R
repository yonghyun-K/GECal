# devtools::build_manual
# devtools::check
# 
# 
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

# devtools::check(manual = T, document = T)
# GECal::GEcalib()
# lm
# sampling::calib()
# 
# survey::svydesign()
# 
# survey::calibrate()
# 
# 
# N = 10000
# x = data.frame(x1 = rnorm(N, 2, 1), x2= runif(N, 0, 4))
# pi = pt((-x[,1] / 2 - x[,2] / 2), 3);
# pi = ifelse(pi >.7, .7, pi)
# 
# delta = rbinom(N, 1, pi)
# Index_S = (delta == 1)
# n = sum(Index_S); #print(n)
# pi_S = pi[Index_S]; d_S = 1 / pi_S
# x_S = x[Index_S,,drop = F]
# # pimat_S = diag(d_S^2 - d_S) / N^2 # 1 / pi_i * (1 - 1 / pi_i)
# 
# w <- GECal::GEcalib(~ ., dweight = d_S, data = x_S, 
#                     const = colSums(cbind(1, x)),
#                     entropy = "ET", method = "DS")$w
# 
# w <- GECal::GEcalib(~ ., dweight = d_S, data = x_S, 
#                     const = colSums(cbind(1, x)),
#                     entropy = "ET", method = "GEC0")$w
# 
# w <- GECal::GEcalib(~ . + g(d_S), dweight = d_S, data = x_S, 
#                     const = colSums(cbind(1, x, log(1 / pi))),
#                     entropy = "ET", method = "GEC")$w
# 
# w <- GECal::GEcalib(~ . + g(d_S), dweight = d_S, data = x_S, 
#                     const = colSums(cbind(1, x, NA)),
#                     entropy = "ET", method = "GEC")$w
# 
# tmp <- GECal::GEcalib(~ . + g(d_S), dweight = d_S, data = x_S, 
#                       const = colSums(cbind(1, x, NA)),
#                       entropy = "ET", method = "GEC")
# names(tmp)
# 
# 
# e = rnorm(N, 0, 1)
# y = x[,1] + x[,2] + e;
# y_S = y[Index_S] # plot(x_S, y_S)
# # data_S = cbind(pi, data)[Index_S,]

