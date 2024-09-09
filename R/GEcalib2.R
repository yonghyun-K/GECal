GEcalib2 = function(formula, dweight, data, const, method = c("GEC", "DS"),
  entropy = c("SL", "EL", "ET", "CE", "HD", "PH"), 
  weight.bound = NULL, weight.scale = 1,
  opt.method = c("nleqslv", "optim", "CVXR")
){
  assign("entropy", entropy, envir = environment())
  environment(g) <- environment()
  environment(formula) <- environment()
  Z <- model.matrix(formula, data)
  return(colnames(Z))
}

# GEcalib2(~Sepal.Length + Sepal.Width + g(Petal.Width), dweight = Petal.Width, data = iris, const = c(1,1),
#          entropy = "SL")
# 
# 
# library(GECal)

