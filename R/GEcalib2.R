GEcalib2 = function(formula, dweight, data = NULL, const, method = c("GEC", "GEC0", "DS"),
  entropy = c("SL", "EL", "ET", "CE", "HD", "PH"), 
  weight.bound = NULL, weight.scale = 1,
  opt.method = c("nleqslv", "optim", "CVXR"),
  del = quantile(dweight, 0.75),
  K_alpha = NULL, is.total = T
){
  entropy <- if (is.numeric(entropy)) {
    entropy  # Assign entropy directly if it's numeric
  } else {
    switch(entropy,
           "EL" = -1, "ET" = 0, "SL" = 1, "HD" = -1 / 2,
           "CE" = "CE", "PH" = "PH",
           stop("Invalid entropy value")  # To handle unexpected values
    )
  }

  if (is.null(data)) {
    assign("entropy", entropy, envir = sys.frame())
    assign("del", del, envir = sys.frame())
    environment(g) <- environment(); environment(formula) <- environment()
    mf <- model.frame(formula, parent.frame())  # Evaluate in parent environment
    dweight0 = dweight
  } else {    
    assign("entropy", entropy, envir = environment())
    assign("del", del, envir = environment())
    environment(g) <- environment(); environment(formula) <- environment()
    mf <- model.frame(formula, data)  # Evaluate in the provided data
    dweight0 = data$dweight
  }

  Xs <- model.matrix(attr(mf, "terms"), mf)
  if(rcond(Xs) < .Machine$double.eps){
    stop("Error: rcond(model.matrix(formula, data)) < .Machine$double.eps")
  }
  
  if(length(const) != ncol(Xs)){
    stop("Error: length(const) != ncol(model.matrix(formula, data))")
  }
  
  # Get the name of the transformed column g(dweight)
  transformed_name <- paste0("g(", deparse(substitute(dweight)), ")")
  
  # Find the location of the column in the model matrix
  col_position <- which(colnames(Xs) == transformed_name)
  
  if(method == "GEC" & length(col_position) == 0){
    stop("Method GEC needs g(dweight) in the formula")
  } 
  
  What = sum(g(dweight0, entropy = entropy, del = del) * dweight0)
  if(is.total & attr(attr(mf, "terms"), "intercept") == 1 & !is.null(K_alpha)){
    N = const[1]
    K_alpha = function(x) (What + N) * log(abs(x / N + 1))
  }else{
    K_alpha = identity
  }

  init = rep(0, length(const))
  if(method == "GEC"){
    init[col_position] = 1
    d = rep(1, nrow(Xs))
    intercept = rep(0, length(d))
  }else if(method == "DS"){
    d = dweight0
    if(entropy == 0){
      intercept = rep(0, length(d))      
    }else{
      intercept = rep(1, length(d))
    }
  }else if(method == "GEC0"){
    d = rep(1, nrow(Xs))
    intercept = g(dweight0, entropy = entropy, del = del)
  }
  
  if(any(is.na(const)) & method == "GEC"){
    if(all(which(is.na(const)) == col_position)){
      nlmres= nlm(targetftn, p = What, d = d, Xs = Xs, init = init,
                  const = const, entropy = entropy, del = del,
                  intercept = intercept, K_alpha = K_alpha)
      # if(nlmres$code != 1 & nlmres$code != 2 & nlmres$code != 3) print(nlmres$estimate)
      if(nlmres$code != 1 & nlmres$code != 2 & nlmres$code != 3) stop(nlmres$code)
      W = nlmres$estimate
      if(nlmres$minimum >= .Machine$double.xmax){
        print(nlmres)
        w = NA
      }else{
        w = targetftn(W, d = d, Xs = Xs, init = init,
                        const = const, entropy = entropy, del = del,
                        intercept = intercept, K_alpha = K_alpha, returnw = T)
      }
    }else{
      stop("NA appears in const outside g(d)")
    }
  }else{
    nleqslv_res = nleqslv::nleqslv(init, f, jac = h, d = d, Xs = Xs, 
                                   const = const, entropy = entropy, del = del,
                                   intercept = intercept, control = list(maxit = 1e5, allowSingular = T),
                                   xscalm = "auto")
    # control = control
    
    if(nleqslv_res$termcd != 1){
      tmpval <- f(nleqslv_res$x, d = d, Xs = Xs, 
                  const = const, entropy = entropy, del = del,
                  intercept = intercept)
      # return(f(init, d = d, Xs = Xs, 
      #          const = const, entropy = entropy, del = del,
      #          intercept = intercept))
      
      if(any(is.nan(tmpval)) | (max(abs(tmpval)) > 1e-5)){
        w = NA      
      }else{
        w = NULL
      }
    }else{
      w = NULL  
    }
    if(is.null(w)){
      w = f(nleqslv_res$x, d = d, Xs = Xs, 
            const = const, entropy = entropy, del = del,
            intercept = intercept, returnw = T)    
    }
  }
  return(w)
}


# GEcalib2(~Sepal.Length + Sepal.Width + g(Petal.Width), dweight = Petal.Width,
#          data = iris, const = rep(1, 4),
#          method = "GEC",entropy = 1)
# 
# GEcalib2(~Sepal.Length + Sepal.Width, dweight = Petal.Width,
#          data = iris, const = rep(1, 3),
#          method = "DS",entropy = 1)
# 
# GEcalib2(~iris$Sepal.Length + iris$Sepal.Width + g(iris$Petal.Width), dweight = iris$Petal.Width,
#          const = rep(1, 4), data = NULL,
#          method = "GEC",entropy = 1)
# 
# GEcalib2(~iris$Sepal.Length + iris$Sepal.Width, dweight = iris$Petal.Width,
#          const = rep(1, 3),
#          method = "DS",entropy = 1)

# Xs=cbind(
#   c(1,1,1,1,1,1,1,1,1,1),
#   c(1,1,1,1,1,0,0,0,0,0),
#   c(1,3,5,7,9,6,7,8,9,10)
# )
# # inclusion probabilities
# piks=1:10/10; piks=piks/ sum(piks); d=1/piks
# # vector of population totals
# total=c(160,124,700); colSums(Xs * d)
# 
# # Calibration weights
# sampling::calib(Xs, total, d = d, method="raking") * d
# GECal:::GEcalib2(~ 0 + Xs, dweight = d, const = total, method = "DS", entropy = "ET")
# GECal:::GEcalib2(~ 0 + Xs, dweight = d, const = total, method = "GEC0", entropy = "ET")
# 
# GECal:::GEcalib2(~ 0 + Xs + g(d), dweight = d, const = c(total, NA), method = "GEC", entropy = "ET",
#                  K_alpha = identity)
# 
# GECal:::GEcalib2(~ 0 + Xs + g(d), dweight = d, const = c(total, NA), method = "GEC", entropy = "ET")
# 
# GECal:::GEcalib2(~ 0 + g(d) + Xs, dweight = d, const = c(sum(g(d, entropy = 0) * d), total),
#                  method = "GEC", entropy = "ET")
# 
# GECal:::GEcalib2(~ 0 + Xs + g(d), dweight = d, const = c(total, sum(g(d, entropy = 0) * d)),
#                  method = "GEC", entropy = "ET")
# 
# GECal:::GEcalib2(~ 0 + Xs + g(d), dweight = d, const = c(total, sum(g(d, entropy = -1/2) * d) + 1),
#                  method = "GEC", entropy = -1/2)
# 
# # colSums(cbind(Xs, g(d, entropy = 0)) * GECal:::GEcalib2(~ 0 + Xs + g(d), dweight = d, const = c(total, sum(g(d, entropy = 0) * d)),
# #                               method = "GEC", entropy = "ET"))
# 
# GECal:::GEcalib2(~ g(d), dweight = d, const = c(150, sum(g(d, entropy = -1) * d)),
#                  method = "GEC", entropy = "EL")
# 
# GECal:::GEcalib2(~ 0 + Xs, dweight = d, const = total, method = "DS", entropy = "EL")
# GECal:::GEcalib2(~ 0 + Xs, dweight = d, const = total, method = "DS", entropy = "SL")
# 
# GECal:::GEcalib2(~ 0 + Xs, dweight = d, const = total, method = "DS", entropy = -0.5)
# 
# sapply(seq(-10, 10, length = 100),
#        function(x) {
#          sum(GECal:::GEcalib2(~ 0 + Xs, dweight = d, const = total, method = "DS", entropy = x));
#        })
# 
# sapply(c(seq(-8, 8, length = 100)),
#        function(x) {
#          sum(GECal:::GEcalib2(~ g(d), dweight = d,
#                               const = c(150, sum(g(d, entropy = x) * d)), method = "GEC", entropy = x));
#        })
# 
# GECal:::GEcalib2(~ g(d), dweight = d,
#                  const = c(150, sum(g(d, entropy = "CE") * d)), method = "GEC", entropy = "CE")
# 
# GECal:::GEcalib2(~ g(d), dweight = d,
#                  const = c(150, sum(g(d, entropy = "PH", del = quantile(d, 0.5)) * d)), 
#                  method = "GEC", entropy = "PH", del = quantile(d, 0.5))
# 
# GECal:::GEcalib2(~ g(d), dweight = d,
#                  const = c(150, sum(g(d, entropy = "PH", del = quantile(d, 1)) * d)), 
#                  method = "GEC", entropy = "PH", del = quantile(d, 1))
# 
# GECal:::GEcalib2(~ g(d), dweight = d,
#                  const = c(150, NA), method = "GEC", entropy = "SL")
# 
# 
# GECal:::GEcalib2(~ g(d), dweight = d,
#                  const = c(150, NA), method = "GEC", entropy = "CE")
# 
# 
# attr(attr(model.frame(~ 0 + Sepal.Length + 1, iris), "terms"), "intercept")
# model.matrix(~ 0 + Sepal.Length + 1,iris)
# head(iris)
# 
# attributes(model.frame(~ 0 + Sepal.Length + 1, iris))
# 



