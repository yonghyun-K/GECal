GEcalib2 = function(formula, dweight, data, const, method = c("GEC", "DS"),
  entropy = c("SL", "EL", "ET", "CE", "HD", "PH"), 
  weight.bound = NULL, weight.scale = 1,
  opt.method = c("nleqslv", "optim", "CVXR"),
  del = quantile(dweight, 0.75)
){
  assign("entropy", entropy, envir = environment())
  environment(g) <- environment()
  environment(formula) <- environment()
  Z <- model.matrix(formula, data)
  if(length(const) != ncol(Z)){
    stop("Error: length(const) != ncol(Z)")
  }
  
  order <- if (is.numeric(entropy)) {
    entropy  # Assign entropy directly if it's numeric
  } else {
    switch(entropy,
           "EL" = -1, "ET" = 0, "SL" = 1,
           "CE" = -1, "HD" = -1 / 2, "PH" = 1 / sqrt(1 + 1 / del^2),
           stop("Invalid entropy value")  # To handle unexpected values
    )
  }
  
  init = rep(0, length(const))
  if(method == "GEC"){
    d = rep(1, nrow(Z))
    init[length(init)] = 1
    # init = c(0, 0, 1)
  }else if(method == "DS"){
    d = dweight
    if(entropy == "EL"){
      init[1] = -1
    }else if(entropy == "ET"){
      # init[1] = 0
    }else if(entropy == "SL"){
      init[1] = 1
    }else if(entropy == "CE"){
      init[1] = -1
    }else if(entropy == "HD"){
      init[1] = -1
    }else if(entropy == "PH"){
      init[1] = 1 / sqrt(1 + 1 / del^2)
    }
  }
  
  return(colnames(Z))
}

# nleqslv()
# 
# lm
# 
# library(GECal)
# 
# lm(Sepal.Length ~ . + sqrt(Sepal.Width), data = iris)
# 
# GEcalib2(~Sepal.Length + Sepal.Width + g(Petal.Width), dweight = Petal.Width, data = iris, const = c(1,1),
#          entropy = "SL")
# 
# GEcalib2(~. + g(Petal.Width), dweight = Petal.Width, data = iris, const = c(1,1),
#          entropy = "SL")


