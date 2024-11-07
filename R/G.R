#' @title Debiasing covariate for GECalib
#' 
#' @description
#' It returns the debiasing covariate, which is equivalent to the first order derivatie of
#' the generalized entropy \eqn{G}.
#' 
#' @param x A vector of design weights
#' @param entropy An optional data frame containing the variables in the model (specified by \code{formula}).
#' @param del The optional vector for threshold (\eqn{\delta}) when \code{entropy == "PH"}.
#' 
#' @return A vector of debiasing covariate. 
#' 
#' @examples
#' set.seed(11)
#' N = 10000
#' x = data.frame(x1 = rnorm(N, 2, 1), x2= runif(N, 0, 4))
#' pi = pt((-x[,1] / 2 - x[,2] / 2), 3);
#' pi = ifelse(pi >.7, .7, pi)
#' 
#' g_EL <- g(1 / pi, entropy = 1)
#' g_ET <- g(1 / pi, entropy = 0)
#' g_EL <- g(1 / pi, entropy = -1)
#' 
#' @export
g = function(x, entropy = NULL, del = NULL){
  if(is.null(entropy)){
    if(exists("entropy", envir = sys.frame())){
      entropy <- get("entropy", envir = sys.frame())      
    }else if(exists("entropy", envir = parent.frame())){
      entropy <- get("entropy", envir = parent.frame())
    }else{
      stop("Specify entropy in g function.")      
    }
  }
  if(entropy == "PH" & is.null(del)){
    if(!exists("del", envir = parent.frame())) stop("Specify del in g function.")
    del <- get("del", envir = parent.frame())
  }

  if (is.numeric(entropy)) {
    # if(!(entropy > 0 & entropy %% 2 == 1) & any(x < 0)) return(rep(Inf, length(x)))
    if(entropy == 0){
      log(x)
    }else{
      x^entropy / entropy
    }
  } else {
    switch(entropy,
           CE = log(1 - 1 / x),
           PH = x / sqrt(1 + (x / del)^2))
  }
}


G = function(x, entropy, del = NULL){
  if (is.numeric(entropy)) {
    if(!(entropy > 0 & entropy %% 2 == 1) & any(x < 0)) return(rep(Inf, length(x)))
    if(entropy == 0){
      x * (log(x) - 1)
    }else if(entropy == -1){
      -log(x)
    }else{
      x^(entropy + 1) / (entropy * (entropy + 1))
    }
  } else {
    switch(entropy,
           CE = (x-1) * log(x-1) - x * log(x),
           PH = del^2 * sqrt(1 + (x / del)^2))
  }
}

G_DS = function(x, entropy, del = NULL){
  if (is.numeric(entropy)) {
    if(!(entropy > 0 & entropy %% 2 == 1) & any(x < 0)) return(rep(Inf, length(x)))
    if(entropy == 0){
      x * log(x) - x + 1
    }else if(entropy == -1){
      -log(x) + x - 1
    }else{
      (x^(entropy + 1) - (entropy + 1) * x + entropy) / (entropy * (entropy + 1))
    }
  }else{
    return(rep(Inf, length(x)))
  }
}

ginv = function(x, entropy = NULL, del = NULL, intercept = rep(0, length(x))){
  if(is.null(entropy))  entropy <- get("entropy", envir = environment())
  if (is.numeric(entropy)) {
    if(!(entropy > 0 & entropy %% 2 == 1) & any((intercept + x) *  entropy < 0)) return(rep(Inf, length(x)))
    if(entropy == 0){
      w = exp(intercept + x)
    }else{
      w = ((intercept + x) * entropy)^(1 / entropy)
    }
    # if(!(entropy > 0 & entropy %% 2 == 1) & any(w < 0)) return(rep(Inf, length(x)))
  } else {
    if(entropy == "PH" & any(abs(intercept + x) >= del)) return(rep(Inf, length(x)))
    w = switch(entropy,
               CE = 1 / (1 - exp((intercept + x))),
               PH = 1 / sqrt(1 / (intercept + x)^2 - 1 / del^2))
    if(entropy == "CE" & any(w <= 1)) return(rep(Inf, length(x)))
  }
  return(w)
}

gprimeinv = function(x, entropy = NULL, del = NULL, intercept = rep(0, length(x))){
  if (is.numeric(entropy)) {
    # if(!(entropy %% 2 == 1) & any(intercept + x * entropy < 0)) return(rep(Inf, length(x)))
    if(entropy == 0){
      w = exp(x)
    }else{
      w = ((intercept + x) * entropy)^(1 / entropy - 1)
    }
    # if(any(w < 0)) return(rep(Inf, length(x)))
  } else {
    if(entropy == "PH" & any(abs(intercept + x) >= del)) return(rep(Inf, length(x)))
    w = switch(entropy,
               CE = exp(x) / (1 - exp((intercept + x)))^2,
               PH = (1 - ((intercept + x) / del)^2)^(-1.5))
    if(entropy == "CE" & any(w <= 1)) return(rep(Inf, length(x)))
  }
  return(w)
}

fprime = function(x, entropy = NULL, del = NULL){ # Inverse of gprime(d)
  if (is.numeric(entropy)) {
    if(!(entropy > 0 & entropy %% 2 == 1) & any(x < 0)) return(rep(Inf, length(x)))
    x^(-entropy+1)
  } else {
    switch(entropy,
           CE = x * (x - 1),
           PH = (1 + (x / del)^2)^(1.5))
  }
}