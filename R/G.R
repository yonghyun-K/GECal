G = function(x, entropy, del = NULL){
  if (is.numeric(entropy)) {
    if(!(entropy > 0 & entropy %% 2 == 1) & any(x < 0)) return(rep(Inf, length(x)))
    if(entropy == 0){
      x * (log(x) - 1)
    }else if(entropy == -1){
      -log(x)
    }else{
      x^(entropy + 1) / (entropy * (entrop + 1))
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
      (x^(entropy + 1) - (entrop + 1) * x + entropy) / (entropy * (entrop + 1))
    }
  }else{
    return(rep(Inf, length(x)))
  }
}

ginv = function(x, entropy = NULL, del = NULL, intercept = rep(0, length(x))){
  if(is.null(entropy))  entropy <- get("entropy", envir = environment())
  if (is.numeric(entropy)) {
    if(!(entropy %% 2 == 1) & any(intercept + x * entropy < 0)) return(rep(Inf, length(x)))
    if(entropy == 0){
      w = exp(x)
    }else{
      w = (intercept + x * entropy)^(1 / entropy)
    }
    # if(!(entropy > 0 & entropy %% 2 == 1) & any(w < 0)) return(rep(Inf, length(x)))
  } else {
    if(entropy == "PH" & any(abs(x) >= del)) return(rep(Inf, length(x)))
    w = switch(entropy,
           CE = 1 / (1 - exp(x)),
           PH = 1 / sqrt(1 / x^2 - 1 / del^2))
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
      w = (intercept + x * entropy)^(1 / entropy - 1)
    }
    # if(any(w < 0)) return(rep(Inf, length(x)))
  } else {
    if(entropy == "PH" & any(abs(x) >= del)) return(rep(Inf, length(x)))
    w = switch(entropy,
               CE = exp(x) / (1 - exp(x))^2,
               PH = (1 - (x / del)^2)^(-1.5))
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

#' @export
g = function(x, entropy = NULL, del = NULL){
  if(is.null(entropy))  entropy <- get("entropy", envir = parent.frame())

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
