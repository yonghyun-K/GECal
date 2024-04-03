#' Generalized Entropy Calibration

#' @import nleqslv

f = function(lambda, d, Xs, total, entropy, del, ..., returnw = F){
  # w = d * ginv(drop(Xs %*% lambda), entropy = entropy)
  if(entropy == "HD" & any(Xs %*% lambda >= 0)) return(rep(Inf, length(lambda)))
  if(entropy == "PH" & any(abs(Xs %*% lambda) >= del)) return(rep(Inf, length(lambda)))
  
  if(entropy == "SL"){
    w = d * drop(Xs %*% lambda)
  }else if(entropy == "EL"){
    w = -d / drop(Xs %*% lambda)
  }else if(entropy == "ET"){
    w = d * exp(drop(Xs %*% lambda))
  }else if(entropy == "CE"){
    w = d / (1 - exp(drop(Xs %*% lambda)))
  }else if(entropy == "HD"){
    w = d / drop(Xs %*% lambda)^2
  }else if(entropy == "PH"){
    w = d / sqrt(1 / drop(Xs %*% lambda)^2 - 1 / del^2)
  }
  if(entropy != "SL" & any(w <= 0)) return(rep(Inf, length(lambda)))
  if(entropy == "CE" & any(w <= 1)) return(rep(Inf, length(lambda)))
  if(returnw == T){
    return(w)
  }else{
    return(colSums(Xs * w) - total)
  }
}

h = function(lambda, d, Xs, total, entropy, del){
  # return(t(Xs) %*% (Xs * d * ginvprime(drop(Xs %*% lambda), entropy = entropy)))
  if(entropy == "SL"){
    return(t(Xs) %*% (Xs * d))
  }else if(entropy == "EL"){
    w = -d / drop(Xs %*% lambda)
    return(t(Xs) %*% (Xs * (w^2 / d)))
  }else if(entropy == "ET"){
    w = d * exp(drop(Xs %*% lambda))
    return(t(Xs) %*% (Xs * w))
  }else if(entropy == "CE"){
    p_Stmp = 1 / (1 - exp(drop(Xs %*% lambda)))
    return(t(Xs) %*% (Xs * d * p_Stmp * (p_Stmp - 1)))
  }else if(entropy == "HD"){
    return(t(Xs) %*% (Xs * (-2 * d / drop(Xs %*% lambda)^3)))
  }else if(entropy == "PH"){
    return(t(Xs) %*% (Xs * (d * (1 - (drop(Xs %*% lambda) / del)^2)^(-1.5))))
  }
}

#' @export
GEcalib = function(Xs, d, total, entropy = c("SL", "EL", "ET", "CE", "HD", "PH"), ...,
                   DS = T, method = "Newton", control = list(maxit = 1e5, allowSingular = T)){
  
  del = quantile(d, 0.80) 
  
  init = rep(0, length(total))
  if(!DS){
    init[length(init)] = 1
    # init = c(0, 0, 1)
  }else{
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
  
  nleqslv_res = nleqslv(init, f, jac = h, d = d, Xs = Xs, 
                        total = total, entropy = entropy, del = del,
                        method = method, control = control)
  if(nleqslv_res$termcd != 1){
    if(max(abs(f(nleqslv_res$x, d = d, Xs = Xs, total = total, entropy = entropy, del = del))) > 1e-5)
      print(nleqslv_res)
      w = NA
  }else{
    w = f(nleqslv_res$x, d = d, Xs = Xs, total = total, entropy = entropy, del = del, returnw = T)             
  }
  return(w)
}
