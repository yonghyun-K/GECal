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