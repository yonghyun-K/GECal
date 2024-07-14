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