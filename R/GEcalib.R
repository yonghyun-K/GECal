f = function(lambda, d_S, Z_S, Zbar, type, del, ..., returnw = F){
  # w_S = d_S * ginv(drop(Z_S %*% lambda), type = type)
  if(type == "HD" & any(Z_S %*% lambda >= 0)) return(rep(Inf, length(lambda)))
  if(type == "PH" & any(abs(Z_S %*% lambda) >= del)) return(rep(Inf, length(lambda)))
  
  if(type == "SL"){
    w_S = d_S * drop(Z_S %*% lambda)
  }else if(type == "EL"){
    w_S = -d_S / drop(Z_S %*% lambda)
  }else if(type == "ET"){
    w_S = d_S * exp(drop(Z_S %*% lambda))
  }else if(type == "CE"){
    w_S = d_S / (1 - exp(drop(Z_S %*% lambda)))
  }else if(type == "HD"){
    w_S = d_S / drop(Z_S %*% lambda)^2
  }else if(type == "PH"){
    w_S = d_S / sqrt(1 / drop(Z_S %*% lambda)^2 - 1 / del^2)
  }
  if(type != "SL" & any(w_S <= 0)) return(rep(Inf, length(lambda)))
  if(type == "CE" & any(w_S <= 1)) return(rep(Inf, length(lambda)))
  if(returnw == T){
    return(w_S)
  }else{
    return(colSums(Z_S * w_S) / N - Zbar)
  }
}

h = function(lambda, d_S, Z_S, Z_St, Zbar, type, del){
  # return(Z_St %*% (Z_S * d_S * ginvprime(drop(Z_S %*% lambda), type = type) / N))
  if(type == "SL"){
    return(Z_St %*% (Z_S * d_S) / N)
  }else if(type == "EL"){
    w_S = -d_S / drop(Z_S %*% lambda)
    return(Z_St %*% (Z_S * (w_S^2 / d_S)) / N)
  }else if(type == "ET"){
    w_S = d_S * exp(drop(Z_S %*% lambda))
    return(Z_St %*% (Z_S * w_S) / N)
  }else if(type == "CE"){
    p_Stmp = 1 / (1 - exp(drop(Z_S %*% lambda)))
    return(Z_St %*% (Z_S * d_S * p_Stmp * (p_Stmp - 1)) / N)
  }else if(type == "HD"){
    return(Z_St %*% (Z_S * (-2 * d_S / drop(Z_S %*% lambda)^3)) / N)
  }else if(type == "PH"){
    return(Z_St %*% (Z_S * (d_S * (1 - (drop(Z_S %*% lambda) / del)^2)^(-1.5))) / N)
  }
}

GEcalib = function(d_S, Z_S, Zbar, type, del, ...){
  nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, Z_S = Z_S, 
                        Z_St = Z_St, Zbar = Zbar, type = type, del = del,
                        method = "Newton", control = list(maxit = 1e5, allowSingular = T))
  # print(type); print(cal); print(nleqslv_res)
  # if(cal == "GEC" & type == "ET") print(nleqslv_res$x)
  if(nleqslv_res$termcd != 1){
    if(max(abs(f(nleqslv_res$x, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del))) > 1e-5)
      w_S = NA
    # stop(c(paste(type, cal, sep = "_"), nleqslv_res))
  }else{
    # if(type == "CE" & cal == "DS") stop(c(paste(type, cal, sep = "_"), nleqslv_res))
    w_S = f(nleqslv_res$x, d_S = d_S, Z_S = Z_S, Zbar = Zbar, type = type, del = del, returnw = T)             
  }
  return(w_S)
}