targetftn = function(W, d_S, Z_S, Z_St, init, Zbar, entropy, cal, del,..., returnw = F){
  # d_S = rep(1, n) / n
  # Z_S = cbind(1, x_S, log(pi_S)); 
  # Z_St = t(Z_S)
  # init = c(log(n / N), 0, -1)
  # Zbar = c(1, colMeans(x), W)
  if(entropy == "ET" | entropy == "SL" | entropy == "PH"){
    if(W < 0) return(.Machine$double.xmax)
  }else if(entropy == "EL" | entropy == "CE" | entropy == "HD"){
    if(W > 0) return(.Machine$double.xmax)
  }
  alphaHT = Zbar[length(Zbar)]
  Zbar[length(Zbar)] <- W
  nleqslv_res = nleqslv(init, f, jac = h, d_S = d_S, Z_S = Z_S, 
                        Z_St = Z_St, Zbar = Zbar, entropy = entropy, del = del,
                        method = "Newton", control = list(maxit = 1e5, allowSingular = T))
  # if(nleqslv_res$termcd != 1 & nleqslv_res$termcd != 2){
  if(nleqslv_res$termcd != 1){
    if(max(abs(f(nleqslv_res$x, d_S = d_S, Z_S = Z_S, Zbar = Zbar, entropy = entropy, del = del))) > 1e-5)
      return(.Machine$double.xmax)
  }
  w_S = f(nleqslv_res$x, d_S = d_S, Z_S = Z_S, Zbar = Zbar, entropy = entropy, del = del, returnw = T)
  # if(any(is.infinite(w_S))) return(.Machine$double.xmax)
  
  if(returnw == F){
    if(cal == "GEC1") return(sum(G(w_S, entropy = entropy, del = del)) - N * W)
    # else if(cal == "GEC2") return(sum(G(w_S, entropy = entropy)) - N * alphaHT * log(abs(W)))
    else if(cal == "GEC2") return(sum(G(w_S, entropy = entropy, del = del)) - N * (alphaHT + 1) * log(abs(W + 1)))
    # else if(cal == "GEC2") return(sum(G(w_S, entropy = entropy)) - N / g(alphaHT + 1, entropy = entropy) * G(abs(W + 1), entropy = entropy)) # Not working well
  }else{
    return(w_S)
  }
  
  # if(returnw == F){
  #   if(entropy == "EL"){
  #     return(-sum(log(w_S)) - N * W)
  #     # return(-sum(log(w_S)) - N^3 / sum(d_S0)^2 * W)
  #     # return(-sum(log(w_S)) - sum(d_S0) * W)
  #     # return(-sum(log(w_S)) - sum(Z_S[,2] * d_S0) / Zbar[2] * W)
  #     # return(-sum(log(w_S)))
  #     # return(-sum(log(w_S)) + n  * log(-W))
  #     # if(W > -1){
  #     #   return(-sum(log(w_S)) - (N - n)  * log(1 + W))
  #     # }else{
  #     #   return(.Machine$double.xmax)
  #     # }
  #   }else if(entropy == "SL"){
  #     return(sum(w_S^2) / 2 - N * W)
  #     # return(sum(w_S^2) / 2 - sum(d_S0) * W)
  #     # return(sum(w_S^2) / 2 - N^2 / sum(d_S0^2) * W^2 / 2)
  #     # if(W > -1){
  #     #   return(sum(w_S^2) / 2 - N * (1 + sum((d_S0)^2) / N) * log(1 + W))
  #     # }else{
  #     #   return(.Machine$double.xmax)
  #     # }
  #   }else if(entropy == "ET"){
  #     return(sum((w_S) * (log(w_S) - 1)) - N * W)
  #     # return(sum((w_S) * (log(w_S) - 1)) - N * (1 + sum(d_S0 * log(d_S0)) / N) * log(1 + W))
  #     # if(W > 0){
  #     #   return(sum((w_S) * (log(w_S) - 1)) - N * (1 + sum(d_S0 * log(d_S0)) / N) * log(1 + W))
  #     # }else{
  #     #   return(.Machine$double.xmax)
  #     # }
  #   }
  #   else{
  #     return(-sum(log(w_S)))
  #   }
  # }else{
  #   return(w_S)
  # }
}