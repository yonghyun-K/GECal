targetftn = function(W, d, Xs, init, const, entropy, del, intercept, K_alpha, 
                     weight.scale, G.scale,..., returnw = F){
  if (is.numeric(entropy)) {
    if(entropy >= 0 & W < 0) return(.Machine$double.xmax)
    else if(entropy < 0 & W > 0) return(.Machine$double.xmax)
  } else {
    if(entropy == "PH" & W < 0) return(.Machine$double.xmax)
    else if(entropy == "CE" & W > 0) return(.Machine$double.xmax)
  }
  
  const[length(const)] <- W
  nleqslv_res = nleqslv(init, f, jac = h, d = d, Xs = Xs, 
                        const = const, entropy = entropy, del = del,
                        weight.scale = weight.scale, G.scale = G.scale,
                        intercept = intercept, 
                        method = "Newton", control = list(maxit = 1e5, allowSingular = T))
  # if(nleqslv_res$termcd != 1 & nleqslv_res$termcd != 2){
  if(nleqslv_res$termcd != 1){
    if(max(abs(f(nleqslv_res$x, d = d, Xs = Xs, 
                 const = const, entropy = entropy, del = del,
                 weight.scale = weight.scale, G.scale = G.scale,
                 intercept = intercept))) > 1e-5)
      return(.Machine$double.xmax)
  }
  w_S = f(nleqslv_res$x, d = d, Xs = Xs, 
          const = const, entropy = entropy, del = del,
          weight.scale = weight.scale, G.scale = G.scale,
          intercept = intercept, returnw = T)
  # if(any(is.infinite(w_S))) return(.Machine$double.xmax)
  
  if(returnw == F){
    # return(sum(G(w_S, type = type, del = del)) - n * W)
    return(sum(G.scale * G(weight.scale * w_S, entropy = entropy, del = del)) - K_alpha(W))
  }else{
    return(w_S)
  }
}