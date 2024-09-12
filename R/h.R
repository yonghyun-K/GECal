h = function(lambda, d, Xs, const, entropy, del, weight.scale, G.scale, intercept){
  return(t(Xs) %*% (Xs * (d  / weight.scale^2  / G.scale * gprimeinv(
    drop(Xs %*% lambda / weight.scale / G.scale), 
  entropy = entropy,
  del = del, intercept = intercept))))
}