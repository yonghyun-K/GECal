h = function(lambda, d, Xs, const, entropy, del, intercept){
  return(t(Xs) %*% (Xs * (d * gprimeinv(drop(Xs %*% lambda), entropy = entropy,
  del = del, intercept = intercept))))
}