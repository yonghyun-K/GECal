f = function(lambda, d, Xs, const, entropy, del,..., returnw = F, intercept = rep(0, nrow(Xs))){
  w = d * ginv(drop(Xs %*% lambda), entropy = entropy, del = del, intercept = intercept)
  if(returnw == T){
    return(w)
  }else{
    return(colSums(Xs * w) - const)
  }
}