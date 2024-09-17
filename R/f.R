f = function(lambda, d, Xs, const, entropy, del, weight.scale, G.scale,..., returnw = FALSE, intercept = rep(0, nrow(Xs))){
  w = d / weight.scale * ginv(drop(Xs %*% lambda / G.scale / weight.scale), 
                              entropy = entropy, del = del, intercept = intercept)
  if(returnw == TRUE){
    return(w)
  }else{
    return(colSums(Xs * w) - const)
  }
}