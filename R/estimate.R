#' @export
estimate <- function(formula, data = NULL, calibration, pimat = NULL){
  response_vars <- all.vars(formula[[2]])
  
  # If data is provided, extract the variables from the data
  if (!is.null(data)) {
    # Extract the columns from the data corresponding to the response variables
    ys <- do.call(cbind, lapply(response_vars, function(var) data[[var]]))
  } else {
    # If no data is provided, evaluate the variables in the global environment
    ys <- do.call(cbind, lapply(response_vars, function(var) eval(parse(text = var))))
  }
  
  Xs <- calibration$Xs
  w <- calibration$w
  dweight <- calibration$dweight
  entropy <- calibration$entropy
  del <- calibration$del
  G.scale <- calibration$G.scale
  weight.scale <- calibration$weight.scale
  method <- calibration$method
  const <- calibration$const
  What <- calibration$What
  K_alpha <- calibration$K_alpha
  
  if(is.null(pimat)) pimat = diag(w * (w - 1)) # dweight * (dweight - 1)
  
  if(method == "DS"){
    gammahat = solve(t(Xs) %*% (Xs * dweight / G.scale),
                     t(Xs) %*% (ys * dweight / G.scale))
  }else if(method == "GEC" | method == "GEC0"){
    gammahat = solve(t(Xs) %*% (Xs * fprime(dweight * weight.scale, entropy = entropy, 
                                            del = del) / G.scale / weight.scale^2),
                     t(Xs) %*% (ys * fprime(dweight * weight.scale, entropy = entropy, 
                                            del = del) / G.scale / weight.scale^2))
  }
  
  if(method == "GEC" & any(is.na(const))){
    hatSigmazz = t(Xs) %*% (Xs * fprime(dweight * weight.scale, entropy = entropy, 
                                        del = del) / G.scale / weight.scale^2)
    xcolums = sequence(ncol(Xs))[-col_position]; gcolums = col_position
    hatSigmaxx = hatSigmazz[xcolums, xcolums]
    hatSigmagx = hatSigmazz[gcolums, xcolums, drop = F]
    hatSigmagg = hatSigmazz[gcolums, gcolums, drop = F]
    hatSigmagg_x = drop(hatSigmagg - hatSigmagx %*% solve(hatSigmaxx, t(hatSigmagx)))
    
    if(identical(K_alpha, identity)){
      Xs[,gcolums] <- c(hatSigmagx %*% solve(hatSigmaxx, t(Xs[,xcolums])))
    }else{
      N = const[1]
      Xs[,gcolums] <- (1 / hatSigmagg_x / (1 / (What + N) + 1 / hatSigmagg_x)) * c(hatSigmagx %*% solve(hatSigmaxx, t(Xs[,xcolums])))
    }
  }
  
  yhat = drop(Xs %*% gammahat)
  # print(head(cbind(ys, yhat)))
  Varhat = drop(t(ys - yhat) %*% pimat %*% (ys - yhat))

  return( list(cov = Varhat,
    estimate = cbind("Estimate" = colSums(ys * w), "Std. Error" = drop(sqrt(diag(Varhat, nrow = ncol(ys)))))))
}
