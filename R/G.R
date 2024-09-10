G = function(x, entropy, del){
  switch(entropy,
         SL = x^2/2,
         EL = -log(x),
         ET = x * (log(x) - 1),
         CE = (x-1) * log(x-1) - x * log(x),
         HD = -4 * sqrt(x),
         PH = del^2 * sqrt(1 + (x / del)^2))
}

#' @export
g = function(x, del = NULL){
  entropy <- get("entropy", envir = environment())
  switch(entropy,
         SL = x,
         EL = - 1 / x,
         ET = log(x),
         CE = log(1 - 1 / x),
         HD = -2 / sqrt(x),
         PH = x / sqrt(1 + (x / del)^2))
}
