gprime1 = function(x, entropy, del){ # Inverse of gprime(d)
  switch(entropy,
         SL = rep(1, length(x)),
         EL = x^2,
         ET = x,
         CE = x * (x - 1),
         HD = 2 * x^(1.5),
         Hb = ifelse(abs(x) < del, 1, NA), # ???
         PH = (1 + (x / del)^2)^(1.5))
}