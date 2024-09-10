#' @title Generalized Entropy Calibration
#' 
#' @description
#' \code{GEcalib} computes the generalized entropy calibration weight proposed by Kwon et.al.(2024).
#' The \code{GEcalib} weight minimizes the negative generalized entropy:
#' \deqn{\sum_{i \in A} G(\omega_i)}
#' subject to the calibration constraints \eqn{\sum_{i \in A} \omega_i \bm z_i = \sum_{i \in U} \bm z_i},
#' where \eqn{A} is the index of sample, \eqn{U} is the index of population, and
#' \eqn{\bm z_i} are the auxiliary variables whose population totals are known.
#' 
#' @import nleqslv
#' @importFrom sampling calib
#'
#' @param Xs A matrix of auxiliary variables.
#' @param total A vector of population totals.
#' @param d An optional vector of initial weights.
#' @param entropy An entropy function(Squared-loss:"SL", Empirical Likelihood: "EL", Exponential Tilting: "ET", Cross-Entropy: "CE", Helinger Distance: "HD", "PH")
#' @param method an optional vector that can be used in \code{nleqslv}.
#' @param control an optional list that can be used in \code{nleqslv}.
#' 
#' @return A vector of calibration weights.
#' 
#' @references
#' Kwon, Y., Kim, J., & Qiu, Y. (2024). Debiased calibration estimation using generalized entropy in survey sampling.
#' Arxiv preprint <https://arxiv.org/abs/2404.01076>
#' 
#' Deville, J. C., and Särndal, C. E. (1992). Calibration estimators in survey sampling.
#' Journal of the American statistical Association, 87(418), 376-382.
#' 
#' @examples
#' Xs=cbind(
#'   c(1,1,1,1,1,1,1,1,1,1),
#'   c(1,1,1,1,1,0,0,0,0,0),
#'   c(1,2,3,4,5,6,7,8,9,10)
#' )
#' # inclusion probabilities
#' piks=rep(0.2,times=10); d=1/piks
#' # vector of population totals
#' total=c(50,24,290)
#' 
#' # Calibration weights
#' # calib(Xs, total, d = d, method="raking") * d
#' GEcalib(Xs, total, d, entropy = "ET")
#' GEcalib(Xs, total, entropy = "ET")
#' 
#' # calib(Xs, total, d = d, method="linear") * d
#' GEcalib(Xs, total, d, entropy = "SL")
#' GEcalib(Xs, total, entropy = "SL")

#' @export
GEcalib = function(Xs, total, d = NULL, entropy = c("SL", "EL", "ET", "CE", "HD", "PH"),
                   method = "Newton", control = list(maxit = 1e5, allowSingular = T)){
  
  del = quantile(d, 0.80) 
  
  init = rep(0, length(total))
  if(is.null(d)){
    init[length(init)] = 1
    # init = c(0, 0, 1)
  }else{
    if(entropy == "EL"){
      init[1] = -1
    }else if(entropy == "ET"){
      # init[1] = 0
    }else if(entropy == "SL"){
      init[1] = 1
    }else if(entropy == "CE"){
      init[1] = -1
    }else if(entropy == "HD"){
      init[1] = -2
    }else if(entropy == "PH"){
      init[1] = 1 / sqrt(1 + 1 / del^2)
    }
  }
  
  if(is.null(d)) d = rep(1, nrow(Xs))
  
  nleqslv_res = nleqslv(init, f, jac = h, d = d, Xs = Xs, 
                        total = total, entropy = entropy, del = del,
                        method = method, control = control)
  if(nleqslv_res$termcd != 1){
    if(max(abs(f(nleqslv_res$x, d = d, Xs = Xs, total = total, entropy = entropy, del = del))) > 1e-5)
      print(nleqslv_res)
      w = NA
  }else{
    w = f(nleqslv_res$x, d = d, Xs = Xs, total = total, entropy = entropy, del = del, returnw = T)             
  }
  return(w)
}
