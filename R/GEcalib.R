#' @title Generalized Entropy Calibration
#' 
#' @description
#' \code{GEcalib} computes the generalized entropy calibration weight proposed by Kwon et.al.(2024).
#' The \code{GEcalib} weight minimizes the negative generalized entropy:
#' \deqn{\sum_{i \in A} G(\omega_i)}
#' subject to the calibration constraints \eqn{\sum_{i \in A} \omega_i \bm z_i = \sum_{i \in U} \bm z_i},
#' where \eqn{A} is the index of sample, \eqn{U} is the index of population, and
#' \eqn{\bm z_i^T = (\bm x_i^T, g(d_i))} are the auxiliary variables whose population totals are known.
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
#' @section Summary:
#' 
#' \tabular{cc}{
#' \strong{GEC} \tab \strong{DS} \cr
#' \eqn{\min_{\bm \omega} \left(-H(\bm \omega)\right) = \sum_{i \in A}G(\omega_i) \quad} \tab 
#' \eqn{\quad \min_{\bm \omega} D(\bm \omega, \bm d) = \sum_{i \in A}d_iG(\omega_i / d_i)} \cr
#' \deqn{G(\omega) = \begin{cases} \frac{1}{r(r+1)} \omega^{r+1} & r \neq 0, -1\\ 
#' \omega \log \omega - \omega & r = 0\text{(ET)} \\ 
#' -\log \omega & r = -1\text{(EL)} \end{cases}} 
#' \tab \deqn{G(\omega) = \begin{cases} \frac{1}{r(r+1)} \left(\omega^{r+1} - (r+1)\omega + r\right) & r \neq 0, -1 \\
#' \omega \log \omega - \omega + 1 & r = 0\text{(ET)} \\
#' -\log \omega + \omega - 1 & r = -1\text{(EL)} \end{cases}} \cr
#' }
#' 
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
GEcalib = function(formula, dweight, data = NULL, const, 
                    method = c("GEC", "GEC0", "DS"),
                    entropy = c("SL", "EL", "ET", "CE", "HD", "PH"), 
                    # weight.bound = NULL, 
                   weight.scale = 1, G.scale = 1,
                    # opt.method = c("nleqslv", "optim", "CVXR"),
                    del = NULL,
                    K_alpha = NULL, is.total = T
){
  entropy <- if (is.numeric(entropy)) {
    entropy  # Assign entropy directly if it's numeric
  } else {
    switch(entropy,
           "EL" = -1, "ET" = 0, "SL" = 1, "HD" = -1 / 2,
           "CE" = "CE", "PH" = "PH",
           stop("Invalid entropy value")  # To handle unexpected values
    )
  }
  
  environment(g) <- environment(); environment(formula) <- environment()
  
  if (is.null(data)) {
    assign("entropy", entropy, envir = sys.frame())
    assign("del", del, envir = sys.frame())
    mf <- model.frame(formula, parent.frame())  # Evaluate in parent environment
    dweight0 = dweight
  } else {    
    assign("entropy", entropy, envir = environment())
    assign("del", del, envir = environment())
    mf <- model.frame(formula, data)  # Evaluate in the provided data
    dweight0 = eval(substitute(dweight), envir = data)
  }
  
  if(!is.null(del)) del = quantile(dweight0, 0.75)
  
  Xs <- model.matrix(attr(mf, "terms"), mf)
  if(rcond(Xs) < .Machine$double.eps){
    stop("Error: rcond(model.matrix(formula, data)) < .Machine$double.eps")
  }
  
  if(length(const) != ncol(Xs)){
    stop("Error: length(const) != ncol(model.matrix(formula, data))")
  }
  
  # Get the name of the transformed column g(dweight)
  transformed_name <- paste0("g(", deparse(substitute(dweight)), ")")
  
  # Find the location of the column in the model matrix
  col_position <- which(colnames(Xs) == transformed_name)
  
  if(method == "GEC" & length(col_position) == 0){
    stop("Method GEC needs g(dweight) in the formula")
  } 
  
  if(length(weight.scale) == 0) stop("weight.scale has to be positive length")
  
  if(any(weight.scale != 1)){  
    if(any(weight.scale <= 0)) stop("weight.scale has to be positive")
    if(method == "GEC" | method == "GEC0"){
      Xs[,col_position] <- weight.scale * g(dweight0 * weight.scale, entropy = entropy, del = del)
    }else if(method == "DS"){
      warning("weight.scale is not supported in DS method")
    }
  }
  Xs[,col_position] <- Xs[,col_position] * G.scale
  
  What = sum(g(dweight0 * weight.scale, entropy = entropy, del = del) * dweight0 * weight.scale * G.scale)
  if(is.total & attr(attr(mf, "terms"), "intercept") == 1 & !is.null(K_alpha)){
    N = const[1]
    K_alpha = function(x) (What + N) * log(abs(x / N + 1))
  }else{
    K_alpha = identity
  }
  
  init = rep(0, length(const))
  if(method == "GEC"){
    init[col_position] = 1
    d = rep(1, nrow(Xs))
    intercept = rep(0, length(d))
  }else if(method == "DS"){
    d = dweight0 * weight.scale
    if(entropy == 0){
      intercept = rep(0, length(d))      
    }else{
      intercept = rep(1, length(d))
    }
  }else if(method == "GEC0"){
    d = rep(1, nrow(Xs))
    intercept = g(dweight0 * weight.scale, entropy = entropy, del = del)
  }
  
  if(any(is.na(const)) & method == "GEC"){
    if(all(which(is.na(const)) == col_position)){
      nlmres= nlm(targetftn, p = What, d = d, Xs = Xs, init = init,
                  const = const, entropy = entropy, del = del, 
                  weight.scale = weight.scale, G.scale = G.scale,
                  intercept = intercept, K_alpha = K_alpha)
      # if(nlmres$code != 1 & nlmres$code != 2 & nlmres$code != 3) print(nlmres$estimate)
      if(nlmres$code != 1 & nlmres$code != 2 & nlmres$code != 3) stop(nlmres$code)
      W = nlmres$estimate
      if(nlmres$minimum >= .Machine$double.xmax){
        print(nlmres)
        w = NA
      }else{
        w = targetftn(W, d = d, Xs = Xs, init = init,
                      const = const, entropy = entropy, del = del,
                      weight.scale = weight.scale, G.scale = G.scale,
                      intercept = intercept, K_alpha = K_alpha, returnw = T)
      }
    }else{
      stop("NA appears in const outside g(d)")
    }
  }else{
    nleqslv_res = nleqslv::nleqslv(init, f, jac = h, d = d, Xs = Xs, 
                                   const = const, entropy = entropy, del = del,
                                   weight.scale = weight.scale, G.scale = G.scale,
                                   intercept = intercept, control = list(maxit = 1e5, allowSingular = T),
                                   xscalm = "auto")
    # control = control
    
    if(nleqslv_res$termcd != 1){
      # print("Xs"); print(head(Xs))
      # print("d"); print(d)
      # print(      return(f(init, d = d, Xs = Xs,
      #                      const = const, entropy = entropy, del = del,
      #                      intercept = intercept)))
      tmpval <- f(nleqslv_res$x, d = d, Xs = Xs, 
                  const = const, entropy = entropy, del = del, 
                  weight.scale = weight.scale, G.scale = G.scale,
                  intercept = intercept)
      
      if(any(is.nan(tmpval)) | (max(abs(tmpval)) > 1e-5)){
        print(nleqslv_res)
        w = NA      
      }else{
        w = NULL
      }
    }else{
      w = NULL  
    }
    if(is.null(w)){
      w = f(nleqslv_res$x, d = d, Xs = Xs, 
            const = const, entropy = entropy, del = del, 
            weight.scale = weight.scale, G.scale = G.scale,
            intercept = intercept, returnw = T)    
    }
  }
  return(list(w = w, method = method, entropy = entropy,
              Xs = Xs, weight.scale = weight.scale, G.scale = G.scale,
              dweight0 = dweight0, del = del, const = const, What = What,
              col_position = col_position, K_alpha = K_alpha))
}

