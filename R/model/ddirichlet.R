#' Title Dirichlet Likelihood
#'
#' @param x vector of values
#' @param alpha Expected values w/ concentration sum(alpha)
#' @param log Whether to give log or not
#'
#' @returns
#' @export
#'
#' @examples
ddirichlet <- function(x, alpha, log = TRUE) {
  logres = lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * log(x))
  if(log == TRUE) return(logres) else return(exp(logres))
} # end function
