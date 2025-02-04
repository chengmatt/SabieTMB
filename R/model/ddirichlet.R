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
  first_term <- sum(lgamma(alpha)) - lgamma(sum(alpha))
  second_term <- sum((alpha - 1) * log(x))
  logres = first_term - second_term
  if(log == TRUE) return(logres) else return(exp(logres))
} # end function
