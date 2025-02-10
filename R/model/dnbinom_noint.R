#' Title Negative binomial that can take non-integer values
#'
#' @param x observations
#' @param mu expectation
#' @param size dispersion
#' @param give_log whether to give log
#'
#' @returns
#' @export
#'
#' @examples
dnbinom_noint <- function(x, mu, size, give_log = TRUE) {
  prob <- size/(size + mu)
  logres <- lgamma(x + size) - lgamma(size) - lgamma(x + 1) +size * log(prob) + x * log(1 - prob)
  if(give_log) return(logres) else return(exp(logres))
}
