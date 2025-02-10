#' Title Poisson that can take integer values
#'
#' @param x observations
#' @param pred predicted
#' @param give_log if giving log likelihood
#'
#' @returns
#' @export
#'
#' @examples
dpois_noint <- function(x, pred, give_log = TRUE) {
  logres <- -pred + x*log(pred) - lgamma(x+1)
  if(give_log == TRUE) return(logres) else return(exp(logres))
}