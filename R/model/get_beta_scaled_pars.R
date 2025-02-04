#' Title Get scaled alpha beta parameters for a scaled beta distribution (for steepness)
#'
#' @param mu mean in normal space
#' @param sd sd in normal space
#' @param p_lb lower bound
#' @param p_ub upper bound
#'
#' @returns
#' @export
#'
#' @examples https://stackoverflow.com/questions/75165770/beta-distribution-with-bounds-at-0-1-0-5
#' https://nmfs-ost.github.io/ss3-doc/SS330_User_Manual_release.html

get_beta_scaled_pars <- function(low, high,mu,sigma) {
  # convert mean and sd to alpha and beta
  scale = high - low
  mean = (mu - low) / scale
  var1 = (sigma / scale) ^ 2
  var = mean * (1 - mean) / var1 - 1
  stopifnot(var < 0) # beta dist is not defined if var < 0 
  # alpha and beta conversion
  a = mean * var
  b = (1 - mean) * var
  return(c(a,b,low,scale))
}



