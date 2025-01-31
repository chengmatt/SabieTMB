#' Title Symmetric Beta Function
#'
#' @param p_val Parameter value
#' @param p_ub Upper Bound of Parameter
#' @param p_lb Lower Bound of Parameter
#' @param p_prsd SD of parameter, higher values have a stronger penalty on bounds, lower values have a more difuse penalty on bounds
#' @param log whether or not to return the log likelihood
#'
#' @returns
#' @export
#'
#' @examples
dbeta_symmetric <- function(p_val, p_ub, p_lb, p_prsd, log = TRUE) {
  # Calculate mu term
  mu <- -p_prsd * log((p_ub + p_lb)/2 - p_lb) - p_prsd * log(0.5)
  # Calculate the prior likelihood components
  term1 <- -mu 
  term2 <- -p_prsd * log(p_val - p_lb + 1e-4)
  term3 <- -p_prsd * log(1 - (p_val - p_lb - 1e-4)/(p_ub - p_lb))
  # Combine terms to get final prior likelihood
  nLL <- term1 + term2 + term3
  if(log == TRUE) return(nLL) else return(exp(nLL))
}
