#' Title Dirichlet Mutlinomial Likelihood
#' From https://github.com/James-Thorson/CCSRA/blob/main/inst/executables/CCSRA_v9.cpp
#' @param obs Vector of observed values in proportions
#' @param pred Vector or predicted values in proportions
#' @param Ntotal Input sample size scalar
#' @param ln_theta Weighting parameter in log space
#' @param give_log Whether or not likelihood is in log space
#'
#' @returns
#' @export
#'
#' @examples
ddirmult = function(obs, pred, Ntotal, ln_theta, give_log = TRUE) {
  # Set up function variables
  n_c = length(obs) # number of categories
  p_exp = pred # expected values container
  p_obs = obs # observed values container
  dirichlet_Parm = exp(ln_theta) * Ntotal # Dirichlet alpha parameters
  
  # set up pdf
  logres = lgamma(Ntotal + 1)
  for(c in 1:n_c) logres = logres - lgamma(Ntotal*p_obs[c]+1) # integration constant
  logres = logres + lgamma(dirichlet_Parm) - lgamma(Ntotal+dirichlet_Parm) # 2nd term in formula
  
  # Summation in 3rd term in formula
  for(c in 1:n_c) {
    logres = logres + lgamma(Ntotal*p_obs[c] + dirichlet_Parm*p_exp[c])
    logres = logres - lgamma(dirichlet_Parm * p_exp[c])
  } # end c
  
  if(give_log == TRUE) return(logres)
  else return(exp(logres))
} # end function

