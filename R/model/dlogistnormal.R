#' Title Logistic Normal Likelihood
#'
#' @param obs Vector of observed values in numbers (can be integers or non-integers - vector length of n_bins)
#' @param pred Vector of predicted values in proportions (vector length of n_bins)
#' @param Sigma Covariance matrix
#' @param give_log whether or not to use log space likelihood
#' @param perc percentage to multiply minimum values by so there aren't zeros
#'
#' @returns
#' @export
#'
#' @examples
dlogistnormal = function(obs, pred, Sigma, perc, zero_opt = "aitchison", give_log = TRUE) {

  # Dealing with zeros
  if(any(obs == 0)) {
    # normalize by adding 5% of the minimum observed value in a given year
    if(zero_opt == "min_perc") obs = (obs + min(obs[obs != 0]) * perc) / sum(obs + min(obs[obs != 0]) * perc)
    if(zero_opt == "aitchison") {
      obs = obs / sum(obs)
      zero_idx = which(obs == 0) # figure out number of zeros
      n_zero = length(zero_idx)
      delta = 0.0005 # rounding errors
      const = (delta * (n_zero + 1) * (length(obs) - n_zero)) / length(obs)^2
      obs[zero_idx] = obs[zero_idx] + const # add constant to zeros
      obs[-zero_idx] = obs[-zero_idx] - const # subtract constant to zeros
      obs = obs / sum(obs) # renomalize 
    } # if use aitchison method to add zeros
    # do logistic transformation on observed values
    tmp_Obs = log(obs[-length(obs)])
    tmp_Obs = tmp_Obs - log(obs[length(obs)])
    # do logistic transformation on expected values
    mu = log(pred[-length(pred)]) # remove last bin since it's known
    mu = mu - log(pred[length(pred)]) # calculate log ratio
  } else {
    obs = obs / sum(obs)
    # do logistic transformation on observed values
    tmp_Obs = log(obs[-length(obs)])
    tmp_Obs = tmp_Obs - log(obs[length(obs)])
    # do logistic transformation on expected values
    mu = log(pred[-length(pred)]) # remove last bin since it's known
    mu = mu - log(pred[length(pred)]) # calculate log ratio
  } # if we don't have any zeros

  res = RTMB::dmvnorm(as.vector(tmp_Obs), as.vector(mu), Sigma = Sigma, log = give_log) # fit multivariate normal on log ratios
  return(res)
}

