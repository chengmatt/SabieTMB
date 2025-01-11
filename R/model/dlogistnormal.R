dlogistnormal <- function(obs, pred, Sigma, give_log = TRUE) {
  fit_obs = obs[-length(obs)] # remove last bin since it's known
  mu = log(pred[-length(pred)]) # remove last bin since it's known
  mu = mu - log(pred[length(pred)]) # calculate log ratio
  logres = RTMB::dmvnorm(as.vector(fit_obs), as.vector(mu), Sigma = Sigma, log = give_log) # fit multivariate normal on log ratios
  return(logres)
}
