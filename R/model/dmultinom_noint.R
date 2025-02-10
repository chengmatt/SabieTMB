dmultinom_noint <- function(x, size, prob, give_log = T) {
  logres <- lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))
  if(give_log) return(logres) else return(exp(logres))
}