#' Title Get process error log likelihoods (positive)
#'

#' @param PE_model 
#' @param PE_pars 
#' @param ln_devs 
#' @param map_sel_devs 
#'
#' @returns numeric value of log likelihood (in positive space) 
#' @export
#'
#' @examples
Get_sel_PE_loglik <- function(PE_model, PE_pars, ln_devs, map_sel_devs) {
  
  require(here)
  source(here("R", "model", "v3", "Get_3d_precision.R")) # constructor algorithim for 3d precision matrix
  
  ll = 0 # initialize likelihood
  
  # find unique selectivity deviations to penalize (sort drops NAs)
  unique_sel_devs = sort(unique(as.vector(map_sel_devs))) 
  
  if(PE_model %in% c(1, 2)) {
    for(dev_idx in 1:length(unique_sel_devs)) {
      
      # figure out where unique sel devs first occur
      idx = which(map_sel_devs == unique_sel_devs[dev_idx], arr.ind = TRUE)[1,] 
      r = idx[1] # get unique region index
      y = idx[2] # get unique year index
      i = idx[3] # get unique age or parmeter index
      s = idx[4] # get unique sex index
      f = idx[5] # get unique fleet index      
      
      if(PE_model == 1) {
        ll = ll + dnorm(ln_devs[r,y,i,s,1], 0, exp(PE_pars[r,i,s,1]), TRUE) 
      } # iid process error
      
      if(PE_model == 2) { # random walk
        if(y == 1) { 
          ll = ll + dnorm(ln_devs[r,1,i,s,1], 0, 50, TRUE) # if y == 1, initialize with large sigma on dnorm
        } else {
          ll = ll + dnorm(ln_devs[r,y,i,s,1], ln_devs[r,y-1,i,s,1], exp(PE_pars[r,i,s,1]), TRUE)
        } # end else
      } # end random walk process error
      
    } # end dev_idx loop
  } # end iid or random walk process error
  
  if(PE_model %in% c(3,4)) {

    if(PE_model == 3) Var_Type = 0 # marginal variance
    if(PE_model == 4) Var_Type = 1 # conditional variance
    
    # Get indexing for constructor algorithim
    n_yrs = dim(map_sel_devs)[2]
    n_ages = dim(map_sel_devs)[3]
    
    # get first unique combinations
    unique_comb = which(map_sel_devs == unique_sel_devs[1], arr.ind = TRUE)[1,]
    # cbind to get all unique combinations here
    for(i in 2:length(unique_sel_devs)) unique_comb = cbind(unique_comb, which(map_sel_devs == unique_sel_devs[i], arr.ind = TRUE)[1,])

    # Next, get unique regionm, sex combinations
    unique_rs = expand.grid(unique(unique_comb[1,]), unique(unique_comb[4,]))
    
    for(idx in 1:nrow(unique_rs)) {
      r = unique_rs[idx,1] # get region index
      s = unique_rs[idx,2] # get sex index

      # Construct precision matrix for 3d gmrf
      if(PE_model %in% c(3,4)) {
        Q = Get_3d_precision(n_ages = n_ages, # number of ages
                             n_yrs = n_yrs,  # number of years
                             pcorr_age = PE_pars[r,1,s,1], # unconstrained partial correaltion by age
                             pcorr_year = PE_pars[r,2,s,1], # unconstrained partial correaltion by year
                             pcorr_cohort = PE_pars[r,3,s,1], # unconstrained partial correaltion by cohort
                             ln_var_value = PE_pars[r,4,s,1], # log variance
                             Var_Type = Var_Type) # variance type, == 0 (marginal), == 1 (conditional)
        
        # apply gmrf likelihood
        ll = ll + RTMB::dgmrf(x = t(ln_devs[r,,,s,1]), mu = 0, Q = Q, log = TRUE)
      } # end if
      
    } # end idx loop
  } # end 3dar1 process error

  return(ll)
} # return log likelihood
