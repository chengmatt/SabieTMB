#' Title Get process error log likelihoods (positive)
#'
#' @param PE_model Process error model == 1, iid on a functional form, == 2, random walk on a functional form, == 3 3d gmrf marginal on a functional form, == 4 3d gmrf conditional on a functional form
#' @param ln_devs array of log deviations, if PE_model %in% c(1,2), array should be dimensioned by n_regions, n_yrs, number of parameters, number of sexes
#' if PE_model %in% c(3,4), array should be dimensioned by n_regions, n_ages, n_years, n_sexes (note reversal on ages and years compared to the rest of code, because of how constructor algorihtim is set up)
#' @param PE_pars Process error parameters, if PE_model %in% c(1,2) then should be dimensioned by n_regions, number of parameters on functional form, and number of sexes
#' if PE_model %in% c(3,4) then should be dimensioned by n_regions, 4, and n_sexes, where 4 is the number of parameters (par1 = corr_ages, par2 = corr_year, par3 = corr_cohort, and par4 = log_var)
#' @param n_yrs Number of years
#' @param n_regions Number of regions
#' @param n_ages Number of ages
#' @param n_sexes Number of sexes
#'
#' @returns numeric value of log likelihood (in positive space) 
#' @export
#'
#' @examples
Get_PE_loglik <- function(PE_model, PE_pars, ln_devs, n_regions, n_yrs, n_ages, n_sexes) {
  
  require(here)
  source(here("R", "model", "v3", "Get_3d_precision.R")) # constructor algorithim for 3d precision matrix
  
  ll = 0 # initialize likelihood
  
  if(PE_model == 1) {
    
    
    # Flag need to make this accommodate pars ... (the dim par doesnt work because no longer list ... )
    # Other thing we need to accommodate is to turn off penalization if sharing parameters ...
    # Anything else? 
    
    
    n_pars = dim(ln_devs)[3] # get number of parameters in functional form
    for(i in 1:n_pars) {
      for(r in 1:n_regions) {
        for(y in 1:n_yrs) {
          for(s in 1:n_sexes) {
            ll = ll + dnorm(ln_devs[r,y,i,s,1], 0, exp(PE_pars[r,i,s,1]), TRUE)
          } # end s loop
        } # end y loop
      } # end r loop
    } # end i loop
  } # end iid process error
  
  if(PE_model == 2) {
    n_pars = dim(ln_devs)[3] # get number of parameters in functional form
    for(i in 1:n_pars) {
      for(r in 1:n_regions) {
        for(y in 1:n_yrs) {
          for(s in 1:n_sexes) {
            if(y == 1) {  # if y == 1, initialize with large sigma on dnorm 
              ll = ll + dnorm(ln_devs[r,1,i,s,1], 0, 50, TRUE)
            } else {
              ll = ll + dnorm(ln_devs[r,y,i,s,1], ln_devs[r,y-1,i,s,1], exp(PE_pars[r,i,s,1]), TRUE)
            } # end else
          } # end s loop
        } # end y loop
      } # end r loop
    } # end i loop  
  } # end rw process error
  
  if(PE_model %in% c(3,4)) {
    
    if(PE_model == 3) Var_Type = 0 # marginal variance
    if(PE_model == 4) Var_Type = 1 # conditional variance
    
    for(r in 1:n_regions) {
      for(s in 1:n_sexes) {
        # Construct precision matrix for 3d gmrf
        Q = Get_3d_precision(n_ages = n_ages, # number of ages
                             n_yrs = n_yrs,  # number of years
                             pcorr_age = PE_pars[r,1,s,1], # unconstrained partial correaltion by age
                             pcorr_year = PE_pars[r,2,s,1], # unconstrained partial correaltion by year
                             pcorr_cohort = PE_pars[r,3,s,1], # unconstrained partial correaltion by cohort
                             ln_var_value = PE_pars[r,4,s,1], # log variance
                             Var_Type = Var_Type) # variance type, == 0 (marginal), == 1 (conditional)
        
        # note that log deviations are dimensioned by region, age, year, sex because of precision constructor is set up
        ll = ll + RTMB::dgmrf(x = ln_devs[r,,,s,1], mu = 0, Q = Q, log = TRUE) 
        
      } # end s loop
    } # end r loop
    
  } # end 3dar1 process error
  
  return(ll)
} # return log likelihood


# PE_model = 3
# PE_pars = array(0.3, c(data$n_regions, 4, data$n_sexes)) # sigmas
# n_regions = data$n_regions
# n_yrs = length(data$years)
# n_ages = length(data$ages)
# n_sexes = data$n_sexes
# # ln_devs = array(0, c(data$n_regions, length(data$years), 2, data$n_sexes))
# ln_devs = array(0, c(data$n_regions, length(data$ages), length(data$years), data$n_sexes))
# 
# Get_PE_loglik(PE_model = 3, PE_pars = PE_pars, ln_devs = ln_devs, n_regions, n_yrs, n_ages, n_sexes)
# 
# rmvnorm(1, rep(0, n_ages * n_yrs), as.matrix(solve(Q)))
