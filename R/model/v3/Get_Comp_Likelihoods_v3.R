#' Title Get_Comp_Likelihoods
#' Gives negative log liklelihood values for composition data for a given year and a given fleet (fishery or survey)
#'
#' @param Exp Expected values (catch at age or survey index at age) indexed for a given year and fleet (structured as a matrix by age and sex)
#' @param Obs Observed values (catch at age or survey index at age) indexed for a given year and fleet (structured as a matrix by age and sex)
#' @param ISS Input sample size indexed for a given year and fleet (structured as a vector w/ sexes)
#' @param Wt_Mltnml Mutlinomial weight (if any) for a given fleet (structured as a vector w/ sexes)
#' @param Comp_Type Composition Parameterization Type (== 0, aggregated comps by sex, == 1, split comps by sex and region (no implicit sex and region ratio information), 
#' == 2, joint comps across sexes but split by region (implicit sex ratio information, but not region information)), == 3, joint comps across sexes and regions (implicit sex and region ratio information)
#' @param Likelihood_Type Composition Likelihood Type (== 0, Multinomial, == 1 Dirichlet Multinomial)
#' @param n_sexes Number of sexes modeled
#' @param age_or_len Age or length comps (== 0, Age, == 1, Length)
#' @param AgeingError Ageing Error matrix
#' @param ln_theta Log theta overdispersion for Dirichlet mutlinomial (scalar or vector depending on if 'Split' or 'Joint')
#' @param n_regions number of regions modeled
#' @param n_bins number of bins modeled
#' @param use Vector of 0s and 1s corresponding to regions (==0, don't have obs and dont' use, ==1, have obs and use)
#'
#' @return Returns negative log likelihood for composition data (age and/or length)
#' @export
#'
#' @examples
Get_Comp_Likelihoods = function(Exp,
                                Obs, 
                                ISS, 
                                Wt_Mltnml,
                                ln_theta, 
                                Comp_Type, 
                                Likelihood_Type, 
                                n_regions, 
                                n_bins,
                                n_sexes, 
                                age_or_len, 
                                AgeingError, 
                                use) {
  
  # Read in functions
  require(here)
  source(here("R", "model", "ddirmult.R")) # dirichlet multinomial
  source(here("R", "model", "dlogistnormal.R")) # logistic normal
  
  comp_nLL = array(0, dim = c(n_regions, n_sexes)) # initialize nLL here
  const = 1e-10 # small constant
  # Filter expectation and observations to regions that have observations
  n_regions_obs_use = sum(use == 1) # get number of regions that have observations

  # Obs = data$ObsFishAgeComps[,1,,1,1]
  # Exp = sabie_rtmb_model$rep$CAA[,1,,1,1]
  
  # Making sure things are correctly formatted (and regions are not dropped)
  Obs = array(Obs, dim = c(n_regions, n_bins, n_sexes)) 
  Exp = array(Exp, dim = c(n_regions, n_bins, n_sexes)) 
  ISS = array(ISS, dim = c(n_regions, n_sexes)) 
  Wt_Mltnml = array(Wt_Mltnml, dim = c(n_regions, n_sexes)) 
  ln_theta = array(ln_theta, dim = c(n_regions, n_sexes)) 
  
  # filter regions that have obs
  Obs = Obs[which(use == 1),,,drop = FALSE] 
  Exp = Exp[which(use == 1),,,drop = FALSE] 

  # Aggregated comps by sex and region
  if(Comp_Type == 0) {
    # Expected Values
    tmp_Exp = Exp / array(data = rowSums(Exp, dims = 1), dim = dim(Exp)) # normalize by sex and region
    tmp_Exp = matrix(colSums(tmp_Exp) / (n_sexes * n_regions), nrow = 1) # take average proportions and transpose
    if(age_or_len == 0) tmp_Exp = tmp_Exp %*% AgeingError # apply ageing error
    tmp_Exp = as.vector((tmp_Exp + const) / sum(tmp_Exp + const)) # renormalize
    
    # Multinomial likelihood
    if(Likelihood_Type == 0) { # Note that this indexes 1 because it's only a single sex and single region
      tmp_Obs = apply(Obs + const, 2:3, sum) / sum(Obs + const) # Normalize observed values 
      ESS = ISS[1,1] * Wt_Mltnml[1,1] # Effective sample size
      comp_nLL[1,1] = -1 * ESS * sum(((tmp_Obs) * log(tmp_Exp))) # ADMB multinomial likelihood
      comp_nLL[1,1] = comp_nLL[1,1] - -1 * ESS * sum(((tmp_Obs) * log(tmp_Obs))) # Multinomial offset (subtract offset from actual likelihood)
    } # end if multinomial likelihood
    
    if(Likelihood_Type == 1) {
      tmp_Obs = apply(Obs + const, 2:3, sum) / sum(Obs + const) # Normalize observed values 
      comp_nLL[1,1] = -1 * ddirmult(tmp_Obs, tmp_Exp, ISS, ln_theta[1,1], TRUE) # Dirichlet Multinomial likelihood
    } # end if dirichlet multinomial 
    
    if(Likelihood_Type == 2) {
      tmp_Obs = apply(Obs, 2:3, sum) / sum(Obs) # Normalize observed values 
      Sigma = diag(length(tmp_Obs)-1) * (exp(ln_theta[1,1])^2 / ISS[1,1])
      comp_nLL[1,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma = Sigma, TRUE) # Logistic Normal likelihood
    } # end if logistic normal

  } # end if aggregated comps across sexes and regions
  
  # 'Split' comps by sex and region (no implicit sex ratio information)
  if(Comp_Type == 1) {
    for(s in 1:n_sexes) {
      for(r in 1:n_regions_obs_use) {
        # Expected Values
        if(age_or_len == 0) tmp_Exp = ((Exp[r,,s] + const) / sum(Exp[r,,s] + const)) %*% AgeingError # Normalize temporary variable (ages)
        if(age_or_len == 1) tmp_Exp = (Exp[r,,s] + const) / sum(Exp[r,,s] + const) # Normalize temporary variable (lengths)

        # Multinomial likelihood
        if(Likelihood_Type == 0) { 
          tmp_Obs = (Obs[r,,s] + const) / sum(Obs[r,,s] + const) # Normalize observed temporary variable
          ESS = ISS[r,s] * Wt_Mltnml[r,s] # Effective sample size
          comp_nLL[r,s] = -1 * ESS * sum(((tmp_Obs) * log(tmp_Exp))) # ADMB multinomial likelihood
          comp_nLL[r,s] = comp_nLL[r,s] - -1 * ESS * sum(((tmp_Obs) * log(tmp_Obs))) # Multinomial offset (subtract offset from actual likelihood)
        } # end if multinomial likelihood
        
        if(Likelihood_Type == 1) {
          tmp_Obs = (Obs[r,,s] + const) / sum(Obs[r,,s] + const) # Normalize observed temporary variable
          comp_nLL[r,s] = -1 * ddirmult(tmp_Obs, tmp_Exp, ISS[r,s], ln_theta[r,s], TRUE) # Dirichlet Multinomial likelihood
        } # end if dirichlet multinomial 
        
        if(Likelihood_Type == 2) {
          tmp_Obs = Obs[r,,s] / sum(Obs[r,,s]) # extract variable and normalize
          Sigma = diag(length(tmp_Obs)-1) * (exp(ln_theta[r,s])^2 / ISS[r,s]) # construct sigma
          comp_nLL[r,s] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma = Sigma, TRUE) # Logistic Normal likelihood 
        } # end if logistic normal 
        
      } # end r loop
    } # end s loop
  } # end if 'Split' comps by sex and region
  
  # Joint by sex, Split by region
  if(Comp_Type == 2) {
    for(r in 1:n_regions_obs_use) {
      # Expected values
      if(age_or_len == 0) { # if ages
        tmp_Exp = t(as.vector((Exp[r,,] + const)/ sum(Exp[r,,] + const))) %*% kronecker(diag(n_sexes), AgeingError) # apply ageing error
        tmp_Exp = as.vector((tmp_Exp + const) / sum(tmp_Exp + const)) # renormalize to make sure sum to 1
      } # if ages
      if(age_or_len == 1) tmp_Exp = as.vector((Exp[r,,] + const) / sum((Exp[r,,] + const))) # Normalize temporary variable (lengths)

      # Multinomial likelihood 
      if(Likelihood_Type == 0) { # Indexing by r for a given region since it's 'Split' by region and 1 for sex since it's 'Joint' for sex
        tmp_Obs = as.vector((Obs[r,,] + const) / sum(Obs[r,,] + const)) # Normalize observed temporary variable
        ESS = ISS[r,1] * Wt_Mltnml[r,1] # Effective sample size
        comp_nLL[r,1] = -1 * ESS * sum(((tmp_Obs) * log(tmp_Exp))) # ADMB multinomial likelihood
        comp_nLL[r,1] = comp_nLL[r,1] - -1 * ESS * sum(((tmp_Obs) * log(tmp_Obs))) # Multinomial offset (subtract offset from actual likelihood)
      } # end if multinomial likelihood
      
      if(Likelihood_Type == 1) {
        tmp_Obs = as.vector((Obs[r,,] + const) / sum(Obs[r,,] + const)) # Normalize observed temporary variable
        comp_nLL[r,1] = -1 * ddirmult(tmp_Obs, tmp_Exp, ISS[r,1], ln_theta[r,1], TRUE) # Dirichlet Multinomial likelihood
      } # end if dirichlet multinomial
      
      if(Likelihood_Type == 2) {
        tmp_Obs = Obs[r,,] / sum(Obs[r,,]) # extract temporary observed variable and normalize
        Sigma = diag(length(tmp_Obs)-1) * (exp(ln_theta[r,1])^2 / ISS[r,1])
        comp_nLL[r,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma = Sigma, perc = 0.05, zero_opt = "aitchison", TRUE) # Logistic Normal likelihood
      } # end if dirichlet multinomial
      
    } # end r loop
  } # end if 'Joint' comps by sex, but 'Split' by region
  
  # Joint by sex and region
  if(Comp_Type == 3) {
    tmp_Exp = aperm(Exp, perm = c(2,3,1)) # Reformat expected values so it's ordered by ages, sexes, and then regions
    tmp_Obs = aperm(Obs, perm = c(2,3,1)) # Reformat observed values so it's ordered by ages, sexes, and then regions
    
    # Expected values
    if(age_or_len == 0) { # if ages
      tmp_Exp = t(as.vector((tmp_Exp + const) / sum(tmp_Exp + const))) %*% kronecker(diag(n_regions_obs_use * n_sexes), AgeingError) # apply ageing error
      tmp_Exp = as.vector((tmp_Exp + const)/ sum(tmp_Exp + const)) # renormalize to make sure sum to 1
      tmp_Exp = as.vector((tmp_Exp + const)/ sum(tmp_Exp + const)) # renormalize to make sure sum to 1
    } # if ages
    
    if(age_or_len == 1) tmp_Exp = as.vector((Exp + const) / sum(Exp + const)) # Normalize temporary variable (lengths)
    
    # Multinomial likelihood 
    if(Likelihood_Type == 0) { # Indexing by 1,1 because Joint by sex and regions
      tmp_Obs = as.vector((tmp_Obs + const) / sum(tmp_Obs + const)) # Normalize observed temporary variable   
      ESS = ISS[1,1] * Wt_Mltnml[1,1] # Effective sample size
      comp_nLL[1,1] = -1 * ESS * sum(((tmp_Obs) * log(tmp_Exp))) # ADMB multinomial likelihood
      comp_nLL[1,1] = comp_nLL[1,1] - -1 * ESS * sum(((tmp_Obs) * log(tmp_Obs))) # Multinomial offset (subtract offset from actual likelihood)
    } # end if multinomial likelihood
    
    if(Likelihood_Type == 1)  {
      tmp_Obs = as.vector((tmp_Obs + const) / sum(tmp_Obs + const)) # Normalize temporary variable   
      comp_nLL[1,1] = -1 * ddirmult(tmp_Obs, tmp_Exp, ISS[1,1], ln_theta[1,1], TRUE) # Dirichlet Multinomial likelihood
    } # end if Dirichlet multinomial
    
    if(Likelihood_Type == 2) { 
      Sigma = diag(length(tmp_Obs)-1) * (exp(ln_theta[1,1])^2 / ISS[1,1])
      comp_nLL[1,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma = Sigma, perc = 0.05, zero_opt = "aitchison", TRUE) # Logistic Normal likelihood
    } # Logistic normal likelihood
    
  } # end if comps are joint by sex and regions
  
  return(comp_nLL) # return negative log likelihood
} # end function

