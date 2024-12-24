#' Title Get_Comp_Likelihoods
#' Gives negative log liklelihood values for composition data for a given year and a given fleet (fishery or survey)
#'
#' @param Exp Expected values (catch at age or survey index at age) indexed for a given year and fleet (structured as a matrix by age and sex)
#' @param Obs Observed values (catch at age or survey index at age) indexed for a given year and fleet (structured as a matrix by age and sex)
#' @param ISS Input sample size indexed for a given year and fleet (structured as a vector w/ sexes)
#' @param Wt_Mltnml Mutlinomial weight (if any) for a given fleet (structured as a vector w/ sexes)
#' @param Comp_Type Composition Parameterization Type (== 0, aggregated comps by sex, == 1, split comps by sex (no implicit sex ratio information), == 2, joint comps across sexes (implicit sex ratio information))
#' @param Likelihood_Type Composition Likelihood Type (== 0, Multinomial, == 1 Dirichlet Multinomial)
#' @param n_sexes Number of sexes modeled
#' @param age_or_len Age or length comps (== 0, Age, == 1, Length)
#' @param AgeingError Ageing Error matrix
#' @param ln_theta Log theta overdispersion for Dirichlet mutlinomial (scalar or vector depending on if 'Split' or 'Joint')
#'
#' @return Returns negative log likelihood for composition data (age and/or length)
#' @export
#'
#' @examples
Get_Comp_Likelihoods = function(Exp, Obs, ISS, Wt_Mltnml, ln_theta, Comp_Type, Likelihood_Type, n_sexes, age_or_len, AgeingError) {
  
  # Read in functions
  require(here)
  source(here("R", "model", "ddirmult.R")) # dirichlet multinomial
  
  comp_nLL = rep(0, n_sexes) # initialize nLL here
  
  # Aggregated comps by sex
  if(Comp_Type == 0) {
    # Expected Values
    tmp_Exp = Exp / matrix(data = colSums(Exp), nrow = nrow(Exp), ncol = ncol(Exp), byrow = TRUE) # Normalize
    if(age_or_len == 0) tmp_Exp = t(as.vector(rowSums(tmp_Exp) / n_sexes)) %*% AgeingError # Aggregate and Apply Ageing Error
    if(age_or_len == 1) tmp_Exp = t(as.vector(rowSums(tmp_Exp) / n_sexes)) # Aggregate
    tmp_Exp = tmp_Exp / sum(tmp_Exp) # renormalize
    tmp_Obs = Obs[,1] / sum(Obs[,1]) # Normalize observed values (indexing for sex 1 since comps are combined)
    
    # Multinomial likelihood
    if(Likelihood_Type == 0) { # Note that this indexes 1 because it's only a single sex 
      ESS = ISS[1] * Wt_Mltnml[1] # Effective sample size
      comp_nLL[1] = -1 * ESS * sum(((tmp_Obs + 0.001) * log(tmp_Exp + 0.001))) # ADMB multinomial likelihood
      comp_nLL[1] = comp_nLL[1] - -1 * ESS * sum(((tmp_Obs + 0.001) * log(tmp_Obs + 0.001))) # Multinomial offset (subtract offset from actual likelihood)
    } # end if multinomial likelihood
      if(Likelihood_Type == 1)  comp_nLL[1] = -1 * ddirmult(tmp_Obs, tmp_Exp, ISS, ln_theta, TRUE) # Dirichlet Multinomial likelihood
  } # end if aggregated comps across sex
  
  # 'Split' comps by sex (no implicit sex ratio information)
  if(Comp_Type == 1) {
    for(s in 1:n_sexes) {
      # Expected Values
      if(age_or_len == 0) tmp_Exp = (Exp[,s] / sum(Exp[,s])) %*% AgeingError # Normalize temporary variable (ages)
      if(age_or_len == 1) tmp_Exp = Exp[,s] / sum(Exp[,s]) # Normalize temporary variable (lengths)
      # Observed Values
      tmp_Obs = Obs[,s] / sum(Obs[,s]) # Normalize temporary variable
      # Multinomial likelihood
      if(Likelihood_Type == 0) { 
        ESS = ISS[s] * Wt_Mltnml[s] # Effective sample size
        comp_nLL[s] = -1 * ESS * sum(((tmp_Obs + 0.001) * log(tmp_Exp + 0.001))) # ADMB multinomial likelihood
        comp_nLL[s] = comp_nLL[s] - -1 * ESS * sum(((tmp_Obs + 0.001) * log(tmp_Obs + 0.001))) # Multinomial offset (subtract offset from actual likelihood)
      } # end if multinomial likelihood
      if(Likelihood_Type == 1)  comp_nLL[s] = -1 * ddirmult(tmp_Obs, tmp_Exp, ISS[s], ln_theta[s], TRUE) # Dirichlet Multinomial likelihood
    } # end s loop
  } # end if 'Split' comps by sex
  
  if(Comp_Type == 2) {
    # Expected values
    if(age_or_len == 0) { # if ages
      tmp_Exp = t(as.vector(Exp / sum(Exp))) %*% kronecker(diag(n_sexes), AgeingError) # apply ageing error
      tmp_Exp = as.vector(tmp_Exp / sum(tmp_Exp)) # renormalize to make sure sum to 1
      } # if ages
    if(age_or_len == 1) tmp_Exp = as.vector(Exp / sum(Exp)) # Normalize temporary variable (lengths)
    # Observed Values
    tmp_Obs = as.vector(Obs / sum(Obs)) # Normalize temporary variable
    # Multinomial likelihood 
    if(Likelihood_Type == 0) { # Indexing by 1 since it's a single vector 
      ESS = ISS[1] * Wt_Mltnml[1] # Effective sample size
      comp_nLL[1] = -1 * ESS * sum(((tmp_Obs + 0.001) * log(tmp_Exp + 0.001))) # ADMB multinomial likelihood
      comp_nLL[1] = comp_nLL[1] - -1 * ESS * sum(((tmp_Obs + 0.001) * log(tmp_Obs + 0.001))) # Multinomial offset (subtract offset from actual likelihood)
    } # end if multinomial likelihood
    if(Likelihood_Type == 1)  comp_nLL[1] = -1 * ddirmult(tmp_Obs, tmp_Exp, ISS, ln_theta, TRUE) # Dirichlet Multinomial likelihood
  } # end if 'Joint' comps by sex
 
  return(comp_nLL) # return negative log likelihood
} # end function