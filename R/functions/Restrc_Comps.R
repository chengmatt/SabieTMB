#' Title Restructure composition values for use in Francis Reweighting
#'
#' @param Exp Expected values (catch at age or survey index at age) indexed for a given year and fleet (structured as a matrix by age and sex)
#' @param Obs Observed values (catch at age or survey index at age) indexed for a given year and fleet (structured as a matrix by age and sex)
#' @param ISS Input sample size indexed for a given year and fleet (structured as a vector w/ sexes)
#' @param Wt_Mltnml Mutlinomial weight (if any) for a given fleet (structured as a vector w/ sexes)
#' @param Comp_Type Composition Parameterization Type (== 0, aggregated comps by sex, == 1, split comps by sex (no implicit sex ratio information), == 2, joint comps across sexes (implicit sex ratio information))
#' @param n_sexes Number of sexes modeled
#' @param age_or_len Age or length comps (== 0, Age, == 1, Length)
#' @param AgeingError Ageing Error matrix
#'
#' @return
#' @export
#'
#' @examples
Restrc_Comps_Francis <- function(Exp, Obs, ISS, Wt_Mltnml, Comp_Type, n_sexes, age_or_len, AgeingError) {
  
  # Storage
  Exp_mat = matrix(0, dim(Exp)[1], dim(Exp)[2])
  Obs_mat = matrix(0, dim(Exp)[1], dim(Exp)[2])
  
  # Aggregated comps by sex
  if(Comp_Type == 0) {
    # Expected Values
    tmp_Exp = Exp / matrix(data = colSums(Exp), nrow = nrow(Exp), ncol = ncol(Exp), byrow = TRUE) # Normalize
    if(age_or_len == 0) tmp_Exp = t(as.vector(rowSums(tmp_Exp) / n_sexes)) %*% AgeingError # Aggregate and Apply Ageing Error
    if(age_or_len == 1) tmp_Exp = t(as.vector(rowSums(tmp_Exp) / n_sexes)) # Aggregate
    tmp_Exp = tmp_Exp / sum(tmp_Exp) # renormalize
    tmp_Obs = Obs[,1] / sum(Obs[,1]) # Normalize observed values (indexing for sex 1 since comps are combined)
    # Input into storage matrix
    Exp_mat[,1] = tmp_Exp
    Obs_mat[,1] = tmp_Obs
  } # end if aggregated comps across sex
  
  # 'Split' comps by sex (no implicit sex ratio information)
  if(Comp_Type == 1) {
    for(s in 1:n_sexes) {
      # Expected Values
      if(age_or_len == 0) tmp_Exp = (Exp[,s] / sum(Exp[,s])) %*% AgeingError # Normalize temporary variable (ages)
      if(age_or_len == 1) tmp_Exp = Exp[,s] / sum(Exp[,s]) # Normalize temporary variable (lengths)
      # Observed Values
      tmp_Obs = Obs[,s] / sum(Obs[,s]) # Normalize temporary variable
      # Input into storage matrix
      Exp_mat[,s] = tmp_Exp
      Obs_mat[,s] = tmp_Obs
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
    # Input into storage matrix
    Exp_mat[] = tmp_Exp
    Obs_mat[] = tmp_Obs
  } # end if 'Joint' comps by sex
  
  return(list(Exp = Exp_mat,
              Obs = Obs_mat))
  
} # end function