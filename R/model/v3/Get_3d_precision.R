#' Title Constructor algorithin for correlations within ages, years, and cohort
#'
#' @param n_ages Number of ages
#' @param n_yrs Number of years
#' @param pcorr_age correlations for age
#' @param pcorr_year correaltions for year
#' @param pcorr_cohort correlaitons for cohort
#' @param ln_var_value log space variance
#' @param Var_Type variance type == 0, marginal (stationary and slower run time), == 1 conditional (non-statationary, faster run time)
#'
#' @returns Sparse precision matrix dimensioned by n_ages * n_years, n_ages * n_years
#' @export
#'
#' @examples
Get_3d_precision <- function(n_ages,n_yrs, pcorr_age, pcorr_year, pcorr_cohort, ln_var_value, Var_Type){
    
    require(Matrix)

    index = expand.grid(seq_len(n_ages), seq_len(n_yrs)) # create index combinations to loop through
    i = j = x = NULL # initialize posiiton to fill in precision matrix
    var_value = exp(ln_var_value) # transform to normal space
    
    for(n in 1:nrow(index)){
      age = index[n,1] # get age index out of all index combinations
      year = index[n,2] # get year index out of all index combinations
      if(age > 1 ){
        i = c(i, n)
        j = c(j, which(index[,1] == (age-1) & index[,2] == year))
        x = c(x, pcorr_year) # year correaltion indexing
      }
      if(year > 1){
        i = c(i, n)
        j = c(j, which(index[,1]==age & index[,2]==(year-1)) )
        x = c(x, pcorr_age) # age correlation indexing
      }
      if( age>1 & year>1 ){
        i = c(i, n)
        j = c(j, which(index[,1]==(age-1) & index[,2] == (year-1)) )
        x = c(x, pcorr_cohort) # cohort correlation indexing
      }
    } # end n loop
    
    # Assemble SAR precision (should be sparse matrices)
    B = sparseMatrix(i = i, 
                     j = j, 
                     x = x, 
                     dims = rep(n_ages*n_yrs,2))
    
    # Identity matrix
    I = sparseMatrix(i = seq_len(n_ages*n_yrs),
                     j = seq_len(n_ages*n_yrs),
                     x = rep(1,n_ages*n_yrs))

    # Solve Omega recursively for stationary variance (accumulator function)
    if(Var_Type == 0) {
      L = solve(I-B) # solve to get accumulator function for stationary variance
      d = rep(NA, nrow(index))
      for(n in 1:nrow(index) ){
        if(n==1){
          d[n] = var_value
        }else{
          cumvar = sum(L[n,seq_len(n-1)] * d[seq_len(n-1)] * L[n,seq_len(n-1)])
          d[n] = (var_value-cumvar) / L[n,n]^2
        }
      } # end n loop
    } # end marginal variance (stationary variance)
    
    if(Var_Type == 1) d = var_value # conditional variance (non-stationary variance)
    
    Omega_inv = sparseMatrix(i = seq_len(n_ages*n_yrs), 
                             j = seq_len(n_ages*n_yrs),
                             x = 1/d, 
                             dims = rep(n_ages*n_yrs,2)) # get inverse of omega
    
    Q = (I-t(B)) %*% Omega_inv %*% (I-B) # solve for precision
    
    return(Q)
  }

