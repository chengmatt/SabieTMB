#' Title Get Selectivity
#'
#' @param Selex_Model integer for selectivity model
#' @param ln_Pars log selectivity parameters
#' @param Age vector of Ages
#'
#' @return
#' @export
#'
#' @examples
Get_Selex = function(Selex_Model, ln_Pars, Age) {
  selex = rep(0, length(Age)) # Temporary container vector
  
  if(Selex_Model == 0) { # logistic selectivity (a50 and slope)
    # Extract out and exponentiate the parameters here
    a50 = exp(ln_Pars[1]); # a50
    k = exp(ln_Pars[2]); # slope
    selex = 1 / (1 + exp(-k * (Age - a50))) 
  }
  
  if(Selex_Model == 1) { # gamma dome-shaped selectivity 
    # Extract out and exponentiate the parameters here
    amax = exp(ln_Pars[1]) # age at max selex
    delta = exp(ln_Pars[2]) # slope parameter
    # Now, calculate/derive power parameter + selex values
    p = 0.5 * (sqrt( amax^2 + (4 * delta^2)) - amax)
    selex = (Age / amax)^(amax/p) * exp( (amax - Age) / p ) 
  }
  
  if(Selex_Model == 2) { # power function selectivity
    # Extract out and exponentiate the parameters here
    power = exp(ln_Pars[1]); # power parameter
    selex = 1 / Age^power
  }
  
  
  if(Selex_Model == 3) { # logistic selectivity (a50 and a95)
    # Extract out and exponentiate the parameters here
    a50 = exp(ln_Pars[1]); # a50
    a95 = exp(ln_Pars[2]); # slope
    selex = 1 / (1+19^((a50-Age)/a95))
  }
  
  return(selex)
} # end function