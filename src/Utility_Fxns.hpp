// Purpose: An assortment of utility functions and likelihoods
// Creator: Matthew LH. Cheng (UAF-CFOS)
// Date: 3/18/23

template <class Type> // Function to call different selectivity parameterizations
// @param age = indexed age within age loop
// @param sel_model = integer of selectivity model 
// == 0, logistic
// @param ln_selpars vector of log selectivity parameters.
Type Get_Selex(Type age, 
               int sel_model, 
               vector<Type> ln_selpars) {
  
  // Create container to return predicted selectivity value
  Type selex = 0;
  
  if(sel_model == 0) { // logistic selectivity
    // Extract out and exponentiate the parameters here
    Type a50 = exp(ln_selpars(0)); // a50
    Type k = exp(ln_selpars(1)); // slope
    selex = Type(1.0) / (Type(1) + exp(-k * (age - a50)));
  }
  
  if(sel_model == 1) { // gamma dome-shaped selectivity 
    // Extract out and exponentiate the parameters here
    Type amax = exp(ln_selpars(0)); // age at max selex
    Type delta = exp(ln_selpars(1)); // slope parameter

    // Now, calculate/derive power parameter + selex values
    Type p = 0.5 * (sqrt( pow(amax, 2) + (4 * pow(delta, 2))) - amax);
    selex = pow( (age / amax), (amax/p)  ) * exp( (amax - age) / p ); 
    
  }
  
  if(sel_model == 2) { // power function selectivity
    // Extract out and exponentiate the parameters here
    Type power = exp(ln_selpars(0)); // power parameter
    selex = 1 / pow(age, power);

  }
  
  return selex;
} // end function

template <class Type> // sum to zero constraint
vector<Type> sum_to_zero(vector<Type> vec, Type dev_NA) {
  Type sum_of_vec = 0; // initialize
  Type adjustment = 0; // initialize
  vector<Type> vec1(vec.size());
  sum_of_vec = vec.sum(); // get sum of vector
  adjustment = -sum_of_vec / (vec.size() - dev_NA); // get adjustement factor
  for(int i = 0; i < vec1.size(); i++) vec1(i) = vec(i) + adjustment;
  return vec1;
} // end function

template <class Type> // Square a variable
Type square(Type x){return x*x;}

template<class Type> // Function for detecting NAs
bool isNA(Type x){return R_IsNA(asDouble(x));}
