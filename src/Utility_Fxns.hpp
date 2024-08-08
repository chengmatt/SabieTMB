// Purpose: An assortment of utility functions and likelihoods
// Creator: Matthew LH. Cheng (UAF-CFOS)
// Date: 3/18/23

// Selectivity Function --------------------------------------------------------
// Function to call different selectivity parameterizations
// @param age = indexed age within age loop
// @param sel_model = integer of selectivity model 
// == 0, logistic
// @param ln_selpars vector of log selectivity parameters.
template <class Type> 
Type Get_Selex(Type age, 
               int sel_model, 
               vector<Type> ln_selpars) {
  
  // Container to return predicted selectivity value
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


// Other Ancillary Functions ---------------------------------------------------
// function to sum to zero constraint
// @param vec Vector to sum to zero,
// @param dev_NA size of vector with NAs (i.e., if zero NAs, dev_NA should be specified as 0)
template <class Type> 
vector<Type> sum_to_zero(vector<Type> vec, Type dev_NA) {
  Type sum_of_vec = 0; // initialize
  Type adjustment = 0; // initialize
  vector<Type> vec1(vec.size());
  sum_of_vec = vec.sum(); // get sum of vector
  adjustment = -sum_of_vec / (vec.size() - dev_NA); // get adjustement factor
  for(int i = 0; i < vec1.size(); i++) vec1(i) = vec(i) + adjustment;
  return vec1;
} // end function

// Square a variable
// @ param x variable to square
template <class Type> 
Type square(Type x){
  return x*x;
  }

// Function for detecting NAs
// @ param x variable to detect NAs for
template<class Type> 
bool isNA(Type x){
  return R_IsNA(asDouble(x));
  }

// Function for taking the arithmetic mean
// @ param vector to take the mean of
template <class Type> 
Type arith_mean(vector<Type> x) {
  int n = x.size();
  Type mean = x.sum()/n;
  return mean;
}