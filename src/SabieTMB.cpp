// Purpose: Generalized Age-Structured Model tailored for Alaska Sablefish
// Generalizes ages, sexes, fishery fleets, and survey fleets
// Generalizes time block structure

#include<TMB.hpp>
#include"Utility_Fxns.hpp"

template<class Type> 
Type objective_function<Type>::operator() () 
{
  
  using namespace density; // Define namespace to use multivariate distributions
  using namespace Eigen; // Define namespace for Eigen functions (i.e., sparse matrix)
  
  // DATA SECTION --------------------------------------------------------------
  
  // Model Dimensions 
  DATA_VECTOR(ages); // vector of ages
  DATA_VECTOR(lens); // vector of lengths
  DATA_VECTOR(yrs); // vector of years
  DATA_INTEGER(n_sexes); // integer of sexes
  DATA_INTEGER(n_fish_fleets); // integer of number of fishery fleets
  DATA_INTEGER(n_srv_fleets); // integer of number of survey fleets

  // Biological Processes
  DATA_INTEGER(do_rec_bias_ramp); // switch for if we want to do bias ramp from Methot and Taylor 2011 (==0 no ramp), (==1 ramp)
  DATA_VECTOR(bias_year); // breakpoints for when to change bias ramp (slot 0 = period w/o bias correction, slot 1 = ascending, slot 2 = descending, slot 3 = no bias correction in last year)
  DATA_INTEGER(sigmaR_switch); // switch for sigmaR (year we want to switch sigmaR - from low rec variability to high rec variability)
  // DATA_VECTOR(mean_rec_period); // Vector for periods to use for mean recruitment
  DATA_ARRAY(WAA); // array of weight-at-age (n_years, n_ages, n_sexes)
  DATA_ARRAY(MatAA); // array of maturity at age (n_years, n_ages, n_sexes)
  DATA_VECTOR(sexratio); // vector of recruitment sex ratio (n_sexes)
  DATA_MATRIX(AgeingError); // matrix of ageing error (n_ages, n_ages) 
  DATA_ARRAY(SizeAgeTrans); // array for size age transition (n_years, n_lens, n_ages, n_sexes)
  
  // Observations (Data)
  // Catches
  DATA_SCALAR(init_F_prop); // initial F proportion for initial age structure (fished age structure)
  DATA_VECTOR(Catch_Constant); // vector of catch constants to add (0.01 for fixed gear, 0.8 for trawl)
  DATA_MATRIX(ObsCatch); // matrix of observed catches (n_years, n_fish_fleets)
  DATA_IMATRIX(UseCatch); // indicator for using catch (n_years, n_fish_fleets), == 0 not using, == 1 using
  
  // Indices
  DATA_MATRIX(ObsFishIdx); // matrix of observed fishery indices (n_years, n_fish_fleets)
  DATA_MATRIX(ObsFishIdx_SE); // matrix of observed fishery indices standard errors (n_years, n_fish_fleets)
  DATA_IMATRIX(UseFishIdx); // indicator for using fishery indices (n_years, n_fish_fleets), == 0 not using, == 1 using
  DATA_MATRIX(ObsSrvIdx); // matrix of observed survey indices (n_years, n_srv_fleets)
  DATA_MATRIX(ObsSrvIdx_SE); // matrix of observed survey indices standard errors (n_years, n_srv_fleets)
  DATA_IMATRIX(UseSrvIdx); // indicator for using survey indices (n_years, n_srv_fleets), == 0 not using, == 1 using
  
  // Fishery Age Compositions
  DATA_ARRAY(ObsFishAgeComps); // array of observed fishery age compositions (n_years, n_ages, n_sexes, n_fish_fleets)
  DATA_ARRAY(ISS_FishAgeComps); // array of input sample sizes for fishery age compositions (n_years, n_sexes, n_fish_fleets)
  DATA_MATRIX(Wt_FishAgeComps); // matrix of weights for for fishery age compositions - can be derived subjectively or via reweighting approaches (n_sexes, n_fish_fleets)
  DATA_MATRIX(UseFishAgeComps); // indicator for using fishery age comps (n_years, n_fish_fleets), == 0 not using, == 1 using
  DATA_MATRIX(AggFishAgeComps); // indicator for aggregating fishery age comps (n_years, n_fish_fleets), == 0 aggregating, == 1 sex-specific
  
  // Fishery Length Compositions
  DATA_ARRAY(ObsFishLenComps); // array of observed fishery length compositions (n_years, n_lens, n_sexes, n_fish_fleets)
  DATA_ARRAY(ISS_FishLenComps); // array of input sample sizes for fishery length compositions (n_years, n_sexes, n_fish_fleets)
  DATA_MATRIX(Wt_FishLenComps); // matrix of weights for for fishery length compositions - can be derived subjectively or via reweighting approaches (n_sexes, n_fish_fleets)
  DATA_MATRIX(UseFishLenComps); // indicator for using fishery length comps (n_years, n_fish_fleets), == 0 not using, == 1 using

  // Survey Age Compositions
  DATA_ARRAY(ObsSrvAgeComps); // array of observed survey age compositions (n_years, n_ages, n_sexes, n_srv_fleets)
  DATA_ARRAY(ISS_SrvAgeComps); // array of input sample sizes for survey age compositions (n_years, n_sexes, n_srv_fleets)
  DATA_MATRIX(Wt_SrvAgeComps); // matrix of weights for for survey age compositions - can be derived subjectively or via reweighting approaches (n_sexes, n_srv_fleets)
  DATA_MATRIX(UseSrvAgeComps); // indicator for using survey age comps (n_years, n_srv_fleets), == 0 not using, == 1 using
  DATA_MATRIX(AggSrvAgeComps); // indicator for aggregating survey age comps (n_years, n_srv_fleets), == 0 aggregating, == 1 sex-specific
  
  // Survey Length Compositions
  DATA_ARRAY(ObsSrvLenComps); // array of observed survey length compositions (n_years, n_lens, n_sexes, n_srv_fleets)
  DATA_ARRAY(ISS_SrvLenComps); // array of input sample sizes for survey length compositions (n_years, n_sexes, n_srv_fleets)
  DATA_MATRIX(Wt_SrvLenComps); // matrix of weights for for survey length compositions - can be derived subjectively or via reweighting approaches (n_sexes, n_srv_fleets)
  DATA_MATRIX(UseSrvLenComps); // indicator for using survey length comps (n_years, n_srv_fleets), == 0 not using, == 1 using
  
  // Model Specifications
  // Specifications for Fishery Processes
  DATA_INTEGER(share_sel); // whether to share sex specific selectivity for index calculations (sablefish specific)
  DATA_IVECTOR(normalize_fish_sel); // normalize fishery selectivity to max at 1 == 0, no normalization, == 1 normalize
  DATA_IMATRIX(fish_sel_blocks); // matrix of specifying fishery selectivity blocks (n_years, n_fish_fleets)
  DATA_IMATRIX(fish_q_blocks); // matrix of specifying fishery catchability blocks (n_years, n_fish_fleets)
  DATA_IVECTOR(fish_idx_type); // vector of specifying fishery index type (n_fish_fleets), == 0, abundance, == biomass
  DATA_IMATRIX(fish_sel_model); // matrix of specifying fishery selectivity models (n_years, n_fish_fleets)
  
  // Specifications for Survey Processes
  DATA_IVECTOR(normalize_srv_sel); // normalize survey selectivity to max at 1 == 0, no normalization, == 1 normalize
  DATA_IMATRIX(srv_sel_blocks); // matrix of specifying survey selectivity blocks (n_years, n_srv_fleets)
  DATA_IMATRIX(srv_q_blocks); // matrix of specifying survey catchability blocks (n_years, n_srv_fleets)
  DATA_IVECTOR(srv_idx_type); // vector of specifying survey index type (n_srv_fleets), == 0, abundance, == biomass
  DATA_IMATRIX(srv_sel_model); // matrix of specifying survey selectivity models (n_years, n_srv_fleets)
  
  // Likelihoods and priors
  DATA_SCALAR(Wt_Catch); // Scalar value for catch weights
  DATA_SCALAR(Wt_FishIdx); // Scalar value for fishery index weights
  DATA_SCALAR(Wt_SrvIdx); // Scalar value for survey index weights
  DATA_SCALAR(Wt_Rec); // Scalar value for recruitment weights
  DATA_SCALAR(Wt_F); // Scalar value for fishiery penalty weights
  DATA_INTEGER(Use_M_prior); // indicator for using natural mortality prior (== 0, don't use, == 1, use)
  DATA_VECTOR(M_prior); // prior for natural mortality dim 1 = mean, dim 2 = sd

  // PARAMETER SECTION ---------------------------------------------------------
  PARAMETER(dummy); // dummy parameter for testing
  
  // Biological Processes
  PARAMETER(ln_R0); // virgin recruitment or mean recruitment in log space
  PARAMETER(ln_sigmaR_early); // early sigma R for initial devs
  PARAMETER(ln_sigmaR_late); // late sigma R for initial devs
  PARAMETER_VECTOR(ln_InitDevs); // initial deviations for first year age structure (n_ages - 2)
  PARAMETER_VECTOR(ln_RecDevs); // recruitment deviations (n_years)
  PARAMETER(ln_M); // base natural mortality
  PARAMETER(M_offset); // natural mortality offset (additive, e.g., additional sex)
  
  // Fishery Processes
  PARAMETER_VECTOR(ln_F_mean); // vector of mean fishing mortality (n_fish_fleets)
  PARAMETER_MATRIX(ln_F_devs); // matrix of fishing mortality deviations (n_fish_fleets, n_fish_fleets)
  PARAMETER_ARRAY(ln_fish_fixed_sel_pars); // array of fishery selectivity pars (n_pars, n_blocks, n_sexes, n_fish_fleets)
  PARAMETER_MATRIX(ln_fish_q); // matrix of fishery catchability parameters (n_blocks, n_fish_fleets)
  
  // Survey Processes
  PARAMETER_ARRAY(ln_srv_fixed_sel_pars); // array of survey selectivity pars (n_pars, n_blocks, n_sexes, n_srv_fleets)
  PARAMETER_MATRIX(ln_srv_q); // matrix of survey catchability parameters (n_blocks, n_srv_fleets)
  
  // CONTAINER SECTION ---------------------------------------------------------
  int n_ages = ages.size(); // number of ages
  int n_yrs = yrs.size();  // number of years
  int n_lens = lens.size(); // number of lengths
  
  // Population Dynamics Processes
  array<Type> NAA(n_yrs + 1, n_ages, n_sexes); // array of numbers-at-age
  array<Type> ZAA(n_yrs, n_ages, n_sexes); // array of total mortality at age
  array<Type> SAA(n_yrs, n_ages, n_sexes); // array of survival at age annual
  array<Type> SAA_mid(n_yrs, n_ages, n_sexes); // array of survival at age midpoint
  vector<Type> SSB(n_yrs); // spawning stock biomass vector
  vector<Type> Rec(n_yrs); // recruitment vector
  array<Type> natmort(n_yrs, n_ages, n_sexes); // array of natural mortality
  Type sigmaR2_early = pow(exp(ln_sigmaR_early),2); // calculate sigmaR2 for early period
  Type sigmaR2_late = pow(exp(ln_sigmaR_late),2); // calculate sigmaR2 for late period
  
  // Fishery Processes
  matrix<Type> Fmort(n_yrs, n_fish_fleets); // matrix of fully selected fishing mortality
  array<Type> FAA(n_yrs, n_ages, n_sexes, n_fish_fleets); // array of fishing mortality at age
  array<Type> TotalFAA(n_yrs, n_ages, n_sexes); // array of total fishing mortality at age
  array<Type> CAA(n_yrs, n_ages, n_sexes, n_fish_fleets); // array of catch at age
  array<Type> CAL(n_yrs, n_lens, n_sexes, n_fish_fleets); // array of catch at length
  matrix<Type> PredCatch(n_yrs, n_fish_fleets); // matrix of total catch by fleet (in weight)
  matrix<Type> PredFishIdx(n_yrs, n_fish_fleets); // matrix of predicted fishery index
  array<Type> fish_sel(n_yrs, n_ages, n_sexes, n_fish_fleets); // array of fishery selectivity
  matrix<Type> fish_q(n_yrs, n_fish_fleets); // matrix of fishery catchability
  array<Type> ESS_FishAgeComps(n_yrs, n_sexes, n_fish_fleets); // array of effective sample sizes for fishery ages (wts * ISS)
  array<Type> ESS_FishLenComps(n_yrs, n_sexes, n_fish_fleets); // array of effective sample sizes for fishery lengths (wts * ISS)

  // Survey Processes
  array<Type> SrvIAA(n_yrs, n_ages, n_sexes, n_srv_fleets); // array of survey index at age
  array<Type> SrvIAL(n_yrs, n_lens, n_sexes, n_srv_fleets); // array of survey index at length
  matrix<Type> PredSrvIdx(n_yrs, n_srv_fleets); // matrix of predicted survey index
  array<Type> srv_sel(n_yrs, n_ages, n_sexes, n_srv_fleets); // array of survey selectivity
  matrix<Type> srv_q(n_yrs, n_srv_fleets); // matrix of survey catchability
  array<Type> ESS_SrvAgeComps(n_yrs, n_sexes, n_srv_fleets); // array of effective sample sizes for survey ages (wts * ISS)
  array<Type> ESS_SrvLenComps(n_yrs, n_sexes, n_srv_fleets); // array of effective sample sizes for survey lengths (wts * ISS)
  
  // Likelihoods
  matrix<Type> Catch_nLL(n_yrs,n_fish_fleets); // Catch likelihoods
  matrix<Type> FishIdx_nLL(n_yrs,n_fish_fleets); // Fishery Index likelihoods
  array<Type> FishAgeComps_nLL(n_yrs,n_sexes,n_fish_fleets); // Fishery Age Compositions likelihoods
  array<Type> FishLenComps_nLL(n_yrs,n_sexes,n_fish_fleets); // Fishery Length Compositions likelihoods
  array<Type> SrvAgeComps_nLL(n_yrs,n_sexes,n_srv_fleets); // Survey Age Compositions likelihoods
  array<Type> SrvLenComps_nLL(n_yrs,n_sexes,n_srv_fleets); // Survey Length Compositions likelihoods
  matrix<Type> SrvIdx_nLL(n_yrs,n_srv_fleets); // Survey Index likelihoods
  vector<Type> Fmort_Pen(n_fish_fleets); // Fishing Mortality Deviation penalty
  vector<Type> Rec_nLL(n_yrs - 1); // recruitment penalty
  vector<Type> Init_Rec_nLL(n_ages - 2); // initial recruitment penalty
  vector<Type> bias_ramp(n_yrs); // bias ramp from Methot and Taylor 2011
  Type M_Pen = 0; // penalty/prior for natural mortality
  Type jnLL = 0; // joint negative log likelihood

  // Initialize at 0
  bias_ramp.setZero();
  SSB.setZero(); 
  CAA.setZero();
  CAL.setZero();
  SrvIAA.setZero();
  SrvIAL.setZero();
  PredCatch.setZero();
  PredFishIdx.setZero();
  PredSrvIdx.setZero();
  fish_q.setZero();
  srv_q.setZero();
  FishIdx_nLL.setZero();
  SrvIdx_nLL.setZero();
  FishAgeComps_nLL.setZero();
  FishLenComps_nLL.setZero();
  SrvAgeComps_nLL.setZero();
  SrvLenComps_nLL.setZero();
  Catch_nLL.setZero();
  Fmort_Pen.setZero();
  Fmort.setZero();
  ESS_FishAgeComps.setZero();
  ESS_SrvAgeComps.setZero();
  ESS_FishLenComps.setZero();
  ESS_SrvLenComps.setZero();
  
  // MODEL SECTION -------------------------------------------------------------
  
  // Process Model(s) ----------------------------------------------------------
  
  // Fishery Selectivity
  for(int y = 0; y < n_yrs; y++) {
    for(int f = 0; f < n_fish_fleets; f++) {
      
      // Extract time-block specification
      int fish_sel_blk_idx = fish_sel_blocks(y,f);
      
      for(int s = 0; s < n_sexes; s++) {
        // extract selectivity parameters
        vector<Type> tmp_fish_sel_vec = ln_fish_fixed_sel_pars.col(f).col(s).col(fish_sel_blk_idx);
        
        for(int a = 0; a < n_ages; a++) {
          Type fish_age_idx = ages(a); // Get age index
          fish_sel(y,a,s,f) = Get_Selex(fish_age_idx, fish_sel_model(y,f), tmp_fish_sel_vec);
        } // end a loop 
        
        // Normalize to max out at 1
        if(normalize_fish_sel(f)) {
          Type tmp_max_fish_sel = max(fish_sel.col(f).col(s).transpose().col(y).vec()); // get max
          for(int a = 0; a < n_ages; a++) fish_sel.col(f).col(s).col(a).col(y) /= tmp_max_fish_sel; // normalize
        } // normalize selectivity to max at 1
          
      } // end s loop
    } // end f loop
  } // end y loop
  
  // Survey Selectivity
  for(int y = 0; y < n_yrs; y++) {
    for(int sf = 0; sf < n_srv_fleets; sf++) {
      
      // Extract time-block specification
      int srv_sel_blk_idx = srv_sel_blocks(y,sf);
      
      for(int s = 0; s < n_sexes; s++) {
        // extract selectivity parameters
        vector<Type> tmp_srv_sel_vec = ln_srv_fixed_sel_pars.col(sf).col(s).col(srv_sel_blk_idx);
        
        for(int a = 0; a < n_ages; a++) {
          Type srv_age_idx = ages(a); // Get age index
          srv_sel(y,a,s,sf) = Get_Selex(srv_age_idx, srv_sel_model(y,sf), tmp_srv_sel_vec);
        } // end a loop 
        
        // Normalize to max out at 1
        if(normalize_srv_sel(sf) == 1) {
          Type tmp_max_srv_sel = max(srv_sel.col(sf).col(s).rotate(1).col(y).vec()); // get max
          for(int a = 0; a < n_ages; a++) srv_sel.col(sf).col(s).col(a).col(y) /= tmp_max_srv_sel; // normalize
        } // normalize selectivity to max at 1

      } // end s loop
    } // end sf loop
  } // end y loop
  
  // Mortality Calculations
  for(int y = 0; y < n_yrs; y++) {
    for(int a = 0; a < n_ages; a++) {
      for(int s = 0; s < n_sexes; s++) {
        
        // Fishing mortality at age calculation
        for(int f = 0; f < n_fish_fleets; f++) {
          if(isNA(ObsCatch(y,f))) { // w/o catch data
            Fmort(y,f) = 0; // set F to zero if not catch data
            FAA(y,a,s,f) = 0; 
          } else {
            Fmort(y,f) = exp(ln_F_mean(f) + ln_F_devs(y,f)); // get fully selected F to report
            FAA(y,a,s,f) = exp(ln_F_mean(f) + ln_F_devs(y,f)) * fish_sel(y,a,s,f); // fishing mortality at age calculation
          } // w / catch data
          TotalFAA(y,a,s) += FAA(y,a,s,f); // get total FAA
        } // end f loop
        
        // Calculate population mortality/survival
        if(s == 0) natmort(y,a,s) = exp(ln_M); // get natural mortality (females)
        if(s == 1) natmort(y,a,s) = exp(ln_M) + M_offset; // get natural mortality (males)
        ZAA(y,a,s) = TotalFAA(y,a,s) + natmort(y,a,s); // Total FAA and natmort to get ZAA
        SAA(y,a,s) = exp(-1.0 * ZAA(y,a,s)); // Calculate survival fraction annually
        SAA_mid(y,a,s) = exp(-0.5 * ZAA(y,a,s)); // Calculate survival fraction during midpoint of year
        
      } // end s loop
    } // end a loop
  } // end y loop
  
  // Calculate Methot and Taylor 2011 bias ramp
  for(int y = 0; y < n_yrs; y++) {
    if(do_rec_bias_ramp == 0) bias_ramp(y) = 1; // don't do bias ramp correction
    if(do_rec_bias_ramp == 1) {
      if(y < bias_year(0) || y == bias_year(3)) bias_ramp(y) = 0; // no bias correction during poor data
      if(y >= bias_year(0) && y < bias_year(1)) bias_ramp(y) = 1 * ((y - bias_year(0)) / (bias_year(1) - bias_year(0))); // ascending limb
      if(y >= bias_year(1) && y < bias_year(2)) bias_ramp(y) = 1; // full bias correction
      if(y >= bias_year(2) && y < bias_year(3)) bias_ramp(y) = 1 * (1 - ((y - bias_year(2)) / (bias_year(3) - bias_year(2)))); // descending limb 
    } // if we want to do bias ramp
  } // end y loop
  
  // Initial age-structure
  Type init_F = init_F_prop * exp(ln_F_mean(0)); // get initial age structure F
  for(int a = 1; a < n_ages; a++) { // (minus the first recruitment deviation)
    for(int s = 0; s < n_sexes; s++) {
      // Initialize age structure
      if(a < n_ages - 1) {
        NAA(0, a, s) = exp(ln_R0  + ln_InitDevs(a-1) -
                      (Type(a) * (natmort(0,a,s) + (init_F * fish_sel(0,a,s,0)) ))) * sexratio(s);
      } // other ages (a-1 for initdevs because of TMB indexing)

      if(a == n_ages - 1) {
        NAA(0,n_ages - 1,s) = exp(ln_R0 - (Type(n_ages - 1) * (natmort(0,n_ages - 1,s) + (init_F * fish_sel(0,n_ages - 1,s,0))) )) /
                              (1 - exp(-(natmort(0,n_ages - 1,s) + (init_F * fish_sel(0,n_ages - 1,s,0)) ) )) * sexratio(s);
      } // plus group

    } // end s loop
  } // end a loop
  
  // Annual Recruitment
  for(int y = 0; y < n_yrs; y++) {
    for(int s = 0; s < n_sexes; s++) {
      // Annual mean recruitment
      if(y < sigmaR_switch) NAA(y,0,s) = exp(ln_R0 + ln_RecDevs(y) - (bias_ramp(y) * sigmaR2_early/2)) * sexratio(s); // early period
      if(y >= sigmaR_switch && y < (n_yrs - 1)) NAA(y,0,s) = exp(ln_R0 + ln_RecDevs(y) - (bias_ramp(y) * sigmaR2_late/2)) * sexratio(s); // late period
      if(y == (n_yrs - 1)) NAA(y,0,s) = exp(ln_R0) * sexratio(s); 
      Rec(y) += NAA(y,0,s); // get vector of predicted recruitments
    } // end s loop
  } // end y loop
  
  // Project Population Forward
  for(int y = 0; y < n_yrs; y++) { // (Exponential mortality model)
    for(int a = 0; a < n_ages; a++) {
      for(int s = 0; s < n_sexes; s++) {

        // Project ages and years forward
        if(a < n_ages-1) NAA(y+1,a+1,s) = NAA(y,a,s) * SAA(y,a,s); // not plus group
        if(a == n_ages-1) NAA(y+1,n_ages-1,s) += (NAA(y,n_ages-1,s) * SAA(y,n_ages-1,s)); // plus group

      } // end s loop
    } // end a loop
  } // end y loop

  // Spawning biomass calculation
  for(int y = 0; y < n_yrs; y++) {
    for(int a = 0; a < n_ages; a++) {
      SSB(y) += NAA(y,a,0) * WAA(y,a,0) * MatAA(y,a,0); 
    } // end a loop
  } // end y loop
  
  // Observation Model(s) ------------------------------------------------------
  
  // Fishery Observation Model
  for(int y = 0; y < n_yrs; y++) {
    for(int f = 0; f < n_fish_fleets; f++) {
      
      int fish_q_blk_idx = fish_q_blocks(y,f); // Extract time-block catchability specification
      fish_q(y,f) = exp(ln_fish_q(fish_q_blk_idx,f)); // Input into fishery catchability container
      
      for(int a = 0; a < n_ages; a++) {
        for(int s = 0; s < n_sexes; s++) {
          
          CAA(y,a,s,f) = FAA(y,a,s,f) / ZAA(y,a,s) * NAA(y,a,s) * (1 - exp(-ZAA(y,a,s))); // Catch at age
          for(int l = 0; l < n_lens; l++) CAL(y,l,s,f) += CAA(y,a,s,f) * SizeAgeTrans(y,l,a,s); // Catch at length
          PredCatch(y,f) += CAA(y,a,s,f) * WAA(y,a,s); // Update to get total catch

          // Get predicted fishery index
          if(fish_idx_type(f) == 0) PredFishIdx(y,f) += NAA(y,a,s) * SAA_mid(y,a,s) * fish_sel(y,a,s,f); // abundance
          if(fish_idx_type(f) == 1) {
            // Sablefish specific - for fitting Japanese LL RPW cpue fishery as is implemented in the ADMB assessment (i.e., only using female selex to calculate this)
            if(fish_q_blk_idx == 0 && share_sel == 0) PredFishIdx(y,f) += NAA(y,a,s) * SAA_mid(y,a,s) * fish_sel(y,a,0,f) * WAA(y,a,s); 
            else PredFishIdx(y,f) += NAA(y,a,s) * SAA_mid(y,a,s) * fish_sel(y,a,s,f) * WAA(y,a,s); // for not first time block
          } // biomass
          
        } // end s loop
      } // end a loop
      
      PredFishIdx(y,f) *= fish_q(y,f); // Update to get predicted index with catchability
      
    } // end f loop
  } // end y loop
  
  // Survey Observation Model
  for(int y = 0; y < n_yrs; y++) {
    for(int sf = 0; sf < n_srv_fleets; sf++) {
      
      int srv_q_blk_idx = srv_q_blocks(y,sf); // Extract time-block catchability specification
      srv_q(y,sf) = exp(ln_srv_q(srv_q_blk_idx,sf)); // Input into survey catchability container
      
      for(int a = 0; a < n_ages; a++) {
        for(int s = 0; s < n_sexes; s++) {
          
          SrvIAA(y,a,s,sf) = NAA(y,a,s) * srv_sel(y,a,s,sf); // Survey index at age
          for(int l = 0; l < n_lens; l++) SrvIAL(y,l,s,sf) += SrvIAA(y,a,s,sf) * SizeAgeTrans(y,l,a,s); // Survey index at length
          
          // Get predicted survey index
          if(srv_idx_type(sf) == 0) PredSrvIdx(y,sf) += SrvIAA(y,a,s,sf) * SAA_mid(y,a,s); // abundance
          if(srv_idx_type(sf) == 1) PredSrvIdx(y,sf) += SrvIAA(y,a,s,sf) * SAA_mid(y,a,s) * WAA(y,a,s); // biomass
          
        } // end s loop
      } // end a loop
      
      PredSrvIdx(y,sf) *= srv_q(y,sf); // Update to get predicted index with catchability
      
    } // end sf loop
  } // end y loop
  
  // Likelihood(s) -------------------------------------------------------------
  // Fishery Likelihoods
  for(int y = 0; y < n_yrs; y++) {
    for(int f = 0; f < n_fish_fleets; f++) {
      
      // Fishery Catches -----------------------------------------------------
      if(!isNA(ObsCatch(y,f))) Catch_nLL(y,f) += UseCatch(y,f) * square(log(ObsCatch(y,f)+Catch_Constant(f)) - log(PredCatch(y,f)+Catch_Constant(f)));

      // Fishery Indices -----------------------------------------------------
      if(!isNA(ObsFishIdx(y,f))) FishIdx_nLL(y,f) += UseFishIdx(y,f) * square(log(ObsFishIdx(y,f)+0.0001) - log(PredFishIdx(y,f)+0.0001)) / (2* square(ObsFishIdx_SE(y,f) / ObsFishIdx(y,f)));

      // Aggregating Fishery Age Compositions --------------------------------
      if(AggFishAgeComps(y,f) == 0) { 
        if(!isNA(ObsFishAgeComps.col(f).col(0).rotate(1).col(y).sum())) { // only evaluate when there aren't NAs in data
          
          // Extract out quantities and reformat
          vector<Type> tmp_ObsFishAgeComps = ObsFishAgeComps.col(f).col(0).rotate(1).col(y); // Get proportions (dim = n_ages, 1)
          matrix<Type> tmp_CAA = CAA.col(f).transpose().col(y).transpose().matrix(); // Extract out CAA to aggregate (dim = n_ages, n_sexes)
          
          for(int s = 0; s < n_sexes; s++) { // normalize proportions within a sex and then sum 
            Type total_pred_age_comps = tmp_CAA.col(s).sum(); // get total predicted comps for a given sex
            tmp_CAA.col(s) /= total_pred_age_comps; // normalize by total within a given sex
          } // end s loop
          
          // Sum across rows (ages), divide to get aggregated sexes, and apply ageing error
          vector<Type> tmp_Agg_CAA = (tmp_CAA.rowwise().sum() / n_sexes).transpose() * AgeingError;
          // Calculate effective sample size for fishery age compositions
          ESS_FishAgeComps(y,0,f) = ISS_FishAgeComps(y,0,f) * Wt_FishAgeComps(0,f); 
          
          // Compute ADMB Multinomial Likelihood
          FishAgeComps_nLL(y,0,f) -= ESS_FishAgeComps(y,0,f) * UseFishAgeComps(y,f) * ((tmp_ObsFishAgeComps + 0.001) * log(tmp_Agg_CAA + 0.001)).sum();

        } // if there are age composition data for a given fleet
      } // if aggregating fishery age comps
      
      // Fishery Length Compositions -----------------------------------------
      for(int s = 0; s < n_sexes; s++) {
        if(!isNA(ObsFishLenComps.col(f).col(s).rotate(1).col(y).sum())) { // only evaluate when there aren't NAs in data

        // Extract out quantities and reformat
        vector<Type> tmp_ObsFishLenComps = ObsFishLenComps.col(f).col(s).rotate(1).col(y); // Get observed proportions (dim = n_lens)
        vector<Type> tmp_CAL = CAL.col(f).col(s).rotate(1).col(y); // Extract out CAL (dim = n_lens)
        Type total_pred_len_comps = tmp_CAL.sum(); // get total predicted comps for a given sex
        tmp_CAL /= total_pred_len_comps; // normalize by total within a given sex

        // Calculate effective sample size for fishery length compositions
        ESS_FishLenComps(y,s,f) = ISS_FishLenComps(y,s,f) * Wt_FishLenComps(s,f);

        // Compute ADMB Multinomial Likelihood
        FishLenComps_nLL(y,s,f) -= ESS_FishLenComps(y,s,f) * UseFishLenComps(y,f) * ((tmp_ObsFishLenComps + 0.001) * log(tmp_CAL + 0.001)).sum();

        } // if len comps for a given fleet
      } // end s loop
      
    } // end f loop
  } // end y loop
  
  // Survey Likelihoods
  for(int y = 0; y < n_yrs; y++) {
    for(int sf = 0; sf < n_srv_fleets; sf++) {
      
      // Survey Indices ------------------------------------------------------
      // ADMB likelihood
      if(!isNA(ObsSrvIdx(y,sf))) SrvIdx_nLL(y,sf) += UseSrvIdx(y,sf) * square(log(ObsSrvIdx(y,sf)+0.0001) - log(PredSrvIdx(y,sf)+0.0001)) / (2 * square(ObsSrvIdx_SE(y,sf) / ObsSrvIdx(y,sf)));

      // Aggregating Survey Age Compositions ---------------------------------
      if(AggSrvAgeComps(y,sf) == 0) { 
        if(!isNA(ObsSrvAgeComps.col(sf).col(0).rotate(1).col(y).sum())) { // only evaluate when there aren't NAs in data
          
          // Extract out quantities and reformat
          vector<Type> tmp_ObsSrvAgeComps = ObsSrvAgeComps.col(sf).col(0).rotate(1).col(y); // Get proportions (dim = n_ages, 1)
          matrix<Type> tmp_SrvIAA = SrvIAA.col(sf).transpose().col(y).transpose().matrix(); // Extract out survey index at age to aggregate (dim = n_ages, n_sexes)
          
          for(int s = 0; s < n_sexes; s++) { // normalize proportions within a sex and then sum 
            Type total_pred_age_comps = tmp_SrvIAA.col(s).sum(); // get total predicted comps for a given sex
            tmp_SrvIAA.col(s) /= total_pred_age_comps; // normalize by total within a given sex
          } // end s loop
          
          // Sum across rows (ages), divide to get aggregated sexes, and apply ageing error
          vector<Type> tmp_Agg_SrvIAA = (tmp_SrvIAA.rowwise().sum() / n_sexes).transpose() * AgeingError;
          // Calculate effective sample size for fishery age compositions
          ESS_SrvAgeComps(y,0,sf) = ISS_SrvAgeComps(y,0,sf) * Wt_SrvAgeComps(0,sf); 
          
          // Compute Multinomial Likelihood
          SrvAgeComps_nLL(y,0,sf) -= ESS_SrvAgeComps(y,0,sf) * UseSrvAgeComps(y,sf) * ((tmp_ObsSrvAgeComps + 0.001) * log(tmp_Agg_SrvIAA + 0.001)).sum(); 
     
        } // if there are age composition data for a given survey fleet
      } // if aggregating survey age comps
    } // end sf loop
  } // end y loop
  
  for(int y = 0; y < n_yrs; y++) {
    for(int sf = 0; sf < n_srv_fleets; sf++) {
      // Survey Length Compositions ------------------------------------------
      for(int s = 0; s < n_sexes; s++) {
        if(!isNA(ObsSrvLenComps.col(sf).col(s).rotate(1).col(y).sum())) { // only evaluate when there aren't NAs in data
          
          // Extract out quantities and reformat
          vector<Type> tmp_ObsSrvLenComps = ObsSrvLenComps.col(sf).col(s).rotate(1).col(y); // Get observed proportions (dim = n_lens)
          vector<Type> tmp_SrvIAL = SrvIAL.col(sf).col(s).rotate(1).col(y); // Extract out survey index at length (dim = n_lens)
          Type total_pred_len_comps = tmp_SrvIAL.sum(); // get total predicted comps for a given sex
          tmp_SrvIAL /= total_pred_len_comps; // normalize by total within a given sex
          
          // Calculate effective sample size for survey length compositions
          ESS_SrvLenComps(y,s,sf) = ISS_SrvLenComps(y,s,sf) * Wt_SrvLenComps(s,sf);
          
          // Compute Multinomial Likelihood
          SrvLenComps_nLL(y,s,sf) -= ESS_SrvLenComps(y,s,sf) * UseSrvLenComps(y,sf) * ((tmp_ObsSrvLenComps + 0.001) * log(tmp_SrvIAL + 0.001)).sum();

        } // if len comps for a given fleet
      } // end s loop
    } // end sf loop
  } // end y loop

  // Priors and Penalties
  // Fishing Mortality (Penalty)
  for(int f = 0; f < n_fish_fleets; f++) {
    for(int y = 0; y < n_yrs; y++) {
      // Only penalize when there are catch data
      if(!isNA(ObsCatch(y,f))) Fmort_Pen(f) += square(ln_F_devs(y,f)); // SSQ penalize
    } // end y loop
  } // end f loop
  
  // Natural Mortality (Prior)
  if(Use_M_prior == 1) M_Pen = square(ln_M - log(M_prior(0))) / (2 * square(M_prior(1)));

  // Recruitment (Penalty)
  for(int a = 0; a < n_ages - 2; a++) Init_Rec_nLL(a) = square(ln_InitDevs(a) / exp(ln_sigmaR_early)); // initial age structure
  for(int y = 0; y < sigmaR_switch; y++) Rec_nLL(y) = square(ln_RecDevs(y)/exp(ln_sigmaR_early)) + bias_ramp(y)*ln_sigmaR_early; // early period
  for(int y = sigmaR_switch; y < (n_yrs - 1); y++) Rec_nLL(y) = square(ln_RecDevs(y)/exp(ln_sigmaR_late)) + bias_ramp(y)*ln_sigmaR_late; // late period
  
  // Apply likelihood weights here and compute joint negative log likelihood
  jnLL = (Wt_Catch* Catch_nLL.sum()) + // Catch likelihoods
         (Wt_FishIdx * FishIdx_nLL.sum()) + // Fishery Index likelihood
         FishAgeComps_nLL.sum() + // Fishery Age likelihood
         FishLenComps_nLL.sum() + // Fishery Length likelihood
         (Wt_SrvIdx * SrvIdx_nLL.sum()) + // Survey Index likelihood
         SrvAgeComps_nLL.sum() + // Survey Age likelihood
         SrvLenComps_nLL.sum() + // Survey Length likelihood
         (Wt_F* Fmort_Pen.sum()) + // Fishery Mortality Penalty
         M_Pen + // Natural Mortality Prior (Penalty)
         (Wt_Rec * 0.5 * Rec_nLL.sum()) + // Recruitment Penalty
         (Wt_Rec * 0.5 * Init_Rec_nLL.sum()); // Initial Age Penalty

  // REPORT SECTION ------------------------------------------------------------
  // Biological Processes
  REPORT(NAA);
  REPORT(ZAA);
  REPORT(natmort); 
  REPORT(bias_ramp);
  
  // Fishery Processes
  REPORT(Fmort); 
  REPORT(FAA);
  REPORT(CAA);
  REPORT(CAL);
  REPORT(PredCatch);
  REPORT(PredFishIdx);
  REPORT(fish_sel);
  REPORT(fish_q);
  
  // Survey Processes
  REPORT(PredSrvIdx);
  REPORT(srv_sel);
  REPORT(srv_q);
  REPORT(SrvIAA);
  REPORT(SrvIAL);
  
  // Likelihoods
  REPORT(Catch_nLL); 
  REPORT(FishIdx_nLL); 
  REPORT(SrvIdx_nLL); 
  REPORT(FishAgeComps_nLL); 
  REPORT(SrvAgeComps_nLL); 
  REPORT(FishLenComps_nLL); 
  REPORT(SrvLenComps_nLL); 
  REPORT(M_Pen);
  REPORT(Fmort_Pen); 
  REPORT(Rec_nLL);
  REPORT(Init_Rec_nLL);
  REPORT(Rec_nLL);
  REPORT(jnLL);
  
  // Effective Sample Sizes
  REPORT(ESS_FishAgeComps); 
  REPORT(ESS_SrvAgeComps); 
  REPORT(ESS_FishLenComps); 
  REPORT(ESS_SrvLenComps); 
  
  // Report for derived quantities
  REPORT(SSB);
  REPORT(Rec); 
  ADREPORT(SSB);
  ADREPORT(Rec);
  
  return jnLL; 
} // end objective function
