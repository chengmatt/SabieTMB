# Data weighting ----------------------------------------------------------
fishage_wts = vector()
# Francis reweighting
for(i in 1:5) {

  if(i == 1) { # set weights to 1 at first iteration
    data$Wt_FishAgeComps[!is.na(data$Wt_FishAgeComps)] <- 1
    data$Wt_FishLenComps[!is.na(data$Wt_FishLenComps)] <- 1
    data$Wt_SrvAgeComps[!is.na(data$Wt_SrvAgeComps)] <- 1
    data$Wt_SrvLenComps[!is.na(data$Wt_SrvLenComps)] <- 1
  } # if i == 1 set weights to 1

  # make AD model function
  sabie_rtmb_model <- RTMB::MakeADFun(sabie_RTMB, parameters = parameters, map = mapping)
  
  # Now, optimize the function
  sabie_optim <- stats::nlminb(sabie_rtmb_model$par, sabie_rtmb_model$fn, sabie_rtmb_model$gr,
                               control = list(iter.max = 1e5, eval.max = 1e5))
  # newton steps
  try_improve <- tryCatch(expr =
                            for(n in 1:2) {
                              g = as.numeric(sabie_rtmb_model$gr(sabie_optim$par))
                              h = optimHess(sabie_optim$par, fn = sabie_rtmb_model$fn, gr = sabie_rtmb_model$gr)
                              sabie_optim$par = sabie_optim$par - solve(h,g)
                              sabie_optim$objective = sabie_rtmb_model$fn(sabie_optim$par)
                            }
                          , error = function(e){e}, warning = function(w){w})

  sabie_rtmb_model$optim <- sabie_optim # Save optimized model results
  sabie_rtmb_model$rep <- sabie_rtmb_model$report(sabie_rtmb_model$env$last.par.best) # Get report

  # Get new francis weights
  new_weights <- francis_rwgt(data = data, model = sabie_rtmb_model)

  data$Wt_FishAgeComps <- new_weights$new_fish_age_wts
  data$Wt_FishLenComps <- new_weights$new_fish_len_wts
  data$Wt_SrvAgeComps <- new_weights$new_srv_age_wts
  data$Wt_SrvLenComps <- new_weights$new_srv_len_wts
} # end i loop

sabie_rtmb_model$sd_rep <- sdreport(sabie_rtmb_model)
