# Data weighting ----------------------------------------------------------
# Francis reweighting
for(i in 1:10) {

  if(i == 1) { # set weights to 1 at first iteration
    data$Wt_FishAgeComps[!is.na(data$Wt_FishAgeComps)] <- 1
    data$Wt_FishLenComps[!is.na(data$Wt_FishLenComps)] <- 1
    data$Wt_SrvAgeComps[!is.na(data$Wt_SrvAgeComps)] <- 1
    data$Wt_SrvLenComps[!is.na(data$Wt_SrvLenComps)] <- 1
  } # if i == 1 set weights to 1

  # make AD model function
  sabie_model <- MakeADFun(data = data, parameters = parameters,
                           map = mapping, random = NULL,
                           DLL = "SabieTMB")
  # Now, optimize the function
  sabie_optim <- stats::nlminb(sabie_model$par, sabie_model$fn, sabie_model$gr,
                               control = list(iter.max = 1e5, eval.max = 1e5))
  # newton steps
  try_improve <- tryCatch(expr =
                            for(i in 1:2) {
                              g = as.numeric(sabie_model$gr(sabie_optim$par))
                              h = optimHess(sabie_optim$par, fn = sabie_model$fn, gr = sabie_model$gr)
                              sabie_optim$par = sabie_optim$par - solve(h,g)
                              sabie_optim$objective = sabie_model$fn(sabie_optim$par)
                            }
                          , error = function(e){e}, warning = function(w){w})

  sabie_model$optim <- sabie_optim # Save optimized model results
  sabie_model$rep <- sabie_model$report(sabie_model$env$last.par.best) # Get report

  # Get new francis weights
  new_weights <- francis_rwgt(data = data, model = sabie_model)

  data$Wt_FishAgeComps <- new_weights$new_fish_age_wts
  data$Wt_FishLenComps <- new_weights$new_fish_len_wts
  data$Wt_SrvAgeComps <- new_weights$new_srv_age_wts
  data$Wt_SrvLenComps <- new_weights$new_srv_len_wts
} # end i loop

sabie_model$sd_rep <- sdreport(sabie_model)