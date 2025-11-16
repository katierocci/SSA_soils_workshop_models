#Title: run_functions.R
#Author: Rose Abramoff
#Date: Nov 16, 2025

##This function takes as arguments:
##1 input data
##2 model equations
##3 parameters
##4 whether or not to calculate the eigenvalues (optional; default = 0; 0 = no, 1 = yes)
##5 initial states of pools (optional; default = 1 for all pools and 0 for CO2 flux)
Solve_Model  <- function(inputdata,derivs,parameters,state=c(StrLitter=1, MetLitter=1, ACTIVE=1, SLOW=1, PASSIVE=1)) {

    SStime = 1000*365 #time at which steady state solution is evaluated

    #Define constant forcing functions
    forc_st <- approxfun(1:SStime, rep(mean(inputdata$forc_st),SStime)) #temperature input function
    forc_sw <- approxfun(1:SStime, rep(mean(inputdata$forc_sw),SStime)) #moisture input function
    forc_npp <- approxfun(1:SStime, rep(mean(inputdata$forc_npp),SStime)) #NPP input function

    derivs_SS_wrapper <- function(step.num,state,parameters,forc_st,forc_sw,forc_npp) {
      output <- derivs(step.num,state,parameters,forc_st,forc_sw,forc_npp)
      return(list(output[[1]][1:5]))
    }
    
    #Solve steady state
    SS_output <- stode(y = state, time = SStime, func = derivs_SS_wrapper, parms = parameters, forc_st=forc_st, forc_sw=forc_sw, forc_npp=forc_npp, positive = TRUE)
    return(SS_output)
}

Solve_Model_for_Row <- function(row, parameters, state = c(StrLitter=1, MetLitter=1, ACTIVE=1, SLOW=1, PASSIVE=1)) {
  SStime <- 1000 * 365
  
  forc_st <- approxfun(1:SStime, rep(row$forc_st, SStime))
  forc_sw <- approxfun(1:SStime, rep(row$forc_sw, SStime))
  forc_npp <- approxfun(1:SStime, rep(row$forc_npp, SStime))
  
  # Override default site-specific value with values from AfSIS
  parameters$input_to_metlitter <- as.numeric(row$input_to_metlitter)
  parameters$param_claysilt <- as.numeric(row$param_claysilt)
  parameters$LigFrac <- as.numeric(row$Ls)
  
  derivs_SS_wrapper <- function(step.num, state, parameters, forc_st, forc_sw, forc_npp) {
    
    output <- derivs_Century(step.num, state, parameters, forc_st, forc_sw, forc_npp)
    
    return(list(output[[1]][1:5]))
  }
  
  state.SS <- stode(
    y = state,
    time = SStime,
    func = derivs_SS_wrapper,
    parms = parameters,
    forc_st = forc_st,
    forc_sw = forc_sw,
    forc_npp = forc_npp,
    positive = TRUE
  )
  
  modeled <- as.data.frame(cbind(row$SSN_row_ID, 
                                 state.SS$y[1], 
                                 state.SS$y[2], 
                                 state.SS$y[3], 
                                 state.SS$y[4], 
                                 state.SS$y[5], 
                                 sum(state.SS$y[1:5])))
  names(modeled) <- c("SSN_row_ID", "StrLitter", "MetLitter", "ACTIVE", "SLOW", "PASSIVE", "Soil_Organic_Carbon_kg_m2")
  
  #re-define type as numeric and convert from g/m2 to kg/m2
  modeled$StrLitter <- as.numeric(modeled$StrLitter)/1000
  modeled$MetLitter <- as.numeric(modeled$MetLitter)/1000
  modeled$ACTIVE <- as.numeric(modeled$ACTIVE)/1000
  modeled$SLOW <- as.numeric(modeled$SLOW)/1000
  modeled$PASSIVE <- as.numeric(modeled$PASSIVE)/1000
  modeled$Soil_Organic_Carbon_kg_m2 <- as.numeric(modeled$Soil_Organic_Carbon_kg_m2)/1000
  
  return(modeled)
}

Solve_Model_for_Fit <- function(inputdata, parameters) {
  row_list <- split(inputdata, seq_len(nrow(inputdata)))
  
  modeled_pools_list <- lapply(row_list, function(row) {
    model_out <- Solve_Model_for_Row(row, parameters)
  })
  
  modeled.pools <- do.call(rbind, modeled_pools_list)
  
  rownames(modeled.pools) <- NULL
  return(modeled.pools)
  
}
