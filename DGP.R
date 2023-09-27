# In this file, we have the functions for the DGP used within the paper.
# Basically, given the parameters, the function will output simulated data for
# Y1: the time series of outcomes of the first unit (treatment unit)
##    Output is stored as a vector
# Y0: a matrix whose columns represents the time series of the outcomes of 
# control units.
##    Output is stored as a matrix

# U1, T1 case.
gen_data <- function(Total, N, T0,sd, tau, time_var_param, beta){
  # There are 7 parameters to fill in:
  
  # Total: total number of time periods to simulate
  # N : number of units (including the treated unit)
  # T0 : first time period of treatment effect
  # sd: Standard deviation of the noise term (here it is normal distributed)
  # tau: constant treatment effect
  # time_var_param: time effect ~ EXP (time_var_param)
  # beta: unit effect ~ EXP(beta)
  
  Y = matrix(rep(0, Total*N), nrow = Total)
  unit_params = rexp(N, rate = beta) 
  for (t in 1:Total){
    # Sampling the time effect for period t:
    lam.t = rexp(1, rate = time_var_param)
    # Sampling the noise term for each unit (e_{it})
    eps = rnorm(N, sd= sd)

    Y[t,] = rep(lam.t, N)*unit_params + eps
  }
  # Adding the treatment effect for the treated unit
  Y[,1] = Y[,1] + c(rep(0, T0-1), rep(tau, Total-T0+1))
  return(list("Y1" = Y[,1], "Y0" = Y[,2:N]))
}

# Below is the case where the unit effect is outside the convex hull of control units
# Here mu_1 = max{mu} + 1
# U2, T1 case
gen_data.M2 <- function(Total, N, T0, sd, tau, time_var_param, beta){
  # There are 7 parameters to fill in (similar to those of previous function)
  Y = matrix(rep(0, Total*N), nrow = Total)
  unit_params = c(0,rexp(N-1, rate = beta))
  unit_params[1] = max(unit_params) + 1
  for (t in 1:Total){
    # Sampling time effects
    lam.t = rexp(1, rate = time_var_param)
    # Sampling the noise terms
    eps = rnorm(N, sd= sd)
   
    Y[t,] = rep(lam.t, N)*unit_params + eps
  }
  Y[,1] = Y[,1] + c(rep(0, T0-1), rep(tau, Total-T0+1))
  return(list("Y1" = Y[,1], "Y0" = Y[,2:N]))
}



#######
## Below we also include random walk with drift time effects
# U1, T2 case
gen_data.RW <- function(Total, N, T0,sd, tau, beta, alpha = 0.1){
  # There are 6 parameters to fill in:
  
  # Total: total number of time periods to simulate
  # N : number of units (including the treated unit)
  # T0 : first time period of treatment effect
  # sd: Standard deviation of the noise term (here it is normal distributed)
  # tau: constant treatment effect
  # beta: unit effect ~ EXP(beta)
  # alpha is the drift, default to 0.1
  
  # Note that before we had time effects that are exponentially distributed, 
  # but here, they will follow a random walk with drift
  
  Y = matrix(rep(0, Total*N), nrow = Total)
  unit_params = rexp(N, rate = beta) 
  print(unit_params)
  #Initializing time effects:
  lam.t = runif(1)
  for (t in 1:Total){
    # Note that we set the noise of the RW not too high, otherwise, 
    # control units will not "follow" or "be affected" by the same
    # trends as the treated unit.
    lam.t = alpha + lam.t + rnorm(1, sd=0.2)
    # Sampling the noise term for each unit (e_{it})
    eps = rnorm(N, sd= sd)
    
    Y[t,] = rep(lam.t, N)*unit_params + eps
  }
  # Adding the treatment effect for the treated unit
  Y[,1] = Y[,1] + c(rep(0, T0-1), rep(tau, Total-T0+1))
  return(list("Y1" = Y[,1], "Y0" = Y[,2:N]))
}


# U2, T2 case data
gen_data.RW.M2 <- function(Total, N, T0,sd, tau, beta, alpha = 0.1){
  # There are 6 parameters to fill in:
  
  # Total: total number of time periods to simulate
  # N : number of units (including the treated unit)
  # T0 : first time period of treatment effect
  # sd: Standard deviation of the noise term (here it is normal distributed)
  # tau: constant treatment effect
  # beta: unit effect ~ EXP(beta)
  # alpha is the drift, default to 0.2
  
  # Note that before we had time effects that are exponentially distributed, 
  # but here, they will follow a random walk with drift
  
  Y = matrix(rep(0, Total*N), nrow = Total)
  unit_params = c(0,rexp(N-1, rate = beta))
  unit_params[1] = max(unit_params) + 1 
  print(unit_params)
  #Initializing time effects:
  lam.t = runif(1)
  for (t in 1:Total){
    # Note that we set the noise of the RW not too high, otherwise, 
    # control units will not "follow" or "be affected" by the same
    # trends as the treated unit.
    lam.t = alpha + lam.t + rnorm(1, sd=0.2)
    # Sampling the noise term for each unit (e_{it})
    eps = rnorm(N, sd= sd)
    
    Y[t,] = rep(lam.t, N)*unit_params + eps
  }
  # Adding the treatment effect for the treated unit
  Y[,1] = Y[,1] + c(rep(0, T0-1), rep(tau, Total-T0+1))
  return(list("Y1" = Y[,1], "Y0" = Y[,2:N]))
}


