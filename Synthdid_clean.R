# This code file contains all the Synthetic Control algorithms as well
# as functions to calculate the corresponding estimators.


# Import the necessary optimization package
library(Rsolnp)

## Synthetic Difference in differences method
synthdid <- function(Y0, Y1, T0){
  
  # Inputs:
  ## Y0: matrix of control unit outcomes
  ## Y1: vector of treatment unit outcomes
  ## T0: time period of intervention
  
  # Output:
  ## Optimal unit weights: unit.w
  ## Optimal time weights: time.w
  
  # Here is the function for the unit weights
  sdid.unit.w <- function( pars = rep(1/(dim(Y0)[2]), dim(Y0)[2])){
    Y0 = as.matrix(Y0)
    
    # We first need to obtain the 1st difference matrix:
    Y = cbind(Y0,Y1)[1:(T0-1),]
    diffY = apply(Y, 2, diff)
    Tpre = T0-1
    Tpost = dim(Y0)[1] - Tpre
    Nco = dim(Y0)[2]
    # First, we calculate zeta
    delta.bar <- 1/(Nco*(Tpre-1))*sum(apply(diffY, 2, sum))
    sig.bar <- 1/(Nco*(Tpre-1))*sum(apply((diffY-delta.bar)^2, 2, sum))
    zeta <- (1*Tpost)^(1/4)*sqrt(sig.bar) #Here, we assume only 1 treated unit.
    
    # Constructing the pretreatment periods
    Y0pre = Y0[1:(T0-1),]
    Y1pre = Y1[1:(T0-1)]
    
    # Now we define the unit weights optimization function:
    unit.w <- function(w){
      # We first need to detrend the matrix: 
      # Note that we need the pretreatment part of the matrix
      Y0_detrend_over_time <- t(t(Y0pre) -  apply(Y0pre, 2, mean))
      Y1_detrend_over_time <- Y1pre - mean(Y1pre)
      
      # Now we can define the function, on which to optimize:
      sse <- sum( (Y0_detrend_over_time %*% w - Y1_detrend_over_time)^2   )
      return(sse + zeta^2 * Tpre * sum(w^2) )
    }
    eq = function(x){
      return(sum(x))
    }
    
    opt = solnp(pars,
                unit.w,
                eqfun = eq,
                eqB = 1,
                LB=rep(0, length(pars)),
                UB=rep(1,length(pars)))
    return(opt$pars)
  }
  
  # Function to optimize time weights
  sdid.time.w <- function(pars = rep(1/(T0-1), T0-1)){
    time.w <- function(lam){
      # We first need to detrend the matrix: 
      # Note that we need the pretreatment part of the matrix
      Y0_detrend_over_unit <- Y0 - apply(Y0, 1, mean)
      
      # Now we can define the function, on which to optimize:
      T = dim(Y0_detrend_over_unit)[1]
      sse <- sum( (t(Y0_detrend_over_unit[1:(T0-1),]) %*% lam - 
                     apply(Y0_detrend_over_unit[T0: T,],2,mean) )^2)
      return(sse)
    }
    
    eq = function(x){
      return(sum(x))
    }
    
    opt = solnp(pars,
                time.w,
                eqfun = eq,
                eqB = 1,
                LB=rep(0, length(pars)),
                UB=rep(1,length(pars)))
    return(opt$pars)
  }
  w <- sdid.unit.w()
  lam <- sdid.time.w()
  return(list("unit.w"=w, "time.w"=lam))
}


######################################################
#######Synthetic with no covariates (original version)
synth.no.cov <- function(Y0, Y1, T0, pars = rep(1/(dim(Y0)[2]), dim(Y0)[2])){
  
  # Inputs:
  ## Y0: matrix of control unit outcomes (could be a dataframe, because we coerce to matrix below)
  ## Y1: vector of treatment unit outcomes
  ## T0: time period of intervention
  ## pars: initial values for optimization
  
  # Output:
  ## Optimal unit weights
  
  
  Y0 = as.matrix(Y0)
  # constructing pre-treatment time series:
  Y0_pre = Y0[1:(T0-1),]
  Y1_pre = Y1[1:(T0-1)]
  
  eq = function(x){
    return(sum(x))
  }
  
  w.fnc <- function(w){
    sum((Y0_pre %*% w - Y1_pre)^2)
  }
  
  opt = solnp(pars,
              w.fnc,
              eqfun = eq,
              eqB = 1,
              LB=rep(0, length(pars)),
              UB=rep(1,length(pars)))
  return(opt$pars)
}


##############################################
#Demeaned Synthetic control without covariates 
synth_demean <- function( Y0, Y1,T0, pars = rep(1/(dim(Y0)[2]), dim(Y0)[2])){
  
  # Inputs:
  ## Y0: matrix of control unit outcomes (could be a dataframe, because we coerce to matrix below)
  ## Y1: vector of treatment unit outcomes
  ## T0: time period of intervention
  ## pars: initial values for optimization
  
  # Output:
  ## Optimal unit weights
  
  Y0 = as.matrix(Y0)
  
  # Constructing the pretreatment periods
  Y0pre = Y0[1:(T0-1),]
  Y1pre = Y1[1:(T0-1)]
  
  # Now we define the unit weights optimization function:
  unit.w <- function(w){
    # We first need to detrend the matrix: 
    # Note that we need the pretreatment part of the matrix
    Y0_detrend_over_time <- t(t(Y0pre) -  apply(Y0pre, 2, mean))
    Y1_detrend_over_time <- Y1pre - mean(Y1pre)
    
    # Now we can define the function, on which to optimize:
    sse <- sum( (Y0_detrend_over_time %*% w - Y1_detrend_over_time)^2   )
    return(sse)
  }
  eq = function(x){
    return(sum(x))
  }
  
  opt = solnp(pars,
              unit.w,
              eqfun = eq,
              eqB = 1,
              LB=rep(0, length(pars)),
              UB=rep(1,length(pars)))
  return(opt$pars)
}

#######################################################################################
### Down below, we write functions to calculate the SC estimators 

# w : weights supplied
# Y1: full vector of treated unit
# Y0: full matrix (or data frame) of control units
# T0: first time of intervention
# indiv: whether we would like to look at individual post-intervention effect
   # indiv default = F => calculate average post treatment effect
# t1 (only pass in if indiv = T): one specific time period of interest to calculate intervention effect.

# 1 )For the SC
est.sc <- function(w, Y1, Y0, T0, indiv = F, t1 = NULL){
  Y0 = as.matrix(Y0)
  if (indiv == F){
    Total = length(Y1)
    return( mean(Y1[T0:Total]) - sum( w * apply(Y0[T0:Total,],2,mean)))
  } else {
    return( Y1[t1] - sum( w*Y0[t1,]) )
  }
}

# 2) For the SD
est.sd <- function(w, Y1, Y0, T0, indiv = F, t1 = NULL){
  Y0 = as.matrix(Y0)
  if (indiv == F){
    Total = length(Y1)
    return( (mean(Y1[T0:Total]) - mean(Y1[1:(T0-1)]) ) -
              sum( w * (apply(Y0[T0:Total,],2,mean) -apply(Y0[1:(T0-1),],2,mean) )   ))
  } else {
    return( (Y1[t1] - mean(Y1[1:(T0-1)])) -
              sum( w * ( Y0[t1,] -apply(Y0[1:(T0-1),],2,mean) ) )  )
  }
} 

# 3) For the SDID
est.sdid <- function(w,lam, Y1, Y0, T0, indiv = F, t1 = NULL){
  Y0 = as.matrix(Y0)
  if (indiv == F){
    Total = length(Y1)
    return( (mean(Y1[T0:Total]) - sum(lam*Y1[1:(T0-1)]) ) -
              sum( w * (apply(Y0[T0:Total,],2,mean) - t(Y0[1:(T0-1),]) %*% lam )   ))
  } else {
    return( (Y1[t1] - sum(lam*Y1[1:(T0-1)]) ) -
              sum( w * ( Y0[t1,] -t(Y0[1:(T0-1),]) %*% lam ) )  )
  }
}




