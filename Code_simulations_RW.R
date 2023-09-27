# In this file, we generate the simulated dataset and save them into .txt files 
# for later use.
# Here, we consider the random walk time effects (T2) case.

#################################################################
# Set the appropriate (local) working directory:
setwd("~/UvA-3rd year/Course_material_sem_2/B5+6 - Thesis/Code")

# Remember to put the following source code in the same folder:
source("Synthdid_clean.R")


#################################################################
# Define the specifications here
post_treatment = 12
T0 = c(12,32,100)
Time = T0 + post_treatment # Total number of time periods
N = c(10,50,150) # Number of units # Change here the number of units 
sig = c(1, sqrt(10))
beta = 1 #Parameter of the exponential distribution
time_effects_param = 1 # Parameter of the exponential distribution for time-effects
iter = 250 # Number of iterations for each specification
#################################################################


# Load in the DGPs
source("DGP.R")

#############################################################


len.sig = length(sig)
len.Time = length(Time)
len.units = length(N)

# We will let a, b, c be the iterators of sig, Time, and units respectively

for (a in 1:len.sig){
  for (b in 1:len.Time){
    for (c in 1:len.units){
      SCestimates = matrix(rep(0, 3*iter), ncol = 3)
      for (i in 1:iter){
        data = gen_data.RW(Total = Time[b], N= N[c],T0= T0[b],
                        sd=sig[a], tau = 5, time_effects_param, beta = beta)
        
        w.sc = synth.no.cov(data$Y0, data$Y1, T0 =T0[b])
        w.sd = synth_demean(data$Y0, data$Y1, T0 =T0[b])
        w.sdid = synthdid(data$Y0, data$Y1, T0 =T0[b])
        
        tau_sc = est.sc(w.sc, data$Y1, data$Y0, T0 =T0[b])
        tau_sd = est.sd(w.sd, data$Y1, data$Y0, T0 =T0[b])
        tau_sdid = est.sdid(w.sdid$unit.w, w.sdid$time.w,data$Y1, data$Y0, T0 =T0[b])
        
        SCestimates[i,] = c(tau_sc, tau_sd, tau_sdid)
      }
      
      fname = paste( "RW_M1,", "T=", paste(Time[b]), ",Units=",
                     paste(N[c]), ",eps_sig=", paste(round(sig[a],2)), ".csv", sep = "")
      
      file.create(fname)
      write.table(x = SCestimates, file = fname)
    }
  }
}


for (a in 1:len.sig){
  for (b in 1:len.Time){
    for (c in 1:len.units){
      SCestimates = matrix(rep(0, 3*iter), ncol = 3)
      for (i in 1:iter){
        data = gen_data.RW.M2(Total = Time[b], N= N[c],T0= T0[b],
                           sd=sig[a], tau = 5, time_effects_param, beta = beta)
        
        w.sc = synth.no.cov(data$Y0, data$Y1, T0 =T0[b])
        w.sd = synth_demean(data$Y0, data$Y1, T0 =T0[b])
        w.sdid = synthdid(data$Y0, data$Y1, T0 =T0[b])
        
        tau_sc = est.sc(w.sc, data$Y1, data$Y0, T0 =T0[b])
        tau_sd = est.sd(w.sd, data$Y1, data$Y0, T0 =T0[b])
        tau_sdid = est.sdid(w.sdid$unit.w, w.sdid$time.w,data$Y1, data$Y0, T0 =T0[b])
        
        SCestimates[i,] = c(tau_sc, tau_sd, tau_sdid)
      }
      
      fname = paste( "RW_M2,", "T=", paste(Time[b]), ",Units=",
                     paste(N[c]), ",eps_sig=", paste(round(sig[a],2)), ".csv", sep = "")
      
      file.create(fname)
      write.table(x = SCestimates, file = fname)
    }
  }
}



