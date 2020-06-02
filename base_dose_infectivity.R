## This runs simulations where dose-infectivity relationships can be 
## accelerating, linear, or decelerating, but where dose-mortality and dose-excretion
## relationships are absent. We record the density of susceptible and infected
## individuals of two host species, as well as pathogen density, as we increase 
## the growth rate of wone of the host species. 

##This creates results from initial figure, as well as supplemetal figures (SI3)

#Read in libraries
library(deSolve)
library(here)

#Create model
dilution.model <- function (t, x, params) {
  
  #State variables
  S1 = x[1]
  I1 = x[2]
  S2 = x[3]
  I2 = x[4]
  P = x[5]
  
  #Paramters
  r1 = params[1]
  r2 = params[2]
  intra = params[3]
  inter1on2 = params[4]
  inter2on1 = params[5]
  feed1 = parameters[6]
  feed2 = parameters[7]
  inf1 = params[8]
  inf2 = params[9]
  K1 = params[10]
  K2 = params[11]
  mort1 = params[12]
  mort2 = params[13]
  release1 = params[14]
  release2 = params[15]
  degrade = params[16]
  
  #Functions as described in manuscript (equations 3-7)
  dS1dt = (S1 + I1) * (r1 - intra * (S1+I1) - inter2on1 * (S2+I2)) - feed1 * inf1 * (P ^ K1) * S1
  dI1dt = feed1 * inf1 * (P ^ K1) * S1 - mort1 * I1
  dS2dt = (S2 + I2) * (r2 - intra * (S2+I2) - inter1on2 * (S1+I1)) - feed2 * inf2 * (P ^ K2) * S2
  dI2dt = feed2 * inf2 * (P ^ K2) * S2 - mort2 * I2
  dPdt = release1 * I1 + release2 * I2 - degrade * P - feed1 * (S1 + I1) * P - feed2 * (S2 + I2) * P
  
  if((S1 + dS1dt) < 0){dS1dt <- -S1}
  if((I1 + dI1dt) < 0){dI1dt <- -I1}
  if((S2 + dS2dt) < 0){dS2dt <- -S2}
  if((I2 + dI2dt) < 0){dI2dt <- -I2}
  if((P + dPdt) < 0){dPdt <- -P}
  
  return(list(c(dS1dt,dI1dt,dS2dt,dI2dt,dPdt)))
}

############## RUN MODEL FOR PARAMETERS IN MAIN TEXT ##################

#Create index for use in recording results
index = 1

#Create matrix for recording results
dilution_curve <- matrix(0,(400*4*3),9) #(400 r1, 4 release, 3 K)

#Runs for 1000 time steps (enough to reach equilibria)
times = seq(0, 1000, by=0.01)

#Loop over dose-infectivity values
for(ex in 1:3){
  
  #Loop over propagule excretion values of competitor
  for(Z in 0:3){
    
    release2 = Z * 0.5 * 100 #0.75 is 
    
    #Loop over growth rates of competing host
    for(i in 1:400){
      
      ##Initial values of state variables
      S10 = 1
      I10 = 1 
      S20 = 1 
      I20 = 1
      P0 = 50 
      
      initial_values = c(S10,I10,S20,I20,P0)
      
      #Parameter values
      r1 = 1
      r2 = (1/200) * i
      intra = 1 
      inter1on2 = 0.5
      inter2on1 = 0.5
      feed1 = 1
      feed2 = 1
      K1 = ex/2
      K2 = ex/2
      mort1 = 0.1
      mort2 = 0.1
      release1 = 100
      degrade = .1
      
      inf1 = mort1/(feed1*(release1*((2*r1-mort2)/4*intra)/(degrade+feed1*((2*r1-mort1)/(2*intra))))^K1)
      inf2 = inf1
      
      
      parameters = c(r1,r2,intra,inter1on2,inter2on1,feed1,feed2,inf1,inf2,K1,K2,mort1,mort2,release1,release2,degrade)
      
      #Run model
      results = lsoda(initial_values, times, dilution.model, parameters)
      
      #Record Results
      dilution_curve[index,1]<-r1
      dilution_curve[index,2]<-inter2on1
      dilution_curve[index,3]<-release1
      dilution_curve[index,4]<-K1
      dilution_curve[index,5:9]<-results[length(times),2:6]
      
      index = index + 1
    }
    
  }  
  
}

# Naming columns for easy identification
colnames(dilution_curve) = c("r1","inter","dilution", "K", "S1", "I1","S2","I2","P")

#Save data
dilution_curve <- as.data.frame(dilution_curve)
write.csv(dilution_curve,"base.csv")


############## RUN MODEL when non-focal host is superior competitor ##################

#Create index for use in recording results
index = 1

#Create matrix for recording results
dilution_curve <- matrix(0,(400*4*3),9) #(400 r1, 4 release, 3 K)

#Runs for 1000 time steps (enough to reach equilibria)
times = seq(0, 1000, by=0.01)

#Loop over dose-infectivity values
for(ex in 1:3){
  
  #Loop over propagule excretion values of competitor
  for(Z in 0:3){
    
    release2 = Z * 0.5 * 100 #0.75 is 
    
    #Loop over growth rates of competing host
    for(i in 1:400){
      
      ##Initial values of state variables
      S10 = 1
      I10 = 1 
      S20 = 1 
      I20 = 1
      P0 = 50 
      
      initial_values = c(S10,I10,S20,I20,P0)
      
      #Parameter values
      r1 = 1
      r2 = (1/200) * i
      intra = 1 
      inter1on2 = 0.5
      inter2on1 = 0.75
      feed1 = 1
      feed2 = 1
      K1 = ex/2
      K2 = ex/2
      mort1 = 0.1
      mort2 = 0.1
      release1 = 100
      degrade = .1
      
      inf1 = mort1/(feed1*(release1*((2*r1-mort2)/4*intra)/(degrade+feed1*((2*r1-mort1)/(2*intra))))^K1)
      inf2 = inf1
      
      
      parameters = c(r1,r2,intra,inter1on2,inter2on1,feed1,feed2,inf1,inf2,K1,K2,mort1,mort2,release1,release2,degrade)
      
      #Run model
      results = lsoda(initial_values, times, dilution.model, parameters)
      
      #Record Results
      dilution_curve[index,1]<-r1
      dilution_curve[index,2]<-inter2on1
      dilution_curve[index,3]<-release1
      dilution_curve[index,4]<-K1
      dilution_curve[index,5:9]<-results[length(times),2:6]
      
      index = index + 1
    }
    
  }  
  
}

# Naming columns for easy identification
colnames(dilution_curve) = c("r1","inter","dilution", "K", "S1", "I1","S2","I2","P")

#Save data
dilution_curve <- as.data.frame(dilution_curve)
write.csv(dilution_curve,"base_high_comp.csv")

############## RUN MODEL when non-focal host weak competitor ##################

#Create index for use in recording results
index = 1

#Create matrix for recording results
dilution_curve <- matrix(0,(400*4*3),9) #(400 r1, 4 release, 3 K)

#Runs for 1000 time steps (enough to reach equilibria)
times = seq(0, 1000, by=0.01)

#Loop over dose-infectivity values
for(ex in 1:3){
  
  #Loop over propagule excretion values of competitor
  for(Z in 0:3){
    
    release2 = Z * 0.5 * 100 #0.75 is 
    
    #Loop over growth rates of competing host
    for(i in 1:400){
      
      ##Initial values of state variables
      S10 = 1
      I10 = 1 
      S20 = 1 
      I20 = 1
      P0 = 50 
      
      initial_values = c(S10,I10,S20,I20,P0)
      
      #Parameter values
      r1 = 1
      r2 = (1/200) * i
      intra = 1 
      inter1on2 = 0.5
      inter2on1 = 0.25
      feed1 = 1
      feed2 = 1
      K1 = ex/2
      K2 = ex/2
      mort1 = 0.1
      mort2 = 0.1
      release1 = 100
      degrade = .1
      
      inf1 = mort1/(feed1*(release1*((2*r1-mort2)/4*intra)/(degrade+feed1*((2*r1-mort1)/(2*intra))))^K1)
      inf2 = inf1
      
      
      parameters = c(r1,r2,intra,inter1on2,inter2on1,feed1,feed2,inf1,inf2,K1,K2,mort1,mort2,release1,release2,degrade)
      
      #Run model
      results = lsoda(initial_values, times, dilution.model, parameters)
      
      #Record Results
      dilution_curve[index,1]<-r1
      dilution_curve[index,2]<-inter2on1
      dilution_curve[index,3]<-release1
      dilution_curve[index,4]<-K1
      dilution_curve[index,5:9]<-results[length(times),2:6]
      
      index = index + 1
    }
    
  }  
  
}

# Naming columns for easy identification
colnames(dilution_curve) = c("r1","inter","dilution", "K", "S1", "I1","S2","I2","P")

#Save data
dilution_curve <- as.data.frame(dilution_curve)
write.csv(dilution_curve,"base_low_comp.csv")

############## RUN MODEL when non-focal host has increased contact rate ##################

#Create index for use in recording results
index = 1

#Create matrix for recording results
dilution_curve <- matrix(0,(400*4*3),9) #(400 r1, 4 release, 3 K)

#Runs for 1000 time steps (enough to reach equilibria)
times = seq(0, 1000, by=0.01)

#Loop over dose-infectivity values
for(ex in 1:3){
  
  #Loop over propagule excretion values of competitor
  for(Z in 0:3){
    
    release2 = Z * 0.5 * 100 #0.75 is 
    
    #Loop over growth rates of competing host
    for(i in 1:400){
      
      ##Initial values of state variables
      S10 = 1
      I10 = 1 
      S20 = 1 
      I20 = 1
      P0 = 50 
      
      initial_values = c(S10,I10,S20,I20,P0)
      
      #Parameter values
      r1 = 1
      r2 = (1/200) * i
      intra = 1 
      inter1on2 = 0.5
      inter2on1 = 0.5
      feed1 = 1
      feed2 = 1.5
      K1 = ex/2
      K2 = ex/2
      mort1 = 0.1
      mort2 = 0.1
      release1 = 100
      degrade = .1
      
      inf1 = mort1/(feed1*(release1*((2*r1-mort2)/4*intra)/(degrade+feed1*((2*r1-mort1)/(2*intra))))^K1)
      inf2 = inf1
      
      
      parameters = c(r1,r2,intra,inter1on2,inter2on1,feed1,feed2,inf1,inf2,K1,K2,mort1,mort2,release1,release2,degrade)
      
      #Run model
      results = lsoda(initial_values, times, dilution.model, parameters)
      
      #Record Results
      dilution_curve[index,1]<-r1
      dilution_curve[index,2]<-inter2on1
      dilution_curve[index,3]<-release1
      dilution_curve[index,4]<-K1
      dilution_curve[index,5:9]<-results[length(times),2:6]
      
      index = index + 1
    }
    
  }  
  
}

# Naming columns for easy identification
colnames(dilution_curve) = c("r1","inter","dilution", "K", "S1", "I1","S2","I2","P")

#Save data
dilution_curve <- as.data.frame(dilution_curve)
write.csv(dilution_curve,"base_high_f.csv")

############## RUN MODEL when non-focal host has lower contact rate ##################

#Create index for use in recording results
index = 1

#Create matrix for recording results
dilution_curve <- matrix(0,(400*4*3),9) #(400 r1, 4 release, 3 K)

#Runs for 1000 time steps (enough to reach equilibria)
times = seq(0, 1000, by=0.01)

#Loop over dose-infectivity values
for(ex in 1:3){
  
  #Loop over propagule excretion values of competitor
  for(Z in 0:3){
    
    release2 = Z * 0.5 * 100 #0.75 is 
    
    #Loop over growth rates of competing host
    for(i in 1:400){
      
      ##Initial values of state variables
      S10 = 1
      I10 = 1 
      S20 = 1 
      I20 = 1
      P0 = 50 
      
      initial_values = c(S10,I10,S20,I20,P0)
      
      #Parameter values
      r1 = 1
      r2 = (1/200) * i
      intra = 1 
      inter1on2 = 0.5
      inter2on1 = 0.5
      feed1 = 1
      feed2 = 0.5
      K1 = ex/2
      K2 = ex/2
      mort1 = 0.1
      mort2 = 0.1
      release1 = 100
      degrade = .1
      
      inf1 = mort1/(feed1*(release1*((2*r1-mort2)/4*intra)/(degrade+feed1*((2*r1-mort1)/(2*intra))))^K1)
      inf2 = inf1
      
      
      parameters = c(r1,r2,intra,inter1on2,inter2on1,feed1,feed2,inf1,inf2,K1,K2,mort1,mort2,release1,release2,degrade)
      
      #Run model
      results = lsoda(initial_values, times, dilution.model, parameters)
      
      #Record Results
      dilution_curve[index,1]<-r1
      dilution_curve[index,2]<-inter2on1
      dilution_curve[index,3]<-release1
      dilution_curve[index,4]<-K1
      dilution_curve[index,5:9]<-results[length(times),2:6]
      
      index = index + 1
    }
    
  }  
  
}

# Naming columns for easy identification
colnames(dilution_curve) = c("r1","inter","dilution", "K", "S1", "I1","S2","I2","P")

#Save data
dilution_curve <- as.data.frame(dilution_curve)
write.csv(dilution_curve,"base.csv")

############## RUN MODEL when non-focal host has high susceptibility ##################

#Create index for use in recording results
index = 1

#Create matrix for recording results
dilution_curve <- matrix(0,(400*4*3),9) #(400 r1, 4 release, 3 K)

#Runs for 1000 time steps (enough to reach equilibria)
times = seq(0, 1000, by=0.01)

#Loop over dose-infectivity values
for(ex in 1:3){
  
  #Loop over propagule excretion values of competitor
  for(Z in 0:3){
    
    release2 = Z * 0.5 * 100 #0.75 is 
    
    #Loop over growth rates of competing host
    for(i in 1:400){
      
      ##Initial values of state variables
      S10 = 1
      I10 = 1 
      S20 = 1 
      I20 = 1
      P0 = 50 
      
      initial_values = c(S10,I10,S20,I20,P0)
      
      #Parameter values
      r1 = 1
      r2 = (1/200) * i
      intra = 1 
      inter1on2 = 0.5
      inter2on1 = 0.5
      feed1 = 1
      feed2 = 1
      K1 = ex/2
      K2 = ex/2
      mort1 = 0.1
      mort2 = 0.1
      release1 = 100
      degrade = .1
      
      inf1 = mort1/(feed1*(release1*((2*r1-mort2)/4*intra)/(degrade+feed1*((2*r1-mort1)/(2*intra))))^K1)
      inf2 = inf1*1.5
      
      
      parameters = c(r1,r2,intra,inter1on2,inter2on1,feed1,feed2,inf1,inf2,K1,K2,mort1,mort2,release1,release2,degrade)
      
      #Run model
      results = lsoda(initial_values, times, dilution.model, parameters)
      
      #Record Results
      dilution_curve[index,1]<-r1
      dilution_curve[index,2]<-inter2on1
      dilution_curve[index,3]<-release1
      dilution_curve[index,4]<-K1
      dilution_curve[index,5:9]<-results[length(times),2:6]
      
      index = index + 1
    }
    
  }  
  
}

# Naming columns for easy identification
colnames(dilution_curve) = c("r1","inter","dilution", "K", "S1", "I1","S2","I2","P")

#Save data
dilution_curve <- as.data.frame(dilution_curve)
write.csv(dilution_curve,"base_high_sus.csv")

############## RUN MODEL when non-focal host has low susceptibility ##################

#Create index for use in recording results
index = 1

#Create matrix for recording results
dilution_curve <- matrix(0,(400*4*3),9) #(400 r1, 4 release, 3 K)

#Runs for 1000 time steps (enough to reach equilibria)
times = seq(0, 1000, by=0.01)

#Loop over dose-infectivity values
for(ex in 1:3){
  
  #Loop over propagule excretion values of competitor
  for(Z in 0:3){
    
    release2 = Z * 0.5 * 100 #0.75 is 
    
    #Loop over growth rates of competing host
    for(i in 1:400){
      
      ##Initial values of state variables
      S10 = 1
      I10 = 1 
      S20 = 1 
      I20 = 1
      P0 = 50 
      
      initial_values = c(S10,I10,S20,I20,P0)
      
      #Parameter values
      r1 = 1
      r2 = (1/200) * i
      intra = 1 
      inter1on2 = 0.5
      inter2on1 = 0.5
      feed1 = 1
      feed2 = 1
      K1 = ex/2
      K2 = ex/2
      mort1 = 0.1
      mort2 = 0.1
      release1 = 100
      degrade = .1
      
      inf1 = mort1/(feed1*(release1*((2*r1-mort2)/4*intra)/(degrade+feed1*((2*r1-mort1)/(2*intra))))^K1)
      inf2 = inf1 * 0.5
      
      
      parameters = c(r1,r2,intra,inter1on2,inter2on1,feed1,feed2,inf1,inf2,K1,K2,mort1,mort2,release1,release2,degrade)
      
      #Run model
      results = lsoda(initial_values, times, dilution.model, parameters)
      
      #Record Results
      dilution_curve[index,1]<-r1
      dilution_curve[index,2]<-inter2on1
      dilution_curve[index,3]<-release1
      dilution_curve[index,4]<-K1
      dilution_curve[index,5:9]<-results[length(times),2:6]
      
      index = index + 1
    }
    
  }  
  
}

# Naming columns for easy identification
colnames(dilution_curve) = c("r1","inter","dilution", "K", "S1", "I1","S2","I2","P")

#Save data
dilution_curve <- as.data.frame(dilution_curve)
write.csv(dilution_curve,"base_low_sus.csv")

############## RUN MODEL when non-focal host has high mortality ##################

#Create index for use in recording results
index = 1

#Create matrix for recording results
dilution_curve <- matrix(0,(400*4*3),9) #(400 r1, 4 release, 3 K)

#Runs for 1000 time steps (enough to reach equilibria)
times = seq(0, 1000, by=0.01)

#Loop over dose-infectivity values
for(ex in 1:3){
  
  #Loop over propagule excretion values of competitor
  for(Z in 0:3){
    
    release2 = Z * 0.5 * 100 #0.75 is 
    
    #Loop over growth rates of competing host
    for(i in 1:400){
      
      ##Initial values of state variables
      S10 = 1
      I10 = 1 
      S20 = 1 
      I20 = 1
      P0 = 50 
      
      initial_values = c(S10,I10,S20,I20,P0)
      
      #Parameter values
      r1 = 1
      r2 = (1/200) * i
      intra = 1 
      inter1on2 = 0.5
      inter2on1 = 0.5
      feed1 = 1
      feed2 = 1
      K1 = ex/2
      K2 = ex/2
      mort1 = 0.1
      mort2 = 0.15
      release1 = 100
      degrade = .1
      
      inf1 = mort1/(feed1*(release1*((2*r1-mort2)/4*intra)/(degrade+feed1*((2*r1-mort1)/(2*intra))))^K1)
      inf2 = inf1
      
      
      parameters = c(r1,r2,intra,inter1on2,inter2on1,feed1,feed2,inf1,inf2,K1,K2,mort1,mort2,release1,release2,degrade)
      
      #Run model
      results = lsoda(initial_values, times, dilution.model, parameters)
      
      #Record Results
      dilution_curve[index,1]<-r1
      dilution_curve[index,2]<-inter2on1
      dilution_curve[index,3]<-release1
      dilution_curve[index,4]<-K1
      dilution_curve[index,5:9]<-results[length(times),2:6]
      
      index = index + 1
    }
    
  }  
  
}

# Naming columns for easy identification
colnames(dilution_curve) = c("r1","inter","dilution", "K", "S1", "I1","S2","I2","P")

#Save data
dilution_curve <- as.data.frame(dilution_curve)
write.csv(dilution_curve,"base_high_mort.csv")

############## RUN MODEL when model has low mortality ##################

#Create index for use in recording results
index = 1

#Create matrix for recording results
dilution_curve <- matrix(0,(400*4*3),9) #(400 r1, 4 release, 3 K)

#Runs for 1000 time steps (enough to reach equilibria)
times = seq(0, 1000, by=0.01)

#Loop over dose-infectivity values
for(ex in 1:3){
  
  #Loop over propagule excretion values of competitor
  for(Z in 0:3){
    
    release2 = Z * 0.5 * 100 #0.75 is 
    
    #Loop over growth rates of competing host
    for(i in 1:400){
      
      ##Initial values of state variables
      S10 = 1
      I10 = 1 
      S20 = 1 
      I20 = 1
      P0 = 50 
      
      initial_values = c(S10,I10,S20,I20,P0)
      
      #Parameter values
      r1 = 1
      r2 = (1/200) * i
      intra = 1 
      inter1on2 = 0.5
      inter2on1 = 0.5
      feed1 = 1
      feed2 = 1
      K1 = ex/2
      K2 = ex/2
      mort1 = 0.1
      mort2 = 0.05
      release1 = 100
      degrade = .1
      
      inf1 = mort1/(feed1*(release1*((2*r1-mort2)/4*intra)/(degrade+feed1*((2*r1-mort1)/(2*intra))))^K1)
      inf2 = inf1
      
      
      parameters = c(r1,r2,intra,inter1on2,inter2on1,feed1,feed2,inf1,inf2,K1,K2,mort1,mort2,release1,release2,degrade)
      
      #Run model
      results = lsoda(initial_values, times, dilution.model, parameters)
      
      #Record Results
      dilution_curve[index,1]<-r1
      dilution_curve[index,2]<-inter2on1
      dilution_curve[index,3]<-release1
      dilution_curve[index,4]<-K1
      dilution_curve[index,5:9]<-results[length(times),2:6]
      
      index = index + 1
    }
    
  }  
  
}

# Naming columns for easy identification
colnames(dilution_curve) = c("r1","inter","dilution", "K", "S1", "I1","S2","I2","P")

#Save data
dilution_curve <- as.data.frame(dilution_curve)
write.csv(dilution_curve,"base_low_mort.csv")

############## RUN MODEL when dose-infectivity relationships of non-focal 
############## host is always decelerating ##################

#Create index for use in recording results
index = 1

#Create matrix for recording results
dilution_curve <- matrix(0,(400*4*3),9) #(400 r1, 4 release, 3 K)

#Runs for 1000 time steps (enough to reach equilibria)
times = seq(0, 1000, by=0.01)

#Loop over dose-infectivity values
for(ex in 1:3){
  
  #Loop over propagule excretion values of competitor
  for(Z in 0:3){
    
    release2 = Z * 0.5 * 100 #0.75 is 
    
    #Loop over growth rates of competing host
    for(i in 1:400){
      
      ##Initial values of state variables
      S10 = 1
      I10 = 1 
      S20 = 1 
      I20 = 1
      P0 = 50 
      
      initial_values = c(S10,I10,S20,I20,P0)
      
      #Parameter values
      r1 = 1
      r2 = (1/200) * i
      intra = 1 
      inter1on2 = 0.5
      inter2on1 = 0.5
      feed1 = 1
      feed2 = 1
      K1 = ex/2
      K2 = 0.5
      mort1 = 0.1
      mort2 = 0.1
      release1 = 100
      degrade = .1
      
      inf1 = mort1/(feed1*(release1*((2*r1-mort2)/4*intra)/(degrade+feed1*((2*r1-mort1)/(2*intra))))^K1)
      inf2 = inf1
      
      
      parameters = c(r1,r2,intra,inter1on2,inter2on1,feed1,feed2,inf1,inf2,K1,K2,mort1,mort2,release1,release2,degrade)
      
      #Run model
      results = lsoda(initial_values, times, dilution.model, parameters)
      
      #Record Results
      dilution_curve[index,1]<-r1
      dilution_curve[index,2]<-inter2on1
      dilution_curve[index,3]<-release1
      dilution_curve[index,4]<-K1
      dilution_curve[index,5:9]<-results[length(times),2:6]
      
      index = index + 1
    }
    
  }  
  
}

# Naming columns for easy identification
colnames(dilution_curve) = c("r1","inter","dilution", "K", "S1", "I1","S2","I2","P")

#Save data
dilution_curve <- as.data.frame(dilution_curve)
write.csv(dilution_curve,"base_k_dec.csv")

############## RUN MODEL when dose-infectivity relationships of non-focal 
############## host is always linear ##################

#Create index for use in recording results
index = 1

#Create matrix for recording results
dilution_curve <- matrix(0,(400*4*3),9) #(400 r1, 4 release, 3 K)

#Runs for 1000 time steps (enough to reach equilibria)
times = seq(0, 1000, by=0.01)

#Loop over dose-infectivity values
for(ex in 1:3){
  
  #Loop over propagule excretion values of competitor
  for(Z in 0:3){
    
    release2 = Z * 0.5 * 100 #0.75 is 
    
    #Loop over growth rates of competing host
    for(i in 1:400){
      
      ##Initial values of state variables
      S10 = 1
      I10 = 1 
      S20 = 1 
      I20 = 1
      P0 = 50 
      
      initial_values = c(S10,I10,S20,I20,P0)
      
      #Parameter values
      r1 = 1
      r2 = (1/200) * i
      intra = 1 
      inter1on2 = 0.5
      inter2on1 = 0.5
      feed1 = 1
      feed2 = 1
      K1 = ex/2
      K2 = 1.0
      mort1 = 0.1
      mort2 = 0.1
      release1 = 100
      degrade = .1
      
      inf1 = mort1/(feed1*(release1*((2*r1-mort2)/4*intra)/(degrade+feed1*((2*r1-mort1)/(2*intra))))^K1)
      inf2 = inf1
      
      
      parameters = c(r1,r2,intra,inter1on2,inter2on1,feed1,feed2,inf1,inf2,K1,K2,mort1,mort2,release1,release2,degrade)
      
      #Run model
      results = lsoda(initial_values, times, dilution.model, parameters)
      
      #Record Results
      dilution_curve[index,1]<-r1
      dilution_curve[index,2]<-inter2on1
      dilution_curve[index,3]<-release1
      dilution_curve[index,4]<-K1
      dilution_curve[index,5:9]<-results[length(times),2:6]
      
      index = index + 1
    }
    
  }  
  
}

# Naming columns for easy identification
colnames(dilution_curve) = c("r1","inter","dilution", "K", "S1", "I1","S2","I2","P")

#Save data
dilution_curve <- as.data.frame(dilution_curve)
write.csv(dilution_curve,"base_k_lin.csv")

############## RUN MODEL when dose-infectivity relationships of non-focal 
############## host is always accelerating ##################

#Create index for use in recording results
index = 1

#Create matrix for recording results
dilution_curve <- matrix(0,(400*4*3),9) #(400 r1, 4 release, 3 K)

#Runs for 1000 time steps (enough to reach equilibria)
times = seq(0, 1000, by=0.01)

#Loop over dose-infectivity values
for(ex in 1:3){
  
  #Loop over propagule excretion values of competitor
  for(Z in 0:3){
    
    release2 = Z * 0.5 * 100 #0.75 is 
    
    #Loop over growth rates of competing host
    for(i in 1:400){
      
      ##Initial values of state variables
      S10 = 1
      I10 = 1 
      S20 = 1 
      I20 = 1
      P0 = 50 
      
      initial_values = c(S10,I10,S20,I20,P0)
      
      #Parameter values
      r1 = 1
      r2 = (1/200) * i
      intra = 1 
      inter1on2 = 0.5
      inter2on1 = 0.5
      feed1 = 1
      feed2 = 1
      K1 = ex/2
      K2 = 1.5
      mort1 = 0.1
      mort2 = 0.1
      release1 = 100
      degrade = .1
      
      inf1 = mort1/(feed1*(release1*((2*r1-mort2)/4*intra)/(degrade+feed1*((2*r1-mort1)/(2*intra))))^K1)
      inf2 = inf1
      
      
      parameters = c(r1,r2,intra,inter1on2,inter2on1,feed1,feed2,inf1,inf2,K1,K2,mort1,mort2,release1,release2,degrade)
      
      #Run model
      results = lsoda(initial_values, times, dilution.model, parameters)
      
      #Record Results
      dilution_curve[index,1]<-r1
      dilution_curve[index,2]<-inter2on1
      dilution_curve[index,3]<-release1
      dilution_curve[index,4]<-K1
      dilution_curve[index,5:9]<-results[length(times),2:6]
      
      index = index + 1
    }
    
  }  
  
}

# Naming columns for easy identification
colnames(dilution_curve) = c("r1","inter","dilution", "K", "S1", "I1","S2","I2","P")

#Save data
dilution_curve <- as.data.frame(dilution_curve)
write.csv(dilution_curve,"base_k_acc.csv")

