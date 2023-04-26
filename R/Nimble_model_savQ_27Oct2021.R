
Cerrado_FlorestaQ_long<-nimbleCode({

  ## Community-level priors
  
  # Hyperpriors for occupancy model
  
  mu.lpsi ~ dnorm(0, sd=2.5)
  sd.mu.lpsi ~ dunif(0,5)
  
  ##covariate - different effects in two geographic regions, shared SD
  for (s in 1:2){   
    mu.betalpsi[s] ~ dnorm(0, sd=2.5)
    mu.betalpsi2[s] ~ dnorm(0, sd=2.5)
  }
  
  sd.betalpsi ~ dunif(0,5)
  sd.betalpsi2~T(dnorm(0, sd=1.5), 0.0001,)

  ##prior for random effect precision
  # sig.ran~dunif(0,5)
  # try informative prior to limit sig.ran
  sig.ran~dgamma(0.5,2.5)
  
  # Hyperpriors for detection model
  
  mu.lp ~ dnorm(0, sd=2.5)
  sd.mu.lp ~ dunif(0,5)
  
  mu.betalp2 ~ dnorm(0, sd=2.5)
  sd.betalp2 ~ dunif(0,5)
  # 
  # mu.betalp3 ~ dnorm(0, sd=2.5)
  # sd.betalp3 ~ dunif(0,5)
  
  ## species-site-level parameters (all from a single random effect): intercept for occu and detection
  ## effect of covariate on occupancy depends on geographic region (shared SD)
  
  for(k in 1:n.spec) {      ##loop over all species (collectively for all sites)           
    
    for (s in 1:n.survey){     ##one for each survey
      
      lpsi[k,s] ~ dnorm(mu.lpsi, sd=sd.mu.lpsi) 
      lp[k,s] ~ dnorm(mu.lp, sd=sd.mu.lp)
      
    }
    
    betalpsi[k,1] ~ dnorm(mu.betalpsi[1], sd=sd.betalpsi) 
    betalpsi[k,2] ~ dnorm(mu.betalpsi[2], sd=sd.betalpsi) #allows for different effects in norther vs southern areas
    betalpsi2[k,1] ~ dnorm(mu.betalpsi2[1], sd=sd.betalpsi2) 
    betalpsi2[k,2] ~ dnorm(mu.betalpsi2[2], sd=sd.betalpsi2) #allows for different effects in norther vs southern areas
    
    betalp2[k] ~ dnorm(mu.betalp2, sd=sd.betalp2)
    #effect habitat type on detection
  }
  
  ## Effect of effort, fixed across species, surveys
  betalp1 ~ dnorm(0, sd=2.5)				


  
  for(k in 1:n.sp.surv){              ##loop over all species stacked for all 5 surveys; vectorized over sites
    
    logit(psi[1:nsite[survey[k]], k]) <- lpsi[species[k], survey[k]] + 
              betalpsi[species[k], region[k]]*COV[survey[k],1:nsite[survey[k]],12] +
       	  betalpsi2[species[k], region[k]]*(COV[survey[k],1:nsite[survey[k]],12]^2)
    
    for (i in 1:nsite[survey[k]]) { ##only loop over number of locations sampled in each survey
      
      eps[k,i]~dnorm(0,sd=sig.ran) #site level random effect on p    
      
      logit(p[i,first[survey[k],i]:last[survey[k],i],k]) <- lp[species[k], survey[k]] + 
                                                            betalp1*EFF[survey[k],i,first[survey[k],i]:last[survey[k],i]]  +
                                                            betalp2[species[k]]*COV[survey[k],i,2] + 
                                                            eps[k,i]
                                                                  

        y[i,first[survey[k],i]:last[survey[k],i],k] ~ dBernM2(detectionProb=p[i,first[survey[k],i]:last[survey[k],i],k], 
                                                              occProb=psi[i,k])
        
        ### for Bp value
        y.new[i,first[survey[k],i]:last[survey[k],i],k] ~ dBernM2(detectionProb=p[i,first[survey[k],i]:last[survey[k],i],k], 
                                                              occProb=psi[i,k])
        
        FT[k,i]<-(sqrt(sum(  y[i,first[survey[k],i]:last[survey[k],i],k]  ))-sqrt(sum( p[i,first[survey[k],i]:last[survey[k],i],k]*psi[i,k] )))^2
        FT.new[k,i]<-(sqrt(sum(y.new[i,first[survey[k],i]:last[survey[k],i],k]))-sqrt(sum( p[i,first[survey[k],i]:last[survey[k],i],k]*psi[i,k] )))^2

  } #end site loop
    
    ## for species-survey specific Bp value 
    FT2[k]<-sum(FT[k,1:nsite[survey[k]]]) ##sum only over non-NA values to avoid problems
    FT2.new[k]<-sum(FT.new[k,1:nsite[survey[k]]]) ##sum only over non-NA values to avoid problems
    bp[k]<-FT2[k]>FT2.new[k]
} #end outer loop
  
  ## for overall Bayeian p-value
  BP<-sum(FT2[1:n.sp.surv])>sum(FT2.new[1:n.sp.surv])
  
})


