
Cerrado_FlorestaG_intHT_long<-nimbleCode({

  ## Community-level priors
  
  # Hyperpriors for occupancy model
  
  mu.lpsi ~ dnorm(0, sd=2.5)
  sd.mu.lpsi ~ dunif(0,5)
  
  ##covariate - different effects for species groups, shared SD
  for (s in 1:3){
    mu.betalpsi1[s] ~ dnorm(0, sd=2.5)
    mu.betalpsi2[s] ~ dnorm(0, sd=2.5)
  }
    ##not for effect of sav cover bc data too sparse
    for(h in 1:3){
	mu.betalpsi3[h] ~ dnorm(0, sd=2.5)
    }
    sd.betalpsi1 ~ T(dnorm(0, sd=2), 0.0001,)
    sd.betalpsi2 ~ dunif(0,5)
    sd.betalpsi3 ~ dunif(0,5)
    

  ##prior for random effect precision
  sig.ran~dunif(0,5)
  
  # Hyperpriors for detection model
  
  mu.lp ~ dnorm(0, sd=2.5)
  sd.mu.lp ~ dunif(0,5)
  
  mu.betalp2 ~ dnorm(0, sd=2.5)
  sd.betalp2 ~ dunif(0,5)
    
  ## species-site-level parameters (all from a single random effect): intercept for occu and detection
  ## effect of covariate on occupancy depends on geographic region (shared SD)
  
  for(k in 1:n.spec) {      ##loop over all species (collectively for all sites)           
    
    for (s in 1:n.survey){     ##one for each survey
      
      lpsi[k,s] ~ dnorm(mu.lpsi, sd=sd.mu.lpsi) 
      lp[k,s] ~ dnorm(mu.lp, sd=sd.mu.lp)
      
    }
    
    #efeito tipo de habitat (categorico)
      betalpsi1[k] ~ dnorm(mu.betalpsi1[group[k]], sd=sd.betalpsi1) #efeito habitat - dummy 1 (mata ciliar)
      betalpsi2[k] ~ dnorm(mu.betalpsi2[group[k]], sd=sd.betalpsi2) #efeito habitat - dummy 2 (savana)

      ##efeito savana scaled por habitat
      for (h in 1:3){
             betalpsi3[k,h] ~ dnorm(mu.betalpsi3[h], sd=sd.betalpsi3)  #faria com shared SD entre os tipos de habitat
      }
    
    betalp2[k] ~ dnorm(mu.betalp2, sd=sd.betalp2)
    #effect habitat type on detection
  }
  
  ## Effect of effort, fixed across species, surveys
  betalp1 ~ dnorm(0, sd=2.5)				


  
  for(k in 1:n.sp.surv){              ##loop over all species stacked for all 5 surveys; vectorized over sites
    
    logit(psi[1:nsite[survey[k]], k]) <- lpsi[species[k], survey[k]] +  			#intercept
              betalpsi1[species[k]]*COV[survey[k],1:nsite[survey[k]],8] +  	 		#habitat dummy 1
              betalpsi2[species[k]]*COV[survey[k],1:nsite[survey[k]],9] +    			#habitat dummy 2
              betalpsi3[species[k],1]*COV[survey[k],1:nsite[survey[k]],6] *  			#interacao, dummy 1
                                                        COV[survey[k],1:nsite[survey[k]],8] + 
              betalpsi3[species[k],2]*COV[survey[k],1:nsite[survey[k]],6] *   			#interacao, dummy 1
                                                        COV[survey[k],1:nsite[survey[k]],9] +
              betalpsi3[species[k],3]*COV[survey[k],1:nsite[survey[k]],6] *   			#interacao, reference cat
                                                        COV[survey[k],1:nsite[survey[k]],7]  
       
    
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


