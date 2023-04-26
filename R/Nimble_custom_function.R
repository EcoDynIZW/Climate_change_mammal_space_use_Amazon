### custom density function to vectorize observations and integrate out latent z

dBernM2 <- nimbleFunction( run = function(x = double(1), ##observation
                                         detectionProb=double(1),   ##p
                                         occProb = double(0), ##psi
                                         log = double(0)) {
  
  returnType(double())

  nocc<-length(x) ##number of occasions
  
  probDetectionHistoryGivenOccupied <- 1
  probDetectionHistoryGivenUnoccupied <- 1
  
  for (k in 1:nocc){
    
  if(x[k] == 1) {
    probDetectionHistoryGivenOccupied <-
      probDetectionHistoryGivenOccupied * detectionProb[k]
    probDetectionHistoryGivenUnoccupied <- 0
  } else {
    probDetectionHistoryGivenOccupied <-
      probDetectionHistoryGivenOccupied * (1-detectionProb[k])
  }
    
}
  ans <- log(occProb *
               probDetectionHistoryGivenOccupied +
               (1-occProb) *
               probDetectionHistoryGivenUnoccupied)
  
  if(log) return(ans)
  return(exp(ans))
  
})


rBernM2 <- nimbleFunction( run = function(n = integer(0), ##default
                                          detectionProb=double(1),   ##p
                                          occProb = double(0)) {
  
  returnType(double(1))
  
  nocc<-length(detectionProb) ##number of occasions
  z<-rbinom(1, 1, occProb)
  x<-rbinom(nocc, 1, detectionProb*z)
  return(x)
  
})



registerDistributions(list(dBernM2 = list(
  BUGSdist = "dBernM2(detectionProb, occProb)",
  Rdist = "dBernM2(detectionProb,occProb)",
  types = c('value = double(1)',
            'detectionProb = double(1)',
            'occProb = double(0)'
  )
))
)




