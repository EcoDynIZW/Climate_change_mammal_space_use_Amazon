##############################################################################
## loads data and runs model described first in main text, with effect of 
## savannah at 100-m scaled, no group-specific effects

library(nimble)
library(MCMCvis)

##load all data
all.data<-dget("data-raw/Rocha_Sollmann_data.R")

## source nimble model code, custom functions
source("R/Nimble_model_savQ_27Oct2021.R")
source("R/Nimble_custom_function.R") 

## create initial values
inits<-function(){list(mu.lpsi=rnorm(1), sd.mu.lpsi=runif(1, 0.5, 2), 
                       mu.betalpsi=rnorm(2,0,1), sd.betalpsi=runif(1, 0.5, 2),
			     mu.betalpsi2=rnorm(2,0,1), sd.betalpsi2=runif(1, 0.5, 2),
                       mu.betalp2=rnorm(1,0,1), sd.betalp2=runif(1, 0.5, 2),
                       mu.lp=rnorm(1), sd.mu.lp=runif(1, 0.5, 2),
                       betalp1=rnorm(1),
			     sig.ran=runif(1, 0.1, 0.5))}

#set parameters to monitor
params<-c("mu.lpsi","sd.mu.lpsi","mu.betalpsi", "sd.betalpsi", "mu.betalpsi2", "sd.betalpsi2",  
          "mu.lp", "sd.mu.lp",'betalp1',"sig.ran",                        
          "mu.betalp2", "sd.betalp2",
          "lpsi", "lp", "betalpsi", "betalpsi2", "betalp2")

#define constants in the model
const<-list(n.spec=all.data$n.spec, n.survey=all.data$n.survey, n.sp.surv=all.data$n.sp.surv, 
            survey=all.data$survey, nsite=all.data$nsite, 
            species=all.data$species, region=all.data$region, COV=all.data$covmat, 
            first=all.data$first, last=all.data$last, EFF=all.data$effort.sc)

#define data
datn<-list(y=all.data$Y)


####################################################################################################
#### run model in Nimble - this is a slow model! 

#(1) set up model
model <- nimbleModel(Cerrado_FlorestaQ_long, constants = const, data=datn, 
                     inits=inits(), check = FALSE)

#(2) Compile model in c++
cmodel <- compileNimble(model)       

# (3) Configure MCMC 
conf.mcmc<-configureMCMC(model, monitors = params, thin=10, monitors2=c("BP", "bp"), thin2=10)

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

# (6) Run 
samp <- runMCMC(cmcmc, niter = 100000, nburnin = 50000, nchains=3, inits = inits) 

## summary table for everything in "params" vector
summ<-MCMCsummary(samp$samples)
#
## summary table for BP values, with higher thinning
Bp<-MCMCsummary(samp$samples2)
 
##For results interpretation, English species name in same order as data
eng.names<- c("Brazilian tapir", "Short-eared dog", "Giant anteater", "Pampas deer", "Giant armadillo", "Capybara",
              "Collared peccary", "Crab-eating fox", "Marsh deer", "Brazilian porcupine ", "Black agouti", 
              "Green acouchi", "Greater grison", "Maned wolf", "Tayra", "Ocelot", "Crab-eating raccoon", "Margay",
              "Red brocket deer", "Southern tamandua", "Jaguarundi", "Common opossum", "Spotted paca", "Puma",
              "Jaguar", "South American coati", "White-lipped peccary", "Naked-tailed armadillo", "Brown brocket deer",
              "Dasypus spp.", "Bush dog")
