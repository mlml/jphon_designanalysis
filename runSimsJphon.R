## script to actually run simulations, resulting in .rds file
##
## CHANGE NCORES FOR YOUR MACHINE
##


source("powerEs.R")

## Change useParallel, nCores for your machine

useParallel <- TRUE
nCores <- 25

if(useParallel){
    library(doMC)
    registerDoMC(nCores)
    getDoParWorkers()
}

## this takes 4-5 days to run using 25 cores on our server.
## To test before doing a full run, set nSims=10.

output <- runSweep(mod=origMod, betaVals=seq(-2, -10),
         nSVals <- seq(6,26,4),
         nIVals <- seq(10,30,5),
         nRepVals <- seq(1,6),
         nSims=1500,
         saveFName = 'runs_16Mar18_nRuns1500.rds')
