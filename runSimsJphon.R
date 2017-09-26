##
## script to actually run simulations, resulting in .rds file
##
## CHANGE NCORES FOR YOUR MACHINE
##

source("powerEs.R")

## - Change useParallel, nCores for your computer

useParallel <- TRUE
nCores <- 25

if(useParallel){
    library(doMC)
    registerDoMC(nCores)
    getDoParWorkers()
}

# this takse a couple days using 25 cores.
## To test before doing a full run, set nSims=10.
runSweep(mod=origMod, betaVals=seq(-2, -10),
         nSVals <- seq(6,26,4),
         nIVals <- seq(10,30,5),
         nRepVals <- seq(1,6),
         nSims=500,
         saveFName = 'runs_26Sep17_nRuns500.rds')

