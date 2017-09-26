###
### Setup, functions for simulations for Roettger et al. data
###
### M. Sonderegger & J. Kirby
## last update 9/2017
##
## Requires: data file, called roettgerEtAlData.csv
##
## To run:
## - Source file, and execute simulations using the appropriate function call -- as in runSimsJphon.R.


library(lme4)
library(arm)
library(plyr)
library(dplyr)
library(ggplot2)
library(simr)


## name of effect of interest
## 'voiceless' for Roettger et al. data
effectOfInterest <- 'voiceless'


## names of factors in the dataset
factorVars <- c("place", "vowel", "accent_type", "prosodic_boundary")

####
## 1. PRELIMINAIRES
####

## load data
E1 = read.csv("roettgerEtAlData.csv",comment.char="")

## REMOVE THIS ROW, there are NAs for variables
E1 <- E1[-481,]

## turn items, subjects into factors
E1$item_pair = as.factor(E1$item_pair)
E1$subject = as.factor(E1$subject)

## MS: centered voicing var
E1$voiceless <- rescale(E1$voicing)

## WHATEVER IS DONE ABOVE, name dataset 'data'

data <- E1

## fit original model
E1.mdl = lmer(vowel_dur ~ voiceless  +        		# critical fixed effect
                  accent_type + prosodic_boundary +				# prosodic control variables
                  place + vowel + 								# phonological control variables
                  norming_voiceless_count + 						# norming
                  (1+voiceless|subject) + (1+voiceless|item_pair),
              data=E1)

## THIS IS THE NAME OF THE MODEL TO BE USED THROUGHOUT
origMod <- E1.mdl



####
## 2. PROPERTIES OF THE DATASET AND ORIGINAL FITTED MODEL
####

## number of subjects, items, observations
nSubject <- nrow(ranef(origMod)$subject)
nItem <- nrow(ranef(origMod)$item_pair)
n <- nrow(model.frame(origMod))


## NOT DONE: VARIANCES


## original effect size
origEffectSize <- fixef(origMod)[effectOfInterest][[1]]

origNRep <- 1

#####
### 3. FUNCTIONS FOR POWER SIMULATIONS
#####

## Function to prep a new dataset, with smaller/larger sample size
## and changed effect size
##
## NB: for simr functions, weirdly, changing the dataset means changing 
## what's in the model when you use getData on it. so, we actually change the *model*.
##
## mod: original model, which will be re-run for new dataset and 
## with new effect size in your simulation
## nS, nI, nRep: sample size params for new dataset
## beta: effect size you want to get power for, for effect of interest
sampleNewDataset <- function(mod, nS = nSubject, nI = nItem, nRep  = origNRep,
                    beta = origEffectSize){
    
    newMod <- mod
    
    ## this is the dataframe used to fit the original model
    data <- model.frame(newMod)
    
    ## this changes effect size of effect we're interested in.
    fixef(newMod)[[effectOfInterest]] <- beta
    
    ## 
    ## begin process of changing sample size for subjects and items
    ##
        
    ## If we want to do fewer subjects or items than in the original dataset,
    ## we need to choose a subset of subjects/items s.t. we don't fall into one of
    ## two traps:
    ## 1. not having all values of factors present (that were in original dataset)
    ## 2. having a model matrix that's not full rank
    
    ## these vars flag if we will need to get a subset of subjects or items
    chooseSubjects <- ifelse(nS < nSubject, TRUE, FALSE)
    chooseItems <- ifelse(nI < nItem, TRUE, FALSE)
    
    ## keep choosing a subset of items until we get one that does *not* have just a 
    ## single value for item-level variables (this gives an error about contrasts not 
    ## being able to be applied to factors with 1 level), and also
    ## doesn't get the model matrix downsized (excluding a fixed-effect term)
    ## because of something like choosing items s.t. place can be in part predicted
    ## based on vowel.
    selected <- FALSE
    while(!selected){
        
        ## choose subjects
        if(chooseSubjects){
            newSubjects <- sample(levels(data$subject), nS, replace = FALSE)
        } else{
            newSubjects <- levels(data$subject)
        } 
        
        ## choose items
        if(chooseItems){
            newItems <- sample(levels(data$item_pair), nI, replace = FALSE)
        } else{
            newItems <- levels(data$item_pair)
        }  
        
        ## check that this isn't a subset where there is only one value for item-level variables
        selected = TRUE
        temp <- subset(data, subject%in%newSubjects & item_pair %in% newItems)
        
        for(var in factorVars){
            if(length(unique(temp[[var]]))==1){
                selected = FALSE
                cat("issue with var", var, "\n")
            }
        }
        ## have to check if the resulting model has fewer fixed effects
        ## than original model, which can happen due to data not being
        ## fully crossed
        if(chooseSubjects || chooseItems){
            tempMod <- lmer(formula(mod), temp)
            if(length(fixef(tempMod))<length(fixef(mod))){
                cat("fewer fixed effects",'\n')
                selected = FALSE
            }
        }
        
    }
    ## these are messages you should see a lot for small enough subject or item numbers
    if(chooseSubjects){cat("chose working subset of subjects",'\n')}
    if(chooseItems){cat("chose working subset of items",'\n')}
    
    ## subset the data to just the subset of subjects and items chosen, and
    ## make this the dataframe returned by the simr function getData
    data <- filter(data, subject%in%newSubjects & item_pair %in% newItems)
    getData(newMod) <- data
    
    ## if we want a *larger* sample size for subjects or items than in original data, just
    ## use simR extension method (which I think just repeats items/subjects as necessary)
    if(nS>nSubject){
            newMod <- extend(newMod, along='subject', n=nS)
    }
    if(nI>nItem){
        newMod <- extend(newMod, along='item_pair', n=nI)
    }
        
    
    ## dealing with changing nRep is easy, and we
    ## only have to worry about the case where nRep > origNRep, because the
    ## latter is 1. NB: I'm not sure if this will work if we have a case where in the 
    ## original data nRep > 1.
    if(nRep!=origNRep){
        newMod <- extend(newMod, within='subject+item_pair', n=nRep)
    }
    
    return(newMod)
}

## Run 
## Function to run a simulation for one set of parameter values
##
## parameters::
## nS = # subjects
## nI = # items (= pairs of voiced/voiceless words)
## nRep = # repetitions (per subject:item)
## beta: effect size (for R et al data: difference between voiced and voiceless) 
##
runSim <- function(mod, nS = nSubject, nI = nItem, nRep  = origNRep,
                   beta = origEffectSize, nSims=nSims){
    results <- data.frame()
    
    ## changes effect size (beta) if applicable, and 
    ## changes nS, nI, or nRep.
    newMod <- sampleNewDataset(mod, nS=nS, nI=nI, nRep=nRep,beta=beta)
    
    ## dummy function to fit new models and put info summarizing them in dataframe
    tempFun <- function(dummy){
        
        ## simr function: simulate the response variable for this new dataset, once
        y <- doSim(newMod)
        
        ## simr function: refit the model using the new dataset with new response variable
        z <- doFit(y, newMod)
        
        ## get the coeff value of interest in simulated model
        newBeta <- fixef(z)[effectOfInterest]
        
        ## p value for z test of  effect of interest
        lrTestP <- as.numeric(doTest(z, fixed(effectOfInterest, 'lr')))
        lrTestModSelect=ifelse(lrTestP<0.05, 'yes', 'no')
        
        ## dataframe just giving beta value, p value, and yes/no for 'is this term sig?'
        results <- data.frame(beta=newBeta, 
                              lrTestP = lrTestP,
                              ## asks: " ", based on LR test?
                              lrTestModSelect = lrTestModSelect)
    }
    
    ## do this nSims time, in parallel
    results <- ldply(seq(1,nSims), tempFun, .parallel = useParallel)
    
    ## now add in parameter values for this run
    results <- data.frame(results, nS=nS, nI=nI, nRep=nRep, n=nrow(getData(newMod)),trueBeta=beta, nSims=nSims)
    
    return(results)
}

runSweep <- function(mod, betaVals, nSVals, nIVals, nRepVals, nSims, saveFName){
    sweepRuns <- data.frame()

    ptm <- proc.time()
    
    for(beta in betaVals){
        cat("beta =", beta, "\n")
        for(nS in nSVals){
            cat("nS =", nS, "\n")
            for(nI in nIVals){
                cat("nI =", nI, "\n")
                
                for(nRep in nRepVals){
                    cat("nRep =", nRep, "\n")
                    
                    sweepRuns <- rbind(sweepRuns, runSim(origMod, nRep=nRep, nS=nS, nI=nI, beta=beta, nSims = nSims))
                }
            }
        }
    }
    t <- (proc.time() - ptm)
    cat(t[3])
    
    saveRDS(sweepRuns, file= saveFName)
    
    return(sweepRuns)
}





