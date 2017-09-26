## 
## script to generate figures in JPhon submitted paper.
## 

library(dplyr)
library(tidyr)
library(ggplot2)

## set to TRUE to actually print figures
printFigures <- FALSE

# 1. FUNCTIONS  ---------------------------------------------------

## first functions to summarize runs which result from runSweep

## helper functions to give estimates and errorbars for binomial and normally-distributed
## data
binomEst <- function(n1, n){
    if(n>0){
        temp <- binom.test(n1, n)
        return(c(temp$conf.int[1], temp$conf.int[2], as.numeric(temp$estimate)))
    } else{
        return(c(NA, NA, NA))
    }
}

normalEst <- function(x){
    m <- mean(x, na.rm=TRUE)
    se <- sd(x, na.rm = TRUE)/sqrt(length(which(!is.na(x))))
    return(c(m, m-1.96*se, m+1.96*se))
}

## function to compute summary df for power -- estimate, and lower and upper bounds
## for 95% CI --
## as a function of all input vars, just using LR model selection strategy
##
powerSummary <- function(x){
    
    powerVals <- function(nSig, n){
        temp <- binom.test(nSig, n)
        return(c(temp$conf.int[1], temp$conf.int[2], as.numeric(temp$estimate)))
    }
    
    summaryDf <- x %>% group_by(nS,nI,nRep, trueBeta) %>% 
        summarise(nSig=sum(lrTestModSelect=='yes'), n=n()) %>% 
        group_by(nS,nI,nRep, trueBeta) %>% 
        do(setNames(data.frame(t(powerVals(.$nSig, .$n))), c('power_lower', 'power_upper', 'power')))

    return(summaryDf)                                         
    
}



## make plottable dataframe for just for just typeM and typeS error,
## broken down by whether conditional on significance, non-significance, or using both.
typeMSSummary <-  function(x){
    
    
    
    ## summary for typeM
    tempDf1 <- x %>% mutate(betaRat = abs(beta)/abs(trueBeta), sameSign=(sign(beta)==sign(trueBeta))) %>%
        group_by(nS,nI,nRep, trueBeta) %>%
        do(
            cbind(
                setNames(data.frame(t(normalEst(.$betaRat))),
                         c('typeM_both_est', 'typeM_both_lower', 'typeM_both_upper')),
                setNames(data.frame(t(normalEst(ifelse(.$lrTestModSelect=='yes',.$betaRat,  NA)))), 
                          c('typeM_sig_est', 'typeM_sig_lower', 'typeM_sig_upper')),
                setNames(data.frame(t(normalEst(ifelse(.$lrTestModSelect=='no',.$betaRat,  NA)))), 
                         c('typeM_ns_est', 'typeM_ns_lower', 'typeM_ns_upper'))
            )
        )
                         
              
    tempDf2 <- x %>% mutate(sameSign=(sign(beta)==sign(trueBeta))) %>%
        group_by(nS,nI,nRep, trueBeta) %>%
        summarise(nSig=sum(lrTestModSelect=='yes'), n=n(),
                  nSigRightSign=sum(sameSign & lrTestModSelect=='yes'),
                  nNsRightSign=sum(sameSign & lrTestModSelect=='no'),
                  power = nSig/n()) %>%
        group_by(nS,nI,nRep, trueBeta, power) %>%
        do(
            cbind(
                setNames(data.frame(t(binomEst(.$n-.$nSigRightSign -.$nNsRightSign, .$n))), c('typeS_both_lower', 'typeS_both_upper', 'typeS_both_est')),
                setNames(data.frame(t(binomEst(.$nSig-.$nSigRightSign, .$nSig))), c('typeS_sig_lower', 'typeS_sig_upper', 'typeS_sig_est')),
                setNames(data.frame(t(binomEst(.$n-.$nSig-.$nNsRightSign, .$n-.$nSig))), c('typeS_ns_lower', 'typeS_ns_upper', 'typeS_ns_est'))
            )
        )
     tempDf <- left_join(tempDf1, tempDf2)   
    
    tempDf <- tempDf %>% gather("key", "val", starts_with("type")) %>%
        separate(key, into=c('var','sig','measure')) %>% 
        unite(measure, var, measure) %>% 
        spread(measure, val)
    return(tempDf)
    
}

## basically the same function, but now conditioned on significance classes
typeMSSummary2 <- function(x){
    
    x$lrTestP.cut <- cut(runs$lrTestP, breaks=c(0,0.01,0.05,0.20,1.0))
    
    ## summary for typeM
    tempDf1 <- x %>% mutate(betaRat = abs(beta)/abs(trueBeta), sameSign=(sign(beta)==sign(trueBeta))) %>%
        group_by(nS,nI,nRep,trueBeta, lrTestP.cut) %>%
        do(
            cbind(
                setNames(data.frame(t(normalEst(.$betaRat))),
                         c('typeM_both_est', 'typeM_both_lower', 'typeM_both_upper')),
                setNames(data.frame(t(normalEst(ifelse(.$lrTestModSelect=='yes',.$betaRat,  NA)))), 
                         c('typeM_sig_est', 'typeM_sig_lower', 'typeM_sig_upper')),
                setNames(data.frame(t(normalEst(ifelse(.$lrTestModSelect=='no',.$betaRat,  NA)))), 
                         c('typeM_ns_est', 'typeM_ns_lower', 'typeM_ns_upper'))
            )
        )
    
    tempDf2 <- x %>% mutate(sameSign=(sign(beta)==sign(trueBeta))) %>%
        group_by(nS,nI,nRep, trueBeta, lrTestP.cut) %>%
        summarise(nSig=sum(lrTestModSelect=='yes'), n=n(),
                  nSigRightSign=sum(sameSign & lrTestModSelect=='yes'),
                  nNsRightSign=sum(sameSign & lrTestModSelect=='no')) %>%
        group_by(nS,nI,nRep, trueBeta,lrTestP.cut) %>%
        do(
            cbind(
                setNames(data.frame(t(binomEst(.$n-.$nSigRightSign -.$nNsRightSign, .$n))), c('typeS_both_lower', 'typeS_both_upper', 'typeS_both_est')),
                setNames(data.frame(t(binomEst(.$nSig-.$nSigRightSign, .$nSig))), c('typeS_sig_lower', 'typeS_sig_upper', 'typeS_sig_est')),
                setNames(data.frame(t(binomEst(.$n-.$nSig-.$nNsRightSign, .$n-.$nSig))), c('typeS_ns_lower', 'typeS_ns_upper', 'typeS_ns_est'))
            )
        )
    
    ## TODO: figure out n for each case (both, ns, sig), and power
    ## across the three cases
    
    tempDf <- left_join(tempDf1, tempDf2)   
    
    tempDf <- tempDf %>% gather("key", "val", starts_with("type")) %>%
        separate(key, into=c('var','sig','measure')) %>% 
        unite(measure, var, measure) %>% 
        spread(measure, val)
    return(tempDf)
    
}

#runs <- readRDS("sweepRuns_v1_nRuns250.rds")

#ggplot(aes(x=nTot, y=typeM), data=msAll2) + geom_boxplot(aes(fill=lrTestP.cut)) + facet_grid(~trueBeta, scales='free_y') + geom_hline(aes(yintercept=1), lty=2)
# 
# msAll2 <- typeMSSummary2(runs)
# 
# msAll2 <- full_join(msAll2, powerAll)
# msAll2 <- msAll2 %>% mutate(nTot=nS*nI*nRep) 
# msAll2$nTot.cut <- cut(msAll2$nTot,breaks = 4)
# msAll2$power.cut <- cut(msAll2$power,breaks = c(-0.01,0.25,0.5,0.75,1))


# 2. SAMPLE VISUALIZATIONS ---------------------------------------------------


## load data for one set of runs
runs <- readRDS("sweepRuns_v2_nRuns500.rds")

## "low power" subset
sub1 <- subset(runs, nS==6 & nI==10 & nRep==1)
## "medium power" subset (festschrift was 8, 12 3)
sub2 <- subset(runs, nS==10 & nI==15 & nRep==2)
## "high power" subset (festscrhift was 16, 24, 6)
sub3 <- subset(runs, nS==18 & nI==25 & nRep==5)

#summaryDf <- summarizeRuns(runs)


subDf1 <- powerSummary(sub1)
subMsDf1 <- typeMSSummary(sub1)
subDf1$powerClass = 'low'
subMsDf1$powerClass = 'low'


subDf2 <- powerSummary(sub2)
subMsDf2 <- typeMSSummary(sub2)
subDf2$powerClass = 'mid'
subMsDf2$powerClass = 'mid'

powerAll <- powerSummary(runs)
msAll <- typeMSSummary(runs)
msAll$sig <- factor(factor(msAll$sig, levels=c('ns', 'sig', 'both')), labels=c('n.s. only', 'sig. only', 'both'))


subDf3 <- powerSummary(sub3)
subMsDf3 <- typeMSSummary(sub3)
subDf3$powerClass = 'high'
subMsDf3$powerClass = 'high'

## make dataframes for convenient plotting
subDf <- rbind(subDf1, subDf2, subDf3)
subDf$powerClass <- factor(subDf$powerClass, levels=c("low", "mid", "high"))

subMsDf <- rbind(subMsDf1, subMsDf2, subMsDf3)
subMsDf$powerClass <- factor(subMsDf$powerClass, levels=c("low", "mid", "high"))

subMsDf$sig <- factor(factor(subMsDf$sig, levels=c('ns', 'sig', 'both')), labels=c('n.s. only', 'sig. only', 'both'))

## conditioned on significance classes:
msAll2 <- typeMSSummary2(runs)

msAll2 <- full_join(msAll2, powerAll)
msAll2 <- msAll2 %>% mutate(nTot=nS*nI*nRep) 
msAll2$nTot.cut <- cut(msAll2$nTot,breaks = 4)
msAll2$power.cut <- cut(msAll2$power,breaks = c(-0.01,0.25,0.5,0.75,1))



# plot true effect size vs power
powerRegionPlot <- ggplot(aes(x=trueBeta, y=power), data=subDf) + 
    geom_line(aes(color=powerClass)) + 
    geom_ribbon(aes(fill=powerClass, ymin=power_lower, ymax=power_upper), alpha=0.25) + geom_hline(yintercept=0.8, lty=2) + ylim(0,1) +
    xlab("True effect size (ms)") + ylab("Power") + geom_hline(yintercept=0.8, lty=2) + scale_color_discrete("Sample size") + scale_fill_discrete("Sample size")

# plot true effect size vs Type M error, different kinds
typeMRegionPlot <- ggplot(aes(x=trueBeta, y=typeM_est), data = subMsDf) + 
    geom_line(aes(color=powerClass)) + 
    geom_ribbon(aes(ymin=typeM_lower, ymax=typeM_upper, fill=powerClass), alpha=0.25) +
    xlab("True effect size") + ylab("Type M error")  +
    geom_hline(yintercept=1, lty=2) + 
    facet_wrap(~sig) + scale_color_discrete("Sample size") + scale_fill_discrete("Sample size")

# plot true effect size vs Type S error, different kinds
typeSRegionPlot <- ggplot(aes(x=trueBeta, y=typeS_est), data = subMsDf) + 
    geom_line(aes(color=powerClass)) + 
    geom_ribbon(aes(ymin=typeS_lower, ymax=typeS_upper, fill=powerClass), alpha=0.25) +
    xlab("True effect size") + ylab("Type S error") + xlab("True effect size (ms)")  + geom_hline(yintercept=0, lty=2) + 
    facet_wrap(~sig) +  scale_color_discrete("Sample size") + scale_fill_discrete("Sample size")

## plot of power versus all variables varied
powerFullPlot <- ggplot(aes(x=trueBeta, y=power, group=nRep), data=powerAll) + 
    geom_line(aes(color=nRep)) + facet_grid(nS~nI, labeller=label_both) + 
    geom_hline(aes(yintercept=0.8), lty=2) + 
    xlab("True effect size (ms)") + ylab("Power")

## plot of type M error versus versus all variables varied: unconditioned, or conditioned on significance
typeMFullPlot <- ggplot(aes(x=trueBeta, y=typeM_est, group=interaction(nRep, sig)), 
                        data=msAll) + 
    #filter(msAll, sig%in%c('both','sig. only')))+ 
    geom_line(aes(color=nRep, lty=sig)) + 
    facet_grid(nS~nI, labeller=label_both) +
    xlab("True effect size (ms)") + ylab("Type M error") +
    geom_hline(aes(yintercept=1), lty=2) 


## plot of type S error versus versus all variables varied: unconditioned, or conditioned on significance
typeSFullPlot <- ggplot(aes(x=trueBeta, y=typeS_est, group=interaction(nRep, sig)), data=msAll) +
#filter(msAll, sig%in%c('both','sig. only'))) + 
    geom_line(aes(color=nRep, lty=sig)) + 
    facet_grid(nS~nI, labeller=label_both) +
    xlab("True effect size (ms)") + ylab("Type S error") +
    geom_hline(aes(yintercept=0), lty=2) 

## power versus Type S error for 0.05<p<0.20 and for p>0.20 results
typeSPowerPlot1 <- ggplot(aes(x=power, y=typeS_est), 
                          data=filter(msAll2, sig=='both' & lrTestP.cut%in%c('(0.05,0.2]', '(0.2,1]'))) + 
    geom_smooth() +
    geom_point(size=0.1, alpha=0.5) +
    geom_hline(aes(yintercept=0), lty=2) + 
    facet_grid(~lrTestP.cut, scales='free_y') +ylab("Type S error")

## power versus Type M error for 0.05<p<0.20 and for p>0.20 results
typeMPowerPlot1 <- ggplot(aes(x=power, y=typeM_est),
                          data=filter(msAll2, sig=='both' & lrTestP.cut%in%c('(0.05,0.2]', '(0.2,1]'))) + 
    #geom_linerange(aes(ymin=typeM_lower, ymax=typeM_upper), size=0.1) +
    geom_point(size=0.1, alpha=0.5) +
    geom_smooth() +  
    facet_grid(~lrTestP.cut, scales='free_y') + 
    geom_hline(aes(yintercept=1), lty=2)+ ylab("Type M error")


## power versus Type S error for 'large' and 'small' effect sizes, conditioned on significance and unconditioned
typeSPowerPlot2 <- ggplot(aes(x=power, y=typeS_est), data=filter(msAll,  trueBeta%in%c(-10,-4))) + 
    geom_point(aes(color=sig), size=0.1) +
    #geom_linerange(aes(ymin=typeS_lower, ymax=typeS_upper, color=sig), size=0.2) +
    geom_smooth(aes(color=sig)) +  
    facet_grid(~trueBeta, scales='free_y') + 
    geom_hline(aes(yintercept=0), lty=2) +ylab("Type S error") +
    xlab("Power")

## power versus Type M error for 'large' and 'small' effect sizes, conditioned on significance and unconditioned
typeMPowerPlot2 <- ggplot(aes(x=power, y=typeM_est), data=filter(msAll,  trueBeta%in%c(-10,-4))) + 
         #geom_linerange(aes(ymin=typeM_lower, ymax=typeM_upper, color=sig), size=0.2) +
    geom_point(aes(color=sig), size=0.1) +
    geom_smooth(aes(color=sig)) +  
    facet_grid(~trueBeta, scales='free_y') + 
    geom_hline(aes(yintercept=1), lty=2) + ylab("Type M error") +
    xlab("Power")

if(printFigures){
    ggsave(typeSRegionPlot, file="../jphon/jphon_draft/typeSRegionPlot.pdf", width=6, height=3)
    ggsave(powerFullPlot, file="../jphon/jphon_draft/powerFullPlot.pdf", width=7, height=5)
    ggsave(typeMFullPlot, file="../jphon/jphon_draft/typeMFullPlot.pdf", width=7, height=5)
    ggsave(typeSFullPlot, file="../jphon/jphon_draft/typeSFullPlot.pdf", width=7, height=5)
    
    ggsave(typeSPowerPlot2, file="../jphon/jphon_draft/typeSPowerPlot2.pdf", width=4,height=3)
    ggsave(typeMPowerPlot2, file="../jphon/jphon_draft/typeMPowerPlot2.pdf", width=4,height=3)
    ggsave(typeSPowerPlot1, file="../jphon/jphon_draft/typeSPowerPlot1.pdf", width=4,height=3)
    ggsave(typeMPowerPlot1, file="../jphon/jphon_draft/typeMPowerPlot1.pdf", width=4,height=3)
    
}

