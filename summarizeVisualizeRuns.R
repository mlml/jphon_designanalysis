## 
## script to generate figures in JPhon paper "Mixed-effects design analysis for experimental phonetics"
## 
## March 2018

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

## set to TRUE to actually print figures
##
## NB for users besides Kirby & Sonderegger: If set to TRUE, must change path where figures are printed from '"../jphon/jphon_draft/' below, to
## something appropriate for your machine
##
printFigures <- TRUE

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
        dplyr::summarise(nSig=sum(lrTestModSelect=='yes'), n=n()) %>% 
        group_by(nS,nI,nRep, trueBeta, n, nSig) %>% 
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
        dplyr::summarise(nSig=sum(lrTestModSelect=='yes'), n=n(),
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
        dplyr::summarise(nSig=sum(lrTestModSelect=='yes'), n=n(),
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
##
## these include the plots in the JPhon paper

## load data for 1.5k runs
runs <- readRDS("runs_16Mar18_nRuns1500.rds")

## summarize power and Type M/S error for full simulation runs
powerAll <- powerSummary(runs)
msAll <- typeMSSummary(runs)
## factor relabeling for plots
msAll$sig <- factor(factor(msAll$sig, levels=c('ns', 'sig', 'both')), labels=c('n.s. only', 'sig. only', 'both'))

## "low power" subset
sub1 <- subset(runs, nS==6 & nI==10 & nRep==1)
## "medium power" subset
sub2 <- subset(runs, nS==10 & nI==15 & nRep==2)
## "high power" subset
sub3 <- subset(runs, nS==18 & nI==25 & nRep==5)

## "low power" regime 
## calculate/summarize power
subDf1 <- powerSummary(sub1)
## same for Type M/S error
subMsDf1 <- typeMSSummary(sub1)
subDf1$powerClass = 'low'
subMsDf1$powerClass = 'low'


## "mid power" regime
subDf2 <- powerSummary(sub2)
subMsDf2 <- typeMSSummary(sub2)
subDf2$powerClass = 'mid'
subMsDf2$powerClass = 'mid'


## "high power" regime
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

## Type M/S summary conditioned on significance classes:
msAll2 <- typeMSSummary2(runs)

msAll2 <- full_join(msAll2, powerAll)
msAll2 <- msAll2 %>% mutate(nTot=nS*nI*nRep) 
msAll2$nTot.cut <- cut(msAll2$nTot,breaks = 4)
msAll2$power.cut <- cut(msAll2$power,breaks = c(-0.01,0.25,0.5,0.75,1))


# plot true effect size vs power
powerRegionPlot <- ggplot(aes(x=trueBeta, y=power), data=subDf) + 
    geom_line(aes(color=powerClass)) + 
    geom_ribbon(aes(fill=powerClass, ymin=power_lower, ymax=power_upper), alpha=0.25) + geom_hline(yintercept=0.8, lty=2) + ylim(0,1) +
    xlab("True effect size (ms)") + ylab("Power") + geom_hline(yintercept=0.8, lty=2) + theme_set(theme_bw())  + scale_color_discrete("Sample size") + scale_fill_discrete("Sample size")


# plot true effect size vs Type M error, different kinds
typeMRegionPlot <- ggplot(aes(x=trueBeta, y=typeM_est), data = subMsDf) + 
    geom_line(aes(color=powerClass)) + 
    geom_ribbon(aes(ymin=typeM_lower, ymax=typeM_upper, fill=powerClass), alpha=0.25) +
    xlab("True effect size (ms)") + ylab("Type M error")  +
    geom_hline(yintercept=1, lty=2) +  theme_set(theme_bw())  + 
    facet_wrap(~sig) + theme_set(theme_bw()) + scale_color_discrete("Sample size") + scale_fill_discrete("Sample size") + theme(strip.text = element_text(size=10), axis.text = element_text(size=10))

# plot true effect size vs Type S error, different kinds

## HACK: first, just exclude all data for which typeS based on a binomial test where number of hits is
## < 10 (there are some cases where there's only a few n.s. results)
## 
## Justification: when you try to estimate p from a binomial sample where there are no hits (so estimated p is 0) and the 
## total n is not large, the confidence intervals are huge.  Thus, these estimated "zero" p values aren't really informative
## anyway, and also cause confusion when plotted for readers. Their absence should be noted in any presentation of figures.
##
tempDf <- filter(subDf, (n-nSig)<25) %>% ungroup() %>% dplyr::select(nS, nI, nRep, trueBeta, powerClass)
tempDf$sig<-'n.s. only'
tempDf$exclude <- 'yes'
subMsDf.sub <- filter(full_join(subMsDf, tempDf), is.na(exclude))
subMsDf.sub$sig <- factor(subMsDf.sub$sig, levels=c('n.s. only', 'sig. only', 'both'))


## now make the plot
typeSRegionPlot <- ggplot(aes(x=trueBeta, y=typeS_est), data = subMsDf.sub) + 
    geom_line(aes(color=powerClass)) + 
    geom_ribbon(aes(ymin=typeS_lower, ymax=typeS_upper, fill=powerClass), alpha=0.25) + theme_set(theme_bw())  + 
    xlab("True effect size") + ylab("Type S error") + xlab("True effect size (ms)")  + geom_hline(yintercept=0, lty=2) + 
    facet_wrap(~sig) +  theme_set(theme_bw()) + scale_color_discrete("Sample size") + scale_fill_discrete("Sample size") + ylim(0,1) + theme(strip.text = element_text(size=10), axis.text = element_text(size=10))


## THREE PLOTS FOR APPENDIX, showing all simulations, each goes on one full page:

## 1. plot of power versus all variables varied
powerFullPlot <- ggplot(aes(x=trueBeta, y=power, group=nRep), data=powerAll) + 
    geom_line(aes(color=nRep)) + facet_grid(nS~nI, labeller=label_both) + 
    geom_hline(aes(yintercept=0.8), lty=2) +  theme_set(theme_bw())  + 
    xlab("True effect size (ms)") + ylab("Power")


## 2. plot of type M error versus versus all variables varied: unconditioned, or conditioned on significance
typeMFullPlot <- ggplot(aes(x=trueBeta, y=typeM_est, group=interaction(nRep, sig)), 
                        data=msAll) + 
    #filter(msAll, sig%in%c('both','sig. only')))+ 
    geom_line(aes(color=nRep, lty=sig)) + 
    facet_grid(nS~nI, labeller=label_both) +
    xlab("True effect size (ms)") + ylab("Type M error") + theme_set(theme_bw())  + 
    geom_hline(aes(yintercept=1), lty=2) 


## 3. plot of type S error versus versus all variables varied: unconditioned, or conditioned on significance
typeSFullPlot <- ggplot(aes(x=trueBeta, y=typeS_est, group=interaction(nRep, sig)), data=msAll) +
    #filter(msAll, sig%in%c('both','sig. only'))) + 
    geom_line(aes(color=nRep, lty=sig)) + 
    facet_grid(nS~nI, labeller=label_both) +
    xlab("True effect size (ms)") + ylab("Type S error") + theme_set(theme_bw())  + 
    geom_hline(aes(yintercept=0), lty=2) 



## THREE PLOTS FOR TEXT  that are more readable subset versions, for the text:

## 1. power
powerAll$nRep <- factor(powerAll$nRep)
powerAll.sub <- filter(powerAll, nRep%in%c(1,3,5) & nI%in%c(10,20,30) & nS%in%c(6,14,22))

powerFullPlot.sub <- ggplot(aes(x=trueBeta, y=power, group=nRep), data=powerAll.sub) + 
    geom_line(aes(color=nRep), size=1) + facet_grid(nS~nI, labeller=label_both) + 
    geom_hline(aes(yintercept=0.8), lty=2) +  theme_set(theme_bw())  + 
    xlab("True effect size (ms)") + ylab("Power")

## 2. type M
msAll$nRep <- factor(msAll$nRep)
#msAll.sub <- filter(msAll, nRep%in%c(1,3,5))
msAll.sub <- filter(msAll, nRep%in%c(1,6) & nI%in%c(10,20,30) & nS%in%c(6,14,22))

typeMFullPlot.sub <- ggplot(aes(x=trueBeta, y=typeM_est, group=interaction(nRep, sig)), 
                        data=msAll.sub) + 
    #filter(msAll, sig%in%c('both','sig. only')))+ 
    geom_line(aes(lty=nRep, color=sig),size=1) + 
    facet_grid(nS~nI, labeller=label_both) +
    xlab("True effect size (ms)") + ylab("Type M error") + theme_set(theme_bw())  + 
    geom_hline(aes(yintercept=1), lty=2) 

## 3. type S
typeSFullPlot.sub <- ggplot(aes(x=trueBeta, y=typeS_est, group=interaction(nRep, sig)), data=msAll.sub) +
    #filter(msAll, sig%in%c('both','sig. only'))) + 
    geom_line(aes(color=sig, lty=nRep), size=1) + 
    facet_grid(nS~nI, labeller=label_both) +
    xlab("True effect size (ms)") + ylab("Type S error") + theme_set(theme_bw())  + 
    geom_hline(aes(yintercept=0), lty=2)  + ylim(-0.01,0.35)


## power versus Type S error for 0.05<p<0.20 and for p>0.20 results
typeSPowerPlot1 <- ggplot(aes(x=power, y=typeS_est), 
                          data=filter(msAll2, sig=='both' & lrTestP.cut%in%c('(0.05,0.2]', '(0.2,1]'))) + 
    geom_point(size=0.05, alpha=0.5) +
    geom_smooth(size=1) +
    geom_hline(aes(yintercept=0), lty=2) + theme_set(theme_bw())  + 
    facet_grid(~lrTestP.cut, scales='free_y') + ylab("Type S error") + ylim(0, 0.4) + theme(panel.spacing = unit(1, "lines"))

## power versus Type M error for 0.05<p<0.20 and for p>0.20 results
typeMPowerPlot1 <- ggplot(aes(x=power, y=typeM_est),
                          data=filter(msAll2, sig=='both' & lrTestP.cut%in%c('(0.05,0.2]', '(0.2,1]'))) + 
    #geom_linerange(aes(ymin=typeM_lower, ymax=typeM_upper), size=0.1) +
    geom_point(size=0.05, alpha=0.5) +
    geom_smooth() +  
    facet_grid(~lrTestP.cut, scales='free_y') +  theme_set(theme_bw())  + 
    geom_hline(aes(yintercept=1), lty=2)+ ylab("Type M error") + theme(panel.spacing = unit(1, "lines"))


## power versus Type S error for 'large' and 'small' effect sizes, conditioned on significance and unconditioned
typeSPowerPlot2 <- ggplot(aes(x=power, y=typeS_est), data=filter(msAll,  trueBeta%in%c(-10,-4))) + 
    geom_point(aes(color=sig), size=0.1) +
    #geom_linerange(aes(ymin=typeS_lower, ymax=typeS_upper, color=sig), size=0.2) +
    geom_smooth(aes(color=sig)) +  
    facet_grid(~trueBeta, scales='free_y') +  theme_set(theme_bw())  + 
    geom_hline(aes(yintercept=0), lty=2) +ylab("Type S error") +
    xlab("Power") + theme(panel.spacing = unit(1, "lines")) + scale_x_continuous(limits=c(0,1))

## power versus Type M error for 'large' and 'small' effect sizes, conditioned on significance and unconditioned
typeMPowerPlot2 <- ggplot(aes(x=power, y=typeM_est), data=filter(msAll,  trueBeta%in%c(-10,-4))) + 
         #geom_linerange(aes(ymin=typeM_lower, ymax=typeM_upper, color=sig), size=0.2) +
    geom_point(aes(color=sig), size=0.1) +
    geom_smooth(aes(color=sig)) +  
    facet_grid(~trueBeta, scales='free_y') +  theme_set(theme_bw())  + 
    geom_hline(aes(yintercept=1), lty=2) + ylab("Type M error") +
    xlab("Power") + theme(panel.spacing = unit(1, "lines")) + scale_x_continuous(limits=c(0,1))

## change path from '../jphon/jphon_draft' for this to run on your machine
if(printFigures){
    ggsave(powerRegionPlot, file="../jphon/jphon_draft/powerRegionPlot.pdf", width=5.5, height=4)
    ggsave(typeMRegionPlot, file="../jphon/jphon_draft/typeMRegionPlot.pdf", width=8, height=3)
    ggsave(typeSRegionPlot, file="../jphon/jphon_draft/typeSRegionPlot.pdf", width=8, height=3)
    
    ggsave(powerFullPlot, file="../jphon/jphon_draft/powerFullPlot.pdf", width=9, height=7)
    ggsave(typeMFullPlot, file="../jphon/jphon_draft/typeMFullPlot.pdf", width=9, height=7)
    ggsave(typeSFullPlot, file="../jphon/jphon_draft/typeSFullPlot.pdf", width=9, height=7)
    
    ggsave(powerFullPlot.sub, file="../jphon/jphon_draft/powerFullPlot_sub1.pdf", width=7, height=5)
    ggsave(typeMFullPlot.sub, file="../jphon/jphon_draft/typeMFullPlot_sub1.pdf", width=7, height=5)
    ggsave(typeSFullPlot.sub, file="../jphon/jphon_draft/typeSFullPlot_sub1.pdf", width=7, height=5)
   
    typeSMPowerPlot1 <- ggarrange(typeSPowerPlot1, typeMPowerPlot1, ncol=1, nrow=2, align="hv")
    ggsave(typeSMPowerPlot1, file="../jphon/jphon_draft/typeSMPowerPlot1.pdf", width=5, height=5)
    
    
    #ggsave(typeSPowerPlot1, file="../jphon/jphon_draft/typeSPowerPlot1.pdf", width=5,height=2.5)
    #ggsave(typeMPowerPlot1, file="../jphon/jphon_draft/typeMPowerPlot1.pdf", width=5,height=2.5)
    
    typeSMPowerPlot2 <- ggarrange(typeSPowerPlot2, typeMPowerPlot2, ncol=1, nrow=2, align="hv")
    ggsave(typeSMPowerPlot2, file="../jphon/jphon_draft/typeSMPowerPlot2.pdf", width=6, height=5)
    #ggsave(typeSPowerPlot2, file="../jphon/jphon_draft/typeSPowerPlot2.pdf", width=5.05,height=2)
    #ggsave(typeMPowerPlot2, file="../jphon/jphon_draft/typeMPowerPlot2.pdf", width=5,height=2)
}

