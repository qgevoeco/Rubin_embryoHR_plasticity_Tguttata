rm(list = ls())
library(lmerTest) #library(lme4)
library(emmeans)
library(DHARMa)  #<-- for Durbin-Watson autocorrelation test on residuals
library(optimx)  #<-- to use alternative optimizers in lme4

#setwd("")   #<-- FIXME set this to working directory with data

#Load data
hrd <- read.table("data_RubinHRD.txt", header = TRUE)

## Within-subject mean Time (Minutes)
hrd$T2wsmn <- hrd$T1wsmn <- NA  
idMeansT1b <- aggregate(T1 ~ EggID, data = hrd, FUN = mean, na.rm = TRUE)
  hrd$T1wsmn <- idMeansT1b$T1[match(hrd$EggID, idMeansT1b$EggID)]
idMeansT2b <- aggregate(T2 ~ EggID, data = hrd, FUN = mean, na.rm = TRUE)
  hrd$T2wsmn <- idMeansT2b$T2[match(hrd$EggID, idMeansT2b$EggID)]

## Within-subject Time (Minutes) deviations 
hrd$T1wsdev <- with(hrd, T1 - T1wsmn)
hrd$T2wsdev <- with(hrd, T2 - T2wsmn)

## Create Time 1 and Time 2 specific datasets
hrd1 <- hrd[which(!is.na(hrd$HR1)), ]
hrd2 <- hrd[which(!is.na(hrd$HR2)), ]

## Make R factors
hrd <- within(hrd, {
  EggID <- as.factor(EggID)
  NestID <- as.factor(NestID)
  Trt <- as.factor(Trt)
  Order1fac <- as.factor(Order1fac)
  Order2fac <- as.factor(Order2fac)
  })
  
hrd1 <- within(hrd1, {
  EggID <- as.factor(EggID)
  NestID <- as.factor(NestID)
  Trt <- as.factor(Trt)
  Order1fac <- as.factor(Order1fac)
  })
  
hrd2 <- within(hrd2, {
  EggID <- as.factor(EggID)
  NestID <- as.factor(NestID)
  Trt <- as.factor(Trt)
  Order2fac <- as.factor(Order2fac)
  })


# Color palette
## set up a color accessible palette to use for plotting throughout
## checked with `colorblindcheck` that implements model of human color vision
clr <- structure(c(blue = "#03244d", cyan = "#00A08A", orange = "#e86823"), 
  class = "palette", name = "aucfriend")
# check:
## colorblindcheck::palette_check(clr, plot = TRUE)

# Create gradient between colors so one unique per ID
## Number of unique IDs in each treatment
#aggregate(EggID ~ Trt, data = hrd, FUN = function(x) length(unique(x)))
nin <- length(levels(hrd$EggID))
clrg <-  structure(grDevices::colorRampPalette(clr)(nin),
    class = "palette",
    name = "aucfriend heat gradient")
# Make transparent versions for plotting in the "background"
clr_trans <- structure(adjustcolor(clr, alpha.f = 0.4),
  class = "palette",
  name = "aucfriend semi-transparent")
clrg_trans <- structure(adjustcolor(clrg, alpha.f = 0.4),
  class = "palette",
  name = "aucfiend heat gradient semi-transparent")
  




# Used to Calculate AICc
## See p. 66 Burnham and Anderson 2002)
myAICc <- function(object){
  lL <- logLik(object, REML = TRUE)
    K <- attr(lL, "df")
    lL <- c(lL)
  n <- nobs(object)
 -2 * lL + 2 * K * (n / (n - K - 1))
}


################################################################################
R.Version()

packageVersion("lmerTest")
packageVersion("lme4")
citation("lme4")
##############################


## Number of individuals with measurements
### 30% development
length(levels(hrd1$EggID))
### 80% development
length(levels(hrd2$EggID))







###############################################################################
###############################################################################
###############################################################################
####       ######       ########    ####  #####################################
### #####   ####  #####  #######    ###   #####################################
#########    ##  #######  ############   ######################################
#########   ##   #######   ##########   #######################################
#####     ####   #######   #########   ########################################
#########   ###   #####    ########   #########################################
#########    ##   #####   ########   ##########################################
### #####   ####  #####  ########   ###    ####################################
####       ######       ########   ####    ####################################
###############################################################################
###############################################################################
###############################################################################

####  *** 30% Development  ***   

# Correlated random slopes & intercepts
hr1_fullFxd_ML <- lmer(scHR1 ~ Trt + T1wsmn:Trt + T1wsdev:Trt +
	(1 + T1 | EggID) + (1 | NestID),
	data = hrd, REML = FALSE, na.action = na.exclude)
summary(hr1_fullFxd_ML)


## Within-subject mean terms give difference between within- and among-subjects mean slopes when modeling within-subject means and regular (time) covariate
### See van de Pol and Wright 2009
hr1_FxdWithAmongDiff_ML <- lmer(scHR1 ~ Trt + T1wsmn:Trt + T1:Trt +
	(1 + T1 | EggID) + (1 | NestID),
	data = hrd, REML = FALSE, na.action = na.exclude)
summary(hr1_FxdWithAmongDiff_ML)



## Now drop among-individual variation in within-individual plasticity
### refit with REML to compare different varcomp models
hr1_fullFxd_REML <- lmer(scHR1 ~ Trt + T1wsmn:Trt + T1wsdev:Trt +
	(1 + T1 | EggID) + (1 | NestID),
	data = hrd, REML = TRUE, na.action = na.exclude)
hr1_fullFxd_Vint_REML <- lmer(scHR1 ~ Trt + T1wsmn:Trt + T1wsdev:Trt +
	(1 | EggID) + (1 | NestID),
	data = hrd, REML = TRUE, na.action = na.exclude)
### Also fit model with random intercepts and slops UNCORRELATED to compare at once
hr1_fullFxd_VintVslo_REML <- lmer(scHR1 ~ Trt + T1wsmn:Trt + T1wsdev:Trt +
	(1 | EggID) + (0 + T1 | EggID) + (1 | NestID),
	data = hrd, REML = TRUE, na.action = na.exclude)
# Issues with convergence warnings: Nest random intercepts ~0
anova(hr1_fullFxd_REML, hr1_fullFxd_Vint_REML, hr1_fullFxd_VintVslo_REML,
  refit = FALSE)

# Best model
summary(hr1_fullFxd_REML)






##########################
## Is among-individual variation in plasticity different among treatments?
### See this post for hints how to do in `lmer()` 
#### XXX (specifically Ben Bolker's answer further down the page):
#### https://stats.stackexchange.com/questions/489059/obtaining-correlation-between-random-effects-separately-for-2-groups/489181#489181

## Treatment-specific random intercepts and slopes CORRELATED
### Used `lme4::allFit()` after default optimizer gave some warnings/issues
hr1_fullFxd_Trt_REML <- lmer(scHR1 ~ Trt + T1wsmn:Trt + T1wsdev:Trt +
  (dummy(Trt, "cont"):T1 | EggID) + 
  (dummy(Trt, "low"):T1 | EggID) +
  (dummy(Trt, "per"):T1 | EggID) + 
  (1 | NestID),
	data = hrd, REML = TRUE, na.action = na.exclude,
  	control = lmerControl(optimizer = "nlminbwrap"))
summary(hr1_fullFxd_Trt_REML)

 


## Treatment-specific random intercepts and slopes UNCORRELATED
hr1_fullFxd_TrtVintVslo_REML <- lmer(scHR1 ~ Trt + T1wsmn:Trt + T1wsdev:Trt +
  (0 + dummy(Trt, "cont") | EggID) + (0 + dummy(Trt, "cont"):T1 | EggID) + 
  (0 + dummy(Trt, "low") | EggID) + (0 + dummy(Trt, "low"):T1 | EggID) +
  (0 + dummy(Trt, "per") | EggID) + (0 + dummy(Trt, "per"):T1 | EggID) +
  (1 | NestID),
	data = hrd, REML = TRUE, na.action = na.exclude,
  	control = lmerControl(optimizer = "nlminbwrap"))
summary(hr1_fullFxd_TrtVintVslo_REML)

# Is there support for correlated random slopes within each treatment
anova(hr1_fullFxd_Trt_REML, hr1_fullFxd_TrtVintVslo_REML, refit = FALSE)





## Treatment-specific random intercepts and homogeneous random slopes UNCORRELATED
hr1_fullFxd_TrtVintHomVslo_REML <- lmer(scHR1 ~ Trt + T1wsmn:Trt + T1wsdev:Trt + 
  (0 + dummy(Trt, "cont") | EggID) + 
  (0 + dummy(Trt, "low") | EggID) +
  (0 + dummy(Trt, "per") | EggID) +
  (0 + T1 | EggID) +
  (1 | NestID),
	data = hrd, REML = TRUE, na.action = na.exclude,
  	control = lmerControl(optimizer = "nlminbwrap"))
summary(hr1_fullFxd_TrtVintHomVslo_REML)




# Treatment-specific versus homogeneous (across treatments) slopes UNCORRELATED
anova(hr1_fullFxd_TrtVintVslo_REML, hr1_fullFxd_TrtVintHomVslo_REML,
  refit = FALSE)


# Best model
summary(hr1_fullFxd_Trt_REML) 


  

anova(hr1_fullFxd_REML, hr1_fullFxd_Trt_REML, hr1_fullFxd_TrtVintVslo_REML,
  hr1_fullFxd_TrtVintHomVslo_REML, refit = FALSE)
  
myAICc(hr1_fullFxd_REML)
myAICc(hr1_fullFxd_TrtVintHomVslo_REML)
myAICc(hr1_fullFxd_TrtVintVslo_REML)
myAICc(hr1_fullFxd_Trt_REML)

anova(hr1_fullFxd_REML, hr1_fullFxd_Trt_REML, refit = FALSE)





lmm1 <- hr1_fullFxd_Trt_REML


  


















###############################################################################
###############################################################################
###############################################################################
####       ######       ########    ####  #####################################
###   ###   ####  #####  #######    ###   #####################################
##    ###    ##  #######  ############   ######################################
###   ###   ##   #######   ##########   #######################################
#####     ####   #######   #########   ########################################
###   ###   ###   #####    ########   #########################################
##    ###    ##   #####   ########   ##########################################
###   ###   ####  #####  ########   ###    ####################################
####       ######       ########   ####    ####################################
###############################################################################
###############################################################################
###############################################################################

####  *** 80% Development  ***   

#Correlated random slopes & intercepts
hr2_fullFxd_ML <- lmer(scHR2 ~ Trt + T2wsmn:Trt + T2wsdev:Trt +
	(1 + T2 | EggID) + (1 | NestID),
	data = hrd, REML = FALSE, na.action = na.exclude)
summary(hr2_fullFxd_ML)


## Within-subject mean terms give difference between within- and among-subjects mean slopes	
hr2_FxdWithAmongDiff_ML <- lmer(scHR2 ~ Trt + T2wsmn:Trt + T2:Trt +
	(1 + T2 | EggID) + (1 | NestID),
	data = hrd, REML = FALSE, na.action = na.exclude)
summary(hr2_FxdWithAmongDiff_ML)

## Now drop among-individual variation in within-individual plasticity
### refit with REML to compare different varcomp models
hr2_fullFxd_REML <- lmer(scHR2 ~ Trt + T2wsmn:Trt + T2wsdev:Trt +
	(1 + T2 | EggID) + (1 | NestID),
	data = hrd, REML = TRUE, na.action = na.exclude)
hr2_fullFxd_Vint_REML <- lmer(scHR2 ~ Trt + T2wsmn:Trt + T2wsdev:Trt +
	(1 | EggID) + (1 | NestID),
	data = hrd, REML = TRUE, na.action = na.exclude)
### Also fit model with random intercepts and slopes UNCORRELATED to compare at once
hr2_fullFxd_VintVslo_REML <- lmer(scHR2 ~ Trt + T2wsmn:Trt + T2wsdev:Trt +
	(1 | EggID) + (0 + T2 | EggID) + (1 | NestID),
	data = hrd, REML = TRUE, na.action = na.exclude)


anova(hr2_fullFxd_REML, hr2_fullFxd_Vint_REML, hr2_fullFxd_VintVslo_REML,
  refit = FALSE)







##########################
## Is among-individual variation in plasticity different among treatments?

## Treatment-specific random intercepts and slopes CORRELATED
### Used `lme4::allFit()` after default optimizer gave some warnings/issues
hr2_fullFxd_Trt_REML <- lmer(scHR2 ~ Trt + T2wsmn:Trt + T2wsdev:Trt + 
  (dummy(Trt, "cont"):T2 | EggID) + 
  (dummy(Trt, "low"):T2 | EggID) +
  (dummy(Trt, "per"):T2 | EggID) +
  (1 | NestID),
	data = hrd, REML = TRUE, na.action = na.exclude,
  	control = lmerControl(optimizer = "nlminbwrap"))
summary(hr2_fullFxd_Trt_REML)

#afit_Full <- lme4::allFit(hr2_fullFxd_Trt_REML)
#summary(afit_Full)
## ALL but one give extremely similar results


## Treatment-specific random intercepts and slopes UNCORRELATED
hr2_fullFxd_TrtVintVslo_REML <- lmer(scHR2 ~ Trt + T2wsmn:Trt + T2wsdev:Trt +
  (0 + dummy(Trt, "cont") | EggID) + (0 + dummy(Trt, "cont"):T2 | EggID) + 
  (0 + dummy(Trt, "low") | EggID) + (0 + dummy(Trt, "low"):T2 | EggID) +
  (0 + dummy(Trt, "per") | EggID) + (0 + dummy(Trt, "per"):T2 | EggID) +
  (1 | NestID),
	data = hrd, REML = TRUE, na.action = na.exclude,
  	control = lmerControl(optimizer = "nlminbwrap"))
summary(hr2_fullFxd_TrtVintVslo_REML)


# Is there support for correlated random slopes within each treatment
anova(hr2_fullFxd_Trt_REML, hr2_fullFxd_TrtVintVslo_REML, refit = FALSE)




## Treatment-specific random intercepts and homogeneous random slopes UNCORRELATED
hr2_fullFxd_TrtVintHomVslo_REML <- lmer(scHR2 ~ Trt + T2wsmn:Trt + T2wsdev:Trt +
  (0 + dummy(Trt, "cont") | EggID) + 
  (0 + dummy(Trt, "low") | EggID) +
  (0 + dummy(Trt, "per") | EggID) +
  (0 + T2 | EggID) + (1 | NestID),
	data = hrd, REML = TRUE, na.action = na.exclude,
  	control = lmerControl(optimizer = "nlminbwrap"))
summary(hr2_fullFxd_TrtVintHomVslo_REML)



# Treatment-specific versus homogeneous (across treatments) slopes
anova(hr2_fullFxd_TrtVintVslo_REML, hr2_fullFxd_VintVslo_REML, hr2_fullFxd_TrtVintHomVslo_REML, hr2_fullFxd_Trt_REML, refit = FALSE)


summary(hr2_fullFxd_Trt_REML) 





anova(hr2_fullFxd_REML, hr2_fullFxd_Trt_REML, hr2_fullFxd_TrtVintVslo_REML,
  hr2_fullFxd_TrtVintHomVslo_REML, refit = FALSE)
  
myAICc(hr2_fullFxd_REML)
myAICc(hr2_fullFxd_TrtVintHomVslo_REML)
myAICc(hr2_fullFxd_TrtVintVslo_REML)
myAICc(hr2_fullFxd_Trt_REML)

anova(hr2_fullFxd_REML, hr2_fullFxd_Trt_REML, refit = FALSE)






lmm2 <- hr2_fullFxd_Trt_REML















################################################################################

#		RESULTS



################################################################################





#######################################################
######	CONFIDENCE INTERVALS - BOOTSTRAP    ###########
#######################################################
# Create a data.frame to get model predictions
## Find the minimum time of measurement in the data and use that
### use minimum observed so not predicting outside of data
## predict across mean (across all individuals) within-subject mean
maxT1 <- max(hrd$T1, na.rm = TRUE)
minT1 <- min(hrd$T1, na.rm = TRUE)
mean1WSmn <- mean(hrd$T1wsmn, na.rm = TRUE)

maxT2 <- max(hrd$T2, na.rm = TRUE)
#XXX check that doing same time across HR1 and HR2 (dev 30% vs 80%)
stopifnot(minT1 == min(hrd$T2, na.rm = TRUE))
minT2 <- min(hrd$T2, na.rm = TRUE)
mean2WSmn <- mean(hrd$T2wsmn, na.rm = TRUE)



np <- 100  #<-- number of prediction points
cidata1 <- data.frame(Trt = rep(levels(hrd$Trt), each = np),
  T1wsmn = rep(mean1WSmn, np * length(levels(hrd$Trt))),
  T1wsdev = rep(seq(minT1 - mean1WSmn, maxT1 - mean1WSmn, length.out = np),
    length(levels(hrd$Trt))))
cidata2 <- data.frame(Trt = rep(levels(hrd$Trt), each = np),
  T2wsmn = rep(mean2WSmn, np * length(levels(hrd$Trt))),
  T2wsdev = rep(seq(minT2 - mean2WSmn, maxT2 - mean2WSmn, length.out = np),
    length(levels(hrd$Trt))))

# make index to predictions/CIs for minimum time of measurement in each Treatment
minTind <- seq(from = 1, by = np, length.out = 3)
  # check
  stopifnot(all(cidata1$T1wsdev[minTind] == (minT1 - mean1WSmn)))
  stopifnot(all(cidata2$T2wsdev[minTind] == (minT2 - mean2WSmn)))
   
  
# Bootstrap confidence intervals
parmOut1 <- function(mod){
  fxd <- fixef(mod)
    # Convert deviations to actual values of each paramter
    # transform back to original HR scale
    hrcenter <- mean(hrd1$HR1)
    hrscale <- sd(hrd1$HR1)
    ### convert back to original HR scale
    wsdevs <- fxd[7:9]
    wsdevs_unsc <- wsdevs * hrscale
  vc <- as.data.frame(VarCorr(mod),
      order = "cov.last")[, c("vcov", "sdcor")][cbind(seq(11),
                                                    c(rep(c(1,1,2), 3), 1, 1))]
    names(vc) <- c(paste0(rep(c("Vu0.", "Vu1.", "r_u0,u1."), 3),
      rep(levels(hrd$Trt), each = 3)), "nest", "res")
  # predict Treatment averages at minimum time of measurement
  pred <- predict(mod, newdata = cidata1, re.form = ~ 0)
    ## Now convert back to original HR scale
    pred_unsc <- (pred * hrscale) + hrcenter
 return(c(fxd, vc, wsdevs_unsc, pred, pred_unsc))      
}

parmOut2 <- function(mod){
  fxd <- fixef(mod)
    # Convert deviations to actual values of each paramter
    # transform back to original HR scale
    hrcenter <- mean(hrd2$HR2)
    hrscale <- sd(hrd2$HR2)
    ### convert back to original HR scale
    wsdevs <- fxd[7:9]
    wsdevs_unsc <- wsdevs * hrscale
  vc <- as.data.frame(VarCorr(mod),
      order = "cov.last")[, c("vcov", "sdcor")][cbind(seq(11),
                                                    c(rep(c(1,1,2), 3), 1, 1))]
    names(vc) <- c(paste0(rep(c("Vu0.", "Vu1.", "r_u0,u1."), 3),
      rep(levels(hrd$Trt), each = 3)), "nest", "res")
  # predict Treatment averages at minimum time of measurement
  pred <- predict(mod, newdata = cidata2, re.form = ~ 0)
    ## Now convert back to original HR scale
    pred_unsc <- (pred * hrscale) + hrcenter
 return(c(fxd, vc, wsdevs_unsc, pred, pred_unsc))      
}








# 30% Development
if(FALSE){  #<-- turned "off" so doesn't take a long time when run entire script
# nsim=1000|ncpus=2 takes about 3 min.
# nsim=10000|ncpus=2 takes about 32 min.
nboot <- 10000
Sys.time()
system.time(boothr1 <- bootMer(hr1_fullFxd_Trt_REML, FUN = parmOut1,
  nsim = nboot, type = "parametric", use.u = FALSE,
  parallel = "multicore", ncpus = 3))
  
confint(boothr1, type = "perc")

save("boothr1",
  file = paste0("hr1_bootstrap", round(nboot/1000, 1), "k.rdata"))
Sys.time()
} else{  #<-- end turning "off" bootstrap
    load("hr1_bootstrap10k.rdata")
  }



# Have a peek at bootstrap distributions for model parameters
## Change `Sys.sleep()` to pause more between plots
sapply(seq(20), FUN = function(x){ plot(boothr1, index = x);
    mtext(dimnames(boothr1$t)[[2L]][x], side = 3, outer = TRUE, line = -2);
    Sys.sleep(0.002); invisible(x)})

dev.off()

hr1Ests <- cbind(Est = parmOut1(hr1_fullFxd_Trt_REML),
  confint(boothr1, type = "perc"))




# 80% Development
if(FALSE){  #<-- turned "off" so doesn't take a long time when run entire script
# nsim=10000|ncpus=2 takes about 32 min.
nboot <- 10000
Sys.time()
system.time(boothr2 <- bootMer(hr2_fullFxd_Trt_REML, FUN = parmOut2,
  nsim = nboot, type = "parametric", use.u = FALSE,
  parallel = "multicore", ncpus = 3))
  
confint(boothr2, type = "perc")

save("boothr2",
  file = paste0("hr2_bootstrap", round(nboot/1000, 1), "k.rdata"))
Sys.time()
} else{  #<-- end turning "off" bootstrap
    load("hr2_bootstrap10k.rdata")
  }



# Have a peek at bootstrap distributions
## Change `Sys.sleep()` to pause more between plots
sapply(seq(20), FUN = function(x){ plot(boothr2, index = x);
    mtext(dimnames(boothr2$t)[[2L]][x], side = 3, outer = TRUE, line = -2);
    Sys.sleep(0.002); invisible(x)})

dev.off()

hr2Ests <- cbind(Est = parmOut2(hr2_fullFxd_Trt_REML),
  confint(boothr2, type = "perc"))


# Make data.frame for plotting predicted values (with CIs)
## And for pulling out specific predictions and CIs
cidata1[, c("pred", "lCL", "uCL")] <- hr1Ests[24:323, ]
  cidata1[, c("pred_unsc", "lCL_unsc", "uCL_unsc")] <- hr1Ests[324:623, ]
cidata2[, c("pred", "lCL", "uCL")] <- hr2Ests[24:323, ]
  cidata2[, c("pred_unsc", "lCL_unsc", "uCL_unsc")] <- hr2Ests[324:623, ]






#######################################################
###### 		Manuscript summary table of effects
#######################################################
round(hr1Ests[1:23, ], 4)   #<-- subset to drop the model prediction rows
round(hr2Ests[1:23, ], 4)
#XXX EXPLANATION: XXX
  # (Intercept) is where control treatment among-individual (T1wsmn) line meets y-axis at Time=0
  # Trtlow and Trtper are the deviations from this where their among-individual treatment mean lines also intersect y-axis (Time=0)

  # "Trtcont:T1wsdev" etc. at the bottom are then the within-subject mean slopes for each Treatment on un-scaled Heart Rate (slopes, so change in un-scaled HR over Time in minutes)
  
  
# XXX predicted scHR and HR at the minimum time (~0.33 minutes) where the treatment within-subject slopes Time=0.33 at the mean of the within-subject means. So this is the average heart rate for the average individual in a treatment when Time~0.33
cidata1[minTind, ]
cidata2[minTind, ]


#(Alternatively, use emmeans to predict treatment averages at minimum time)
## Compare to last bit (just `minTind` entries)
rg1 <- ref_grid(hr1_fullFxd_Trt_REML,
  at = list(T1wsmn = mean1WSmn, T1wsdev = (minT1 - mean1WSmn)))
emmeans(rg1, specs = ~ Trt + T1wsmn:Trt + T1wsdev:Trt)

rg2 <- ref_grid(hr2_fullFxd_Trt_REML,
  at = list(T2wsmn = mean2WSmn, T2wsdev = (minT2 - mean2WSmn)))
emmeans(rg2, specs = ~ Trt + T2wsmn:Trt + T2wsdev:Trt)


# Type III ANOVA results
anova(hr1_fullFxd_Trt_REML)
anova(hr2_fullFxd_Trt_REML)



# Model Estimate summaries
summary(hr1_fullFxd_Trt_REML)
  logLik(hr1_fullFxd_Trt_REML, REML = TRUE)
summary(hr2_fullFxd_Trt_REML)
  logLik(hr2_fullFxd_Trt_REML, REML = TRUE)










#######################################################
#######    Hypothesis tests / Likelihood ratio tests
#######################################################

## 30% Development
Model1_hr1 <- hr1_fullFxd_Trt_REML
Model2_hr1 <- hr1_fullFxd_TrtVintVslo_REML
Model3_hr1 <- hr1_fullFxd_TrtVintHomVslo_REML

logLik(Model1_hr1, REML = TRUE)
  nobs(Model1_hr1)
logLik(Model2_hr1, REML = TRUE)
  nobs(Model2_hr1)
logLik(Model3_hr1, REML = TRUE)
  nobs(Model3_hr1)


# Significant correlations between treatment-specific intercept and slopes?
anova(Model1_hr1, Model2_hr1, refit = FALSE)
# Significant treatment-specific random slopes (UNcorrelated random int-slopes)?
anova(Model1_hr1, Model3_hr1, refit = FALSE)
  

## 80% Development
Model1_hr2 <- hr2_fullFxd_Trt_REML
Model2_hr2 <- hr2_fullFxd_TrtVintVslo_REML
Model3_hr2 <- hr2_fullFxd_TrtVintHomVslo_REML

logLik(Model1_hr2, REML = TRUE)
  nobs(Model1_hr2)
logLik(Model2_hr2, REML = TRUE)
  nobs(Model2_hr2)
logLik(Model3_hr2, REML = TRUE)
  nobs(Model3_hr2)


# Significant correlations between treatment-specific intercept and slopes?
anova(Model1_hr2, Model2_hr2, refit = FALSE)
# Significant treatment-specific random slopes (UNcorrelated random int-slopes)?
anova(Model1_hr2, Model3_hr2, refit = FALSE)
 














################################################################################

#         FIGURES                 ##############################################

################################################################################



#######################################
# Fixed effects model predictions plot
# (random effects not included)
#######################################
pdf(file = "Fig3_AverageWithinSubject.pdf",
  width = 12, height = 5)
#x11(w = 12, h = 5)
par(mfrow = c(1, 2), mar = c(6, 6, 5, 1.1), cex.lab = 1.3, cex.axis = 1.2)
###############################################################
# Treatment-mean predicted value
###############################################################
#####   30% Development   ###################
Tnumb <- 1
Tcov <- paste0("T", Tnumb)
pForm <- as.formula(paste0("scHR", Tnumb, "~ ", Tcov))

plot(pForm, data = hrd1, type = "n", axes = FALSE,
    xlim = c(0, 1.2),
    ylim = c(-1.50, 1.5),
    main = "",
    xlab = "Time (min.)",
    ylab = "30% Development\nHeart rate (standardized bpm)")
  axis(1)
  axis(2)

  # Confidence Intervals FIRST, so they are "behind" the prediction lines
  ## Go through each treatment
  cnt <- 1
  for(t in levels(hrd1$Trt)){
    ## Confidence Interval limits
    lines(lCL ~ I(T1wsdev + T1wsmn), data = cidata1, subset = cidata1$Trt == t,
      col = clr[cnt], lwd = 3.5, lty = "dotted")
    ## Confidence Interval limits
    lines(uCL ~ I(T1wsdev + T1wsmn), data = cidata1, subset = cidata1$Trt == t,
      col = clr[cnt], lwd = 3.5, lty = "dotted")
    cnt <- cnt + 1
  }

  # Prediction  line
  ## Go through each treatment
  cnt <- 1
  for(t in levels(hrd1$Trt)){
    lines(pred ~ I(T1wsdev + T1wsmn), data = cidata1, subset = cidata1$Trt == t,
      col = clr[cnt], lwd = 5)
    cnt <- cnt + 1
  }
  
mtext(text = expression(bold(A)),
  side = 3, line = 2.5, adj = -0.15, cex = 1.5)
      
###########################

legend("topright", inset = c(-0.04, -0.3), xpd = TRUE,
  lwd = c(4, 4, 4, 3),
  lty = c(rep("solid", 3), "dotted"),
  col = c(clr, "grey30"),
  legend = c("Control", "Low", "Periodic", "95% CI"))
###########################


#####   80% Development   ###################
Tnumb <- 2
Tcov <- paste0("T", Tnumb)
pForm <- as.formula(paste0("scHR", Tnumb, "~ ", Tcov))

plot(pForm, data = hrd2, type = "n", axes = FALSE,
    xlim = c(0, 1.2),
    ylim = c(-1.50, 1.5),
    main = "",
    xlab = "Time (min.)",
    ylab = "80% Development\nHeart rate (standardized bpm)")
  axis(1)
  axis(2)

  # Confidence Intervals FIRST, so they are "behind" the prediction lines
  ## Go through each treatment
  cnt <- 1
  for(t in levels(hrd2$Trt)){
    ## Confidence Interval limits
    lines(lCL ~ I(T2wsdev + T2wsmn), data = cidata2, subset = cidata2$Trt == t,
      col = clr[cnt], lwd = 3.5, lty = "dotted")
    ## Confidence Interval limits
    lines(uCL ~ I(T2wsdev + T2wsmn), data = cidata2, subset = cidata2$Trt == t,
      col = clr[cnt], lwd = 3.5, lty = "dotted")
    cnt <- cnt + 1
  }

  # Prediction  line
  ## Go through each treatment
  cnt <- 1
  for(t in levels(hrd2$Trt)){
    lines(pred ~ I(T2wsdev + T2wsmn), data = cidata2, subset = cidata2$Trt == t,
      col = clr[cnt], lwd = 5)
    cnt <- cnt + 1
  }

mtext(text = expression(bold(B)),
  side = 3, line = 2.5, adj = -0.15, cex = 1.5)


dev.off() 
















###############################################################
# Individual predicted reaction norms
###############################################################
# Plot regression lines, by individual, to data

#####   30% Development   ###################
sdata <- hrd  #<-- define temporary data for plotting script
dev <- 30
sdata$pred1 <- predict(hr1_fullFxd_Trt_REML, re.form = NULL)


pdf(file = "Fig4_individRxnNorms_30dev.pdf",
  width = 12, height = 5)
#x11(w = 12, h = 5)
##########
par(mfrow = c(1, 3), mar = c(6, 6, 4, 1.1), cex.lab = 1.5, cex.axis = 1.4)
##########
Tnumb <- ifelse(dev == 30, 1, 2)
Tcov <- paste0("T", Tnumb)
pForm <- as.formula(paste0("scHR", Tnumb, "~ ", Tcov))

### Go through each individual within each treatment:
  cnt <- 1
  for(t in levels(sdata$Trt)){
    plot(pForm, data = sdata, type = "n", axes = FALSE,
      xlim = c(0, 1.2),
      ylim = c(-2.8, 3.45),
      main = "",
      xlab = "Time (min.)",
      ylab = "30% Development\nHeart rate (standardized bpm)")
    axis(1)
    axis(2)
    uid <- unique(as.character(sdata$EggID[which(sdata$Trt == t)]))

    
    for(i in uid){
      ### Subset data for i-th ID/individual
      sub_datf <- subset(sdata, EggID == i)
      #### Order subsetted data by increasing Time
      ##### Needed to plot lines
      sub_datf <- sub_datf[order(sub_datf[, Tcov]), ]
      ### Plot Individual observations
      points(x = sub_datf[, Tcov],
        y = sub_datf[, paste0("scHR", Tnumb)],
        col = clrg_trans[cnt], pch = 21, lwd = 2, cex = 0.8)
      ### Add predicted regression
      lines(x = sub_datf[, Tcov],
            y = sub_datf[, "pred1"],
            col = clrg[cnt], lwd = 1)
      cnt <- cnt + 1
    }  #<-- end for i in unique IDs

    if(t == "low"){
      legend("top",
        pch = c(21, NA), pt.lwd = 2,
        lwd = c(NA, 3),
        col = "grey80",
        legend = c("observed", "predicted"))
    } 
     
  }  #<-- end for t in Treatments
  mtext(text = expression(bold(A)),
    side = 3, line = 1, adj = -3.0, cex = 1.5)
  mtext(text = expression(bold(B)),
    side = 3, line = 1, adj = -1.6, cex = 1.5)
  mtext(text = expression(bold(C)),
    side = 3, line = 1, adj = -0.2, cex = 1.5)




dev.off() 



#####   80% Development   ###################
sdata <- hrd  #<-- define temporary data for plotting script
dev <- 80
sdata$pred2 <- predict(hr2_fullFxd_Trt_REML, re.form = NULL)


pdf(file = "Fig5_individRxnNorms_80dev.pdf",
  width = 12, height = 5)
#x11(w = 12, h = 5)
##########
par(mfrow = c(1, 3), mar = c(6, 6, 4, 1.1), cex.lab = 1.5, cex.axis = 1.4)
##########
Tnumb <- ifelse(dev == 30, 1, 2)
Tcov <- paste0("T", Tnumb)
pForm <- as.formula(paste0("scHR", Tnumb, "~ ", Tcov))

### Go through each individual within each treatment:
  cnt <- 1
  for(t in levels(sdata$Trt)){
    plot(pForm, data = sdata, type = "n", axes = FALSE,
      xlim = c(0, 1.2),
      ylim = c(-2.8, 3.45),
      main = "",
      xlab = "Time (min.)",
      ylab = "80% Development\nHeart rate (standardized bpm)")
    axis(1)
    axis(2)
    uid <- unique(as.character(sdata$EggID[which(sdata$Trt == t)]))

    
    for(i in uid){
      ### Subset data for i-th ID/individual
      sub_datf <- subset(sdata, EggID == i & !is.na(sdata$pred2))
      #### Order subsetted data by increasing Time
      ##### Needed to plot lines
      sub_datf <- sub_datf[order(sub_datf[, Tcov]), ]
      ### Plot Individual observations
      points(x = sub_datf[, Tcov],
        y = sub_datf[, paste0("scHR", Tnumb)],
        col = clrg_trans[cnt], pch = 21, lwd = 2, cex = 0.8)
      ### Add predicted regression
      lines(x = sub_datf[, Tcov],
            y = sub_datf[, "pred2"],
            col = clrg[cnt], lwd = 1)
      cnt <- cnt + 1
    }  #<-- end for i in unique IDs

    if(t == "low"){
      legend("top",
        pch = c(21, NA), pt.lwd = 2,
        lwd = c(NA, 3),
        col = "grey80",
        legend = c("observed", "predicted"))
    } 
     
  }  #<-- end for t in Treatments
  mtext(text = expression(bold(A)),
    side = 3, line = 1, adj = -3.0, cex = 1.5)
  mtext(text = expression(bold(B)),
    side = 3, line = 1, adj = -1.6, cex = 1.5)
  mtext(text = expression(bold(C)),
    side = 3, line = 1, adj = -0.2, cex = 1.5)




dev.off() 












################################################################################
################################################################################
#####       ######          ####################################################
####   ###  ##########  ########################################################
###    ###############  ########################################################
####   ###############  ########################################################
######     ###########  ########################################################
##########   #########  ########################################################
##########    ########  ########################################################
####   ###   #########  ########################################################
#####       ######          ####################################################
################################################################################
################################################################################



################################################################################

#			SI 1

# DESCRIPTIVE STATISTICS

################################################################################


# TABLE SI 1
# Range, min, mean, median, and max of number measures per individual
msreCntsByEggID1 <- aggregate(HR1 ~ Trt + EggID, data = hrd1, FUN = length)
msreCntsByEggID2 <- aggregate(HR2 ~ Trt + EggID, data = hrd2, FUN = length)

## Summary for all embryos (not by Dev or treatment)
### presented in main results section
mean(c(msreCntsByEggID1[, 3], msreCntsByEggID2[, 3]))
median(c(msreCntsByEggID1[, 3], msreCntsByEggID2[, 3]))

## Now pull out summary statistics by Development stage and Treatment
### 30%
cbind(aggregate(HR1 ~ Trt, data = msreCntsByEggID1, FUN = length),
  min = aggregate(HR1 ~ Trt, data = msreCntsByEggID1, FUN = min)[, 2],
  mean = aggregate(HR1 ~ Trt, data = msreCntsByEggID1, FUN = mean)[, 2],
  median = aggregate(HR1 ~ Trt, data = msreCntsByEggID1, FUN = median)[, 2],
  max = aggregate(HR1 ~ Trt, data = msreCntsByEggID1, FUN = max)[, 2])
### 80%
cbind(aggregate(HR2 ~ Trt, data = msreCntsByEggID2, FUN = length),
  min = aggregate(HR2 ~ Trt, data = msreCntsByEggID2, FUN = min)[, 2],
  mean = aggregate(HR2 ~ Trt, data = msreCntsByEggID2, FUN = mean)[, 2],
  median = aggregate(HR2 ~ Trt, data = msreCntsByEggID2, FUN = median)[, 2],
  max = aggregate(HR2 ~ Trt, data = msreCntsByEggID2, FUN = max)[, 2])












################################################################################

#			SI 2

# Splines to assess suitability of LINEAR REACTION NORMS

################################################################################


# Plot splines, by individual, to data to evalute use of a linear reaction norm
## 30% Development
sdata <- hrd1  #<-- define temporary "smooth data" for plotting script
dev <- 30



pdf(file = "FigSI2-1_SmoothedSplines_30dev.pdf",
  width = 12, height = 5)
#x11(w = 12, h = 5)
##########
par(mfrow = c(1, 3), mar = c(6, 6, 4, 1.1), cex.lab = 1.5, cex.axis = 1.4)
##########
Tnumb <- ifelse(dev == 30, 1, 2)
Tcov <- paste0("T", Tnumb)
pForm <- as.formula(paste0("scHR", Tnumb, "~ ", Tcov))

### Go through each individual within each treatment:
  cnt <- 1
  for(t in levels(sdata$Trt)){
    plot(pForm, data = sdata, type = "n", axes = FALSE,
      xlim = c(0, 1.2),
      ylim = c(-2.8, 3.45),
      main = "",
      xlab = "Time (min.)",
      ylab = "30% Development\nHeart rate (standardized bpm)")
    axis(1)
    axis(2)
    uid <- unique(as.character(sdata$EggID[which(sdata$Trt == t)]))
    for(i in uid){
      ### Subset data for i-th ID/individual
      sub_datf <- subset(sdata, EggID == i)
      #### Order subsetted data by increasing Time
      ##### Needed to plot lines
      sub_datf <- sub_datf[order(sub_datf[, Tcov]), ]
      ### Plot Individual observations
      points(x = sub_datf[, Tcov],
        y = sub_datf[, paste0("scHR", Tnumb)],
        col = clrg_trans[cnt], pch = 21, lwd = 2, cex = 0.8)
      ### Add smoothed spline
      #### (only if have >3 measurement)
      if(nrow(sub_datf) > 3){
        lines(smooth.spline(x = sub_datf[, Tcov],
            y = sub_datf[, paste0("scHR", Tnumb)],
            cv = TRUE),  #<-- uses generalized cross-validation
          col = clrg[cnt], lwd = 1)
      }
      cnt <- cnt + 1
    }  #<-- end for i in unique IDs

    if(t == "low"){
      legend("top",
        pch = c(21, NA), pt.lwd = 2,
        lwd = c(NA, 3),
        col = "grey80",
        legend = c("observed", "smoothed spline"))
    } 
     
  }  #<-- end for t in Treatments
  mtext(text = expression(bold(A)),
    side = 3, line = 1, adj = -3.0, cex = 1.5)
  mtext(text = expression(bold(B)),
    side = 3, line = 1, adj = -1.6, cex = 1.5)
  mtext(text = expression(bold(C)),
    side = 3, line = 1, adj = -0.2, cex = 1.5)




dev.off() 









## 80% Development
sdata <- hrd2  #<-- define temporary "smooth data" for plotting script
dev <- 80




pdf(file = "FigSI2-2_SmoothedSplines_80dev.pdf",
  width = 12, height = 5)
#x11(w = 12, h = 5)
##########
par(mfrow = c(1, 3), mar = c(6, 6, 4, 1.1), cex.lab = 1.5, cex.axis = 1.4)
##########
Tnumb <- ifelse(dev == 30, 1, 2)
Tcov <- paste0("T", Tnumb)
pForm <- as.formula(paste0("scHR", Tnumb, "~ ", Tcov))

### Go through each individual within each treatment:
  cnt <- 1
  for(t in levels(sdata$Trt)){
    plot(pForm, data = sdata, type = "n", axes = FALSE,
      xlim = c(0, 1.2),
      ylim = c(-2.8, 3.45),
      main = "",
      xlab = "Time (min.)",
      ylab = "80% Development\nHeart rate (standardized bpm)")
    axis(1)
    axis(2)
    uid <- unique(as.character(sdata$EggID[which(sdata$Trt == t)]))
    for(i in uid){
      ### Subset data for i-th ID/individual
      sub_datf <- subset(sdata, EggID == i)
      #### Order subsetted data by increasing Time
      ##### Needed to plot lines
      sub_datf <- sub_datf[order(sub_datf[, Tcov]), ]
      ### Plot Individual observations
      points(x = sub_datf[, Tcov],
        y = sub_datf[, paste0("scHR", Tnumb)],
        col = clrg_trans[cnt], pch = 21, lwd = 2, cex = 0.8)
      ### Add smoothed spline
      #### (only if have >3 measurement)
      if(nrow(sub_datf) > 3){
        lines(smooth.spline(x = sub_datf[, Tcov],
            y = sub_datf[, paste0("scHR", Tnumb)],
            cv = FALSE),  #<-- uses generalized cross-validation
          col = clrg[cnt], lwd = 1)
      }
      cnt <- cnt + 1
    }  #<-- end for i in unique IDs

    if(t == "low"){
      legend("top",
        pch = c(21, NA), pt.lwd = 2,
        lwd = c(NA, 3),
        col = "grey80",
        legend = c("observed", "smoothed spline"))
    } 
     
  }  #<-- end for t in Treatments
  mtext(text = expression(bold(A)),
    side = 3, line = 1, adj = -3.0, cex = 1.5)
  mtext(text = expression(bold(B)),
    side = 3, line = 1, adj = -1.6, cex = 1.5)
  mtext(text = expression(bold(C)),
    side = 3, line = 1, adj = -0.2, cex = 1.5)




dev.off() 














################################################################################

#			SI 4

# ORDER AND MASS EFFECTS

################################################################################

# Test that order and mass had no effects:
### Use Type III sums of squares and Satterthwaite's method
#XXX NOTE ued "LayMass" here
hr1_fullFxd_ordMass_Trt_REML <- lmer(scHR1 ~ Trt + T1wsmn:Trt + T1wsdev:Trt +
  Order1fac + LayMass +
  (dummy(Trt, "cont"):T1 | EggID) + 
  (dummy(Trt, "low"):T1 | EggID) +
  (dummy(Trt, "per"):T1 | EggID) + 
  (1 | NestID),
	data = hrd, REML = TRUE, na.action = na.exclude,
	control = lmerControl(optimizer="optimx",
        	optCtrl = list(method="L-BFGS-B")))

summary(hr1_fullFxd_ordMass_Trt_REML)
# Table S4.1
anova(hr1_fullFxd_ordMass_Trt_REML)
  

# Test that order and mass had no effects:
### Use Type III sums of squares and Satterthwaite's method
#XXX NOTE used "Mass2" here
hr2_fullFxd_ordMass_Trt_REML <- lmer(scHR2 ~ Trt + T2wsmn:Trt + T2wsdev:Trt +
  Order2fac + Mass2 +
  (dummy(Trt, "cont"):T2 | EggID) + 
  (dummy(Trt, "low"):T2 | EggID) +
  (dummy(Trt, "per"):T2 | EggID) + 
  (1 | NestID),
	data = hrd, REML = TRUE, na.action = na.exclude,
#  	control = lmerControl(optimizer = "nlminbwrap"))
	control = lmerControl(optimizer="optimx",
        	optCtrl = list(method="L-BFGS-B")))

summary(hr2_fullFxd_ordMass_Trt_REML)
# Table S4.2
anova(hr2_fullFxd_ordMass_Trt_REML)
  







################################################################################

#			SI 5

# MODEL RESIDUAL INSPECTION

################################################################################
## Grab standardized residuals and put them in the data
hrd$stdRes1 <- resid(hr1_fullFxd_Trt_REML, scaled = TRUE)
hrd$stdRes2 <- resid(hr2_fullFxd_Trt_REML, scaled = TRUE)
## Now same for fitted values
hrd$ftd1 <- fitted(hr1_fullFxd_Trt_REML)
hrd$ftd2 <- fitted(hr2_fullFxd_Trt_REML)
# Now give residuals their own object and REMOVE NAs
stdResHR1 <- na.omit(hrd$stdRes1)
stdResHR2 <- na.omit(hrd$stdRes2)


# Q-Q plot of standardized residuals
pdf(file = "FigSI5-1_QQplot_StdizdResiduals.pdf",
  width = 8, height = 5)
#x11(w = 8, h = 5)
##########
par(mfrow = c(1, 2), mar = c(5, 6, 4, 1.1), cex.lab = 1.5, cex.axis = 1.4)
 qqnorm(stdResHR1, main = "Normal Q-Q Plot\n30% Development",
   lwd = 2, col = clr[as.integer(hrd$Trt)])
  qqline(stdResHR1)
 qqnorm(stdResHR2, main = "Normal Q-Q Plot\n80% Development",
   lwd = 2, col = clr[as.integer(hrd$Trt)])
  qqline(stdResHR2)



dev.off() 


# Standardized residuals versus fitted values by Treatment
mains <- c("Control", "Low", "Periodic")
tlvls <- levels(hrd$Trt)

pdf(file = "FigSI5-2_StdizdResiduals-VS-Fitted.pdf",
  width = 12, height = 10)
#x11(w = 12, h = 10)
par(mfrow = c(2, 3), mar = c(6, 6, 7, 1),
  cex.main = 2, cex.lab = 1.9, cex.axis = 1.4)
for(d in 1:2){
  tform <- as.formula(paste0("stdRes", d, " ~ ftd", d))
  for(t in tlvls){
    plot(tform, data = hrd, subset = hrd$Trt == t, type = "n",
      xlab = "Fitted values", ylab = "Standardized residuals",
      main = paste0(ifelse(d==1, "3", "8"), "0% Development", "\n",
        mains[which(tlvls == t)], "\n"))
      abline(h = 0, lwd = 2, col = "grey70")
      points(tform, data = hrd, subset = hrd$Trt == t,
        pch = 21, col = clr[1], lwd = 1)
  }    
} 



dev.off() 









# Standardized residuals versus TIME by Treatment
pdf(file = "FigSI5-3_StdizdResiduals-VS-Time.pdf",
  width = 12, height = 10)
#x11(w = 12, h = 10)
par(mfrow = c(2, 3), mar = c(6, 6, 7, 1),
  cex.main = 2, cex.lab = 1.9, cex.axis = 1.4)
for(d in 1:2){
  tform <- as.formula(paste0("stdRes", d, " ~ T", d))
  for(t in tlvls){
    plot(tform, data = hrd, subset = hrd$Trt == t, type = "n",
      xlab = "Time", ylab = "Standardized residuals",
      main = paste0(ifelse(d==1, "3", "8"), "0% Development", "\n",
        mains[which(tlvls == t)], "\n"))
      abline(h = 0, lwd = 2, col = "grey70")
      points(tform, data = hrd, subset = hrd$Trt == t,
        pch = 21, col = clr[1], lwd = 1)
      # add a smoothed spline
      lines(smooth.spline(x = hrd[which(hrd$Trt == t &
            !is.na(hrd[, paste0("T", d)])), paste0("T", d)],
          y = hrd[which(hrd$Trt == t &
            !is.na(hrd[, paste0("T", d)])), paste0("stdRes", d)],
          cv = TRUE),  #<-- TRUE leave-one-out / FALSE generalized cross-validation
        col = clr[3], lwd = 1)
  }    
} 



dev.off() 




# Durbin-Watson test for Autocorrelation in the residuals
## Need to re-fit model to data without NAs for DHARMa functions
lmm1b <- update(lmm1, data = hrd1)
res1 <- simulateResiduals(lmm1b)
  res1b <- recalculateResiduals(res1, group = hrd1$T1)
lmm2b <- update(lmm2, data = hrd2)
res2 <- simulateResiduals(lmm2b)
  res2b <- recalculateResiduals(res2, group = hrd2$T2)


# 30% Development
pdf(file = "FigSI5-4_Autocorr_StdizdResiduals_30dev.pdf",
  width = 8, height = 5)
#x11(w = 8, h = 5)
par(mar = c(5, 6, 5, 1),
  cex.main = 2, cex.lab = 1.9, cex.axis = 1.4)
testTemporalAutocorrelation(res1b, tim = unique(hrd1$T1))

dev.off()


# 80% Development
pdf(file = "FigSI5-5_Autocorr_StdizdResiduals_80dev.pdf",
  width = 8, height = 5)
#x11(w = 8, h = 5)
par(mar = c(5, 6, 5, 1),
  cex.main = 2, cex.lab = 1.9, cex.axis = 1.4)
testTemporalAutocorrelation(res2b, tim = unique(hrd2$T2))

dev.off()














#  				FIN
