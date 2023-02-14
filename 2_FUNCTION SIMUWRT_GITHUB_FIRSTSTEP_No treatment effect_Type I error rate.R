################################################################################################
# Comparison of statistical methods to analyze a within-person randomized trial 
# evaluating a topical drug with a binary outcome: a simulation study

# by Sophie Leducq
# June 2022
################################################################################################

################################################################################################
# Summary of methods used for analysing of each simulated dataset  
# Included in first step (type I error rate)
################################################################################################

################################################################################################
# Packages
################################################################################################

library (MultiOrd)
library (joineR)
library (plyr)
library (dplyr)
library (lme4)
library (ICC)
library (geepack)
library (stats)
library (descr)
library (MASS)
library (nlme)
library (geeM)
library (lm.beta)
library (broom)
library(broom.mixed)
library(lmerTest)
library (sqldf)
library(tidyr)
library (gee)
library (glmmTMB)


################################################################################################
# INPUT
# K: Total number of patients/stratum (in the 2 groups)
# P_0: Proportion in control group
# P_1: Proportion in treatment group
# ALPHA: Alpha risk
# rho: Within-patient intragroup correlation
# coefeta: coefeta*rho = Within-patient intergroup correlation
# m: Mean number of lesions per patient/stratum
# maxm: Maximum number of lesions per patient/stratum (The distribution is truncated to have a 
# maximum number of 10, 12 and 18 lesions for m equal to 4, 6 and 10, respectively) 
################################################################################################


simulation_WRT <-  function(i, K, P_0, P_1, ALPHA, rho, coefeta, m, maxm) {
  
  
  ################################################################################################
  # Generation of number of lesions per patient (ie stratum) using a Poisson distribution 
  # with parameter m
  # nk: number of lesions per patient
  # nk must not exceeed maxm
  ################################################################################################
  
  #balanced and unbalanced ratio between the 2 groups
  
  nk <- rpois (K, lambda= (m)-2)+2 #The distribution is truncated to have a minimum number of two lesions 
  nk
  
  for (i in 1:K) if(nk[i]>maxm) nk[i] <- maxm #nk must not exceeed maxv, 
  #The distribution is truncated to have a maximum number of 10, 12 and 18 lesions for m equal to 4, 6 and 10, respectively
  
  dataframe <- data.frame(SIZESTRATUM= nk)
  dataframe$IMPAIR <- 1 #I create an uneven variable (named 'IMPAIR') to identify uneven stratum
  for (i in 1:K) {
    if(dataframe$SIZESTRATUM[i] %% 2 == 0) dataframe$IMPAIR[i] <- 0
    else dataframe$IMPAIR[i]<- 1
  }
  #if the number of lesions generated is an odd number, we generate one more outcome and then randomly discarded one of these outcomes after
  for (i in 1:K) {
    if(dataframe$IMPAIR[i]== 1) dataframe$SIZESTRATUM[i] <- (dataframe$SIZESTRATUM[i]+1)
    else dataframe$SIZESTRATUM[i]<- dataframe$SIZESTRATUM[i]
  }
  
  dataframe$SIZESTRATUM<- dataframe$SIZESTRATUM/2 #1 arm/group here
  
  dataframe$SUBJID <- 1:K
  
  dataframeE <- uncount(dataframe, SIZESTRATUM) #to get my rows by stratum and keep my variable uneven
  dataframeC <- uncount(dataframe, SIZESTRATUM)
  
  dataframeE$TTT <- 1
  dataframeC$TTT <- 2 #add the variable treatment
  
  data <- rbind (dataframeE, dataframeC)
  datafinal <- data[order(data$SUBJID), ]
  #I have a data with one observation per line with variability in stratum size 
  #SUBJID = stratum = patient
  #TTT = intervention or control treatment
  
  ################################################################################################
  # Definition of eta
  ################################################################################################
  
  eta <- rho * coefeta
  
  ################################################################################################
  # Creation of correlation matrix regarding number of lesions per stratum
  # Division of datafinal in sub-table
  ################################################################################################
  
  datafinal$somme <- seq (1,1,1)
  for (i in 1:K){
    datafinal$sizestratum[datafinal$SUBJID==i]<- sum(datafinal$somme[data$SUBJID==i])
  } #creation of a variable that specify the stratum size /number of lesions per patient 
  
  
  datastratum2 <- datafinal [is.element (datafinal$sizestratum, 2),]
  datastratum4 <- datafinal [is.element (datafinal$sizestratum, 4),]
  datastratum6 <- datafinal [is.element (datafinal$sizestratum, 6),]
  datastratum8 <- datafinal [is.element (datafinal$sizestratum, 8),]
  datastratum10 <- datafinal [is.element (datafinal$sizestratum, 10),]
  datastratum12 <- datafinal [is.element (datafinal$sizestratum, 12),]
  datastratum14 <- datafinal [is.element (datafinal$sizestratum, 14),]
  datastratum16 <- datafinal [is.element (datafinal$sizestratum, 16),]
  datastratum18 <- datafinal [is.element (datafinal$sizestratum, 18),]
  
  #correlation matrix in case of 2 lesions per stratum
  corr.mat_2 <- toeplitz(c(1,eta)) 
  
  #correlation matrix in case of 4 lesions per stratum
  matrice1 <- toeplitz(c(1,rho))
  matrice2 <- toeplitz(c(eta,eta))
  matrice3 <- toeplitz(c(eta,eta))
  matrice4 <- toeplitz(c(1,rho))
  matricea <- rbind (matrice1, matrice2)
  matriceb <- rbind (matrice3, matrice4)
  corr.mat_4 <- cbind (matricea, matriceb)
  
  #correlation matrix in case of 6 lesions per stratum
  matrice1 <- toeplitz(c(1,rho, rho))
  matrice2 <- toeplitz(c(eta,eta,eta))
  matrice3 <- toeplitz(c(eta,eta,eta))
  matrice4 <- toeplitz(c(1,rho, rho))
  matricea <- rbind (matrice1, matrice2)
  matriceb <- rbind (matrice3, matrice4)
  corr.mat_6 <- cbind (matricea, matriceb)
  
  #correlation matrix in case of 8 lesions per stratum
  matrice1 <- toeplitz(c(1,rho, rho, rho))
  matrice2 <- toeplitz(c(eta,eta,eta, eta))
  matrice3 <- toeplitz(c(eta,eta,eta, eta))
  matrice4 <- toeplitz(c(1,rho, rho,rho))
  matricea <- rbind (matrice1, matrice2)
  matriceb <- rbind (matrice3, matrice4)
  corr.mat_8 <- cbind (matricea, matriceb)
  
  #correlation matrix in case of 10 lesions per stratum
  matrice1 <- toeplitz(c(1,rho, rho, rho, rho))
  matrice2 <- toeplitz(c(eta,eta,eta, eta, eta))
  matrice3 <- toeplitz(c(eta,eta,eta, eta, eta))
  matrice4 <- toeplitz(c(1,rho, rho,rho,rho))
  matricea <- rbind (matrice1, matrice2)
  matriceb <- rbind (matrice3, matrice4)
  corr.mat_10 <- cbind (matricea, matriceb)
  
  #correlation matrix in case of 12 lesions per stratum
  matrice1 <- toeplitz(c(1,rho, rho, rho, rho,rho))
  matrice2 <- toeplitz(c(eta,eta,eta, eta, eta, eta))
  matrice3 <- toeplitz(c(eta,eta,eta, eta, eta, eta))
  matrice4 <- toeplitz(c(1,rho, rho,rho, rho, rho))
  matricea <- rbind (matrice1, matrice2)
  matriceb <- rbind (matrice3, matrice4)
  corr.mat_12 <- cbind (matricea, matriceb)
  
  #correlation matrix in case of 14 lesions per stratum
  matrice1 <- toeplitz(c(1,rho, rho, rho, rho,rho, rho))
  matrice2 <- toeplitz(c(eta,eta,eta, eta, eta, eta, eta))
  matrice3 <- toeplitz(c(eta,eta,eta, eta, eta, eta, eta))
  matrice4 <- toeplitz(c(1,rho, rho,rho, rho, rho, rho))
  matricea <- rbind (matrice1, matrice2)
  matriceb <- rbind (matrice3, matrice4)
  corr.mat_14 <- cbind (matricea, matriceb)
  
  #correlation matrix in case of 16 lesions per stratum
  matrice1 <- toeplitz(c(1,rho, rho, rho, rho,rho,rho,rho))
  matrice2 <- toeplitz(c(eta,eta,eta, eta, eta, eta, eta, eta))
  matrice3 <- toeplitz(c(eta,eta,eta, eta, eta, eta, eta, eta))
  matrice4 <- toeplitz(c(1,rho, rho,rho, rho, rho,rho,rho))
  matricea <- rbind (matrice1, matrice2)
  matriceb <- rbind (matrice3, matrice4)
  corr.mat_16 <- cbind (matricea, matriceb)
  
  #correlation matrix in case of 18 lesions per stratum
  matrice1 <- toeplitz(c(1,rho, rho, rho, rho,rho,rho,rho,rho))
  matrice2 <- toeplitz(c(eta,eta,eta, eta, eta, eta,eta, eta, eta))
  matrice3 <- toeplitz(c(eta,eta,eta, eta, eta, eta,eta, eta, eta))
  matrice4 <- toeplitz(c(1,rho, rho,rho, rho, rho,rho,rho,rho ))
  matricea <- rbind (matrice1, matrice2)
  matriceb <- rbind (matrice3, matrice4)
  corr.mat_18 <- cbind (matricea, matriceb)
  
  
  
  ################################################################################################
  # Definition of proportion in the 2 groups
  ################################################################################################
  
  prop.vec.bin_2 <- c(rep(P_0,1), rep(P_1,1))
  
  
  ################################################################################################
  # Generation of binary correlated data
  # Outcomes are generated considering multivariate correlated binary data with marginal probabilities
  # following the Emrich and Piedmonte's method
  ################################################################################################
  
  #for stratum size = 2
  nObs <- (length(datastratum2$SUBJID)/2)
  if (nObs>0) { 
    data_2 <- generate.binary(nObs , prop.vec.bin_2, corr.mat_2)
    data_2 <- as.data.frame(data_2)
    
    colnames (data_2)<- c("1", "2")
    data_2$SUBJID <- 1:((length(datastratum2$SUBJID)/2)) #add SUBJID
    
    data_2 <-  to.unbalanced(data_2, id.col = 3, times = 1:2, Y.col = 1:2)#Transposition of the data to have one line per value/several lines per patient
    
    names (data_2)[3]<- "Y"
    
    
    data_2$TRAIT <- rep (c(0,1),(length(datastratum2$SUBJID)/2)) #Creation of variable treatment
    
    data_2$IMPAIR <- datastratum2$IMPAIR
    
  }else {data_2 <- data.frame(NULL)}
  
  
  #for stratum size = 4
  prop.vec.bin_4 <- c(rep(P_0,2), rep(P_1,2))
  
  nObs <- (length(datastratum4$SUBJID)/4)
  if (nObs>0) {
    data_4 <- generate.binary(nObs , prop.vec.bin_4, corr.mat_4)
    data_4 <- as.data.frame(data_4)
    
    colnames (data_4)<- c("1", "2","3", "4")
    length1 <- (length(datastratum2$SUBJID)/2)
    length2 <- ((length(datastratum2$SUBJID)/2)+(length(datastratum4$SUBJID)/4))
    data_4$SUBJID <- seq((length1 +1) , length2, 1 )
    
    data_4 <-  to.unbalanced(data_4, id.col = 5, times = 1:4, Y.col = 1:4)
    
    names (data_4)[3]<- "Y"
    
    data_4$TRAIT <- rep (c(0,0,1,1),(length(datastratum4$SUBJID)/4))
    
    data_4$IMPAIR <- datastratum4$IMPAIR
    
  }else {data_4 <- data.frame(NULL)}
  
  
  #for stratum size = 6
  prop.vec.bin_6 <- c(rep(P_0,3), rep(P_1,3))
  
  nObs <- (length(datastratum6$SUBJID)/6)
  if (nObs>0) {
    data_6 <- generate.binary(nObs , prop.vec.bin_6, corr.mat_6)
    data_6 <- as.data.frame(data_6)
    
    colnames (data_6)<- c("1", "2","3", "4", "5", "6")
    length1 <- (length(datastratum4$SUBJID)/4+(length(datastratum2$SUBJID)/2))
    length2 <- ((length(datastratum2$SUBJID)/2)+(length(datastratum4$SUBJID)/4)+(length(datastratum6$SUBJID)/6))
    data_6$SUBJID <- seq((length1 +1) , length2, 1 )
    
    data_6 <-  to.unbalanced(data_6, id.col = 7, times = 1:6, Y.col = 1:6)
    
    names (data_6)[3]<- "Y"
    
    data_6$TRAIT <- rep (c(0,0,0,1,1,1),(length(datastratum6$SUBJID)/6))
    
    data_6$IMPAIR <- datastratum6$IMPAIR
    
  }else {data_6 <- data.frame(NULL)}
  
  
  #for stratum size = 8
  prop.vec.bin_8 <- c(rep(P_0,4), rep(P_1,4))
  
  nObs <- (length(datastratum8$SUBJID)/8)
  if (nObs>0){
    data_8 <- generate.binary(nObs , prop.vec.bin_8, corr.mat_8)
    data_8 <- as.data.frame(data_8)
    
    colnames (data_8)<- c("1", "2","3", "4", "5", "6", "7", "8")
    length1 <- ((length(datastratum6$SUBJID)/6)+(length(datastratum4$SUBJID)/4)+(length(datastratum2$SUBJID)/2))
    length2 <- ((length(datastratum2$SUBJID)/2)+(length(datastratum4$SUBJID)/4)+(length(datastratum6$SUBJID)/6) +(length(datastratum8$SUBJID)/8) )
    data_8$SUBJID <- seq((length1 +1) , length2, 1 )
    
    data_8 <-  to.unbalanced(data_8, id.col = 9, times = 1:8, Y.col = 1:8)
    
    names (data_8)[3]<- "Y"
    
    data_8$TRAIT <- rep (c(0,0,0,0,1,1,1,1),(length(datastratum8$SUBJID)/8))
    
    data_8$IMPAIR <- datastratum8$IMPAIR
    
  }else {data_8 <- data.frame(NULL)}
  
  
  #for stratum size = 10
  prop.vec.bin_10 <- c(rep(P_0,5), rep(P_1,5))
  nObs <- (length(datastratum10$SUBJID)/10)
  if (nObs>0){
    data_10 <- generate.binary(nObs , prop.vec.bin_10, corr.mat_10)
    data_10 <- as.data.frame(data_10)
    
    
    colnames (data_10)<- c("1", "2","3", "4", "5", "6", "7", "8", "9", "10")
    length1 <-  ((length(datastratum2$SUBJID)/2)+(length(datastratum4$SUBJID)/4)+(length(datastratum6$SUBJID)/6) +(length(datastratum8$SUBJID)/8) )
    length2 <- ((length(datastratum2$SUBJID)/2)+(length(datastratum4$SUBJID)/4)+(length(datastratum6$SUBJID)/6) +(length(datastratum8$SUBJID)/8) + (length(datastratum10$SUBJID)/10))
    data_10$SUBJID <- seq((length1 +1) , length2, 1 )
    
    data_10 <-  to.unbalanced(data_10, id.col = 11, times = 1:10, Y.col = 1:10)
    
    names (data_10)[3]<- "Y"
    
    data_10$TRAIT <- rep (c(0,0,0,0,0,1,1,1,1,1),(length(datastratum10$SUBJID)/10))
    
    data_10$IMPAIR <- datastratum10$IMPAIR
    
  }else {data_10 <- data.frame(NULL)}
  
  
  #for stratum size = 12
  prop.vec.bin_12 <- c(rep(P_0,6), rep(P_1,6))
  nObs <- (length(datastratum12$SUBJID)/12)
  if (nObs>0){
    data_12 <- generate.binary(nObs , prop.vec.bin_12, corr.mat_12)
    data_12 <- as.data.frame(data_12)
    
    colnames (data_12)<- c("1", "2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
    length1 <-  ((length(datastratum2$SUBJID)/2)+(length(datastratum4$SUBJID)/4)+(length(datastratum6$SUBJID)/6) +(length(datastratum8$SUBJID)/8) + (length(datastratum10$SUBJID)/10))
    length2 <- ((length(datastratum2$SUBJID)/2)+(length(datastratum4$SUBJID)/4)+(length(datastratum6$SUBJID)/6) +(length(datastratum8$SUBJID)/8) + (length(datastratum10$SUBJID)/10) + (length(datastratum12$SUBJID)/12))
    data_12$SUBJID <- seq((length1 +1) , length2, 1 )
    
    data_12 <-  to.unbalanced(data_12, id.col = 13, times = 1:12, Y.col = 1:12)
    
    names (data_12)[3]<- "Y"
    
    data_12$TRAIT <- rep (c(0,0,0,0,0,0,1,1,1,1,1,1),(length(datastratum12$SUBJID)/12))
    
    data_12$IMPAIR <- datastratum12$IMPAIR
    
  }else {data_12 <- data.frame(NULL)}
  
  
  #for stratum size = 14
  prop.vec.bin_14 <- c(rep(P_0,7), rep(P_1,7))
  nObs <- (length(datastratum14$SUBJID)/14)
  if (nObs>0){
    data_14 <- generate.binary(nObs , prop.vec.bin_14, corr.mat_14)
    data_14 <- as.data.frame(data_14)
    
    colnames (data_14)<- c("1", "2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14")
    length1 <-  ((length(datastratum2$SUBJID)/2)+(length(datastratum4$SUBJID)/4)+(length(datastratum6$SUBJID)/6) +(length(datastratum8$SUBJID)/8) + (length(datastratum10$SUBJID)/10)+ (length(datastratum12$SUBJID)/12))
    length2 <- ((length(datastratum2$SUBJID)/2)+(length(datastratum4$SUBJID)/4)+(length(datastratum6$SUBJID)/6) +(length(datastratum8$SUBJID)/8) + (length(datastratum10$SUBJID)/10) + (length(datastratum12$SUBJID)/12)+ (length(datastratum14$SUBJID)/14))
    data_14$SUBJID <- seq((length1 +1) , length2, 1 )
    
    data_14 <-  to.unbalanced(data_14, id.col = 15, times = 1:14, Y.col = 1:14)
    
    names (data_14)[3]<- "Y"
    
    data_14$TRAIT <- rep (c(0,0,0,0,0,0,0,1,1,1,1,1,1,1),(length(datastratum14$SUBJID)/14))
    
    data_14$IMPAIR <- datastratum14$IMPAIR
    
  }else {data_14 <- data.frame(NULL)}
  
  
  #for stratum size = 16
  prop.vec.bin_16 <- c(rep(P_0,8), rep(P_1,8))
  nObs <- (length(datastratum16 $SUBJID)/16 )
  if (nObs>0){
    data_16  <- generate.binary(nObs , prop.vec.bin_16 , corr.mat_16 )
    data_16  <- as.data.frame(data_16 )
    
    colnames (data_16 )<- c("1", "2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16")
    length1 <-  ((length(datastratum2$SUBJID)/2)+(length(datastratum4$SUBJID)/4)+(length(datastratum6$SUBJID)/6) +(length(datastratum8$SUBJID)/8) + (length(datastratum10$SUBJID)/10)+ (length(datastratum12$SUBJID)/12)+ (length(datastratum14$SUBJID)/14))
    length2 <- ((length(datastratum2$SUBJID)/2)+(length(datastratum4$SUBJID)/4)+(length(datastratum6$SUBJID)/6) +(length(datastratum8$SUBJID)/8) + (length(datastratum10$SUBJID)/10) + (length(datastratum12$SUBJID)/12)+ (length(datastratum14$SUBJID)/14)+ (length(datastratum16$SUBJID)/16))
    data_16 $SUBJID <- seq((length1 +1) , length2, 1 )
    
    data_16 <-  to.unbalanced(data_16, id.col = 17, times = 1:16 , Y.col = 1:16 )
    
    names (data_16)[3]<- "Y"
    
    data_16$TRAIT <- rep (c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1),(length(datastratum16 $SUBJID)/16 ))
    
    data_16$IMPAIR <- datastratum16$IMPAIR
    
  }else {data_16 <- data.frame(NULL)}
  
  
  #for stratum size = 18
  prop.vec.bin_18 <- c(rep(P_0,9), rep(P_1,9))
  nObs <- (length(datastratum18 $SUBJID)/18 )
  if (nObs>0){
    data_18  <- generate.binary(nObs , prop.vec.bin_18 , corr.mat_18 )
    data_18  <- as.data.frame(data_18 )
    
    colnames (data_18 )<- c("1", "2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18")
    length1 <-  ((length(datastratum2$SUBJID)/2)+(length(datastratum4$SUBJID)/4)+(length(datastratum6$SUBJID)/6) +(length(datastratum8$SUBJID)/8) + (length(datastratum10$SUBJID)/10)+ (length(datastratum12$SUBJID)/12)+ (length(datastratum14$SUBJID)/14)+ (length(datastratum16$SUBJID)/16))
    length2 <- ((length(datastratum2$SUBJID)/2)+(length(datastratum4$SUBJID)/4)+(length(datastratum6$SUBJID)/6) +(length(datastratum8$SUBJID)/8) + (length(datastratum10$SUBJID)/10) + (length(datastratum12$SUBJID)/12)+ (length(datastratum14$SUBJID)/14)+ (length(datastratum16$SUBJID)/16) + (length(datastratum18$SUBJID)/18))
    data_18 $SUBJID <- seq((length1 +1) , length2, 1 )
    
    data_18 <-  to.unbalanced(data_18, id.col = 19, times = 1:18 , Y.col = 1:18 )
    
    names (data_18)[3]<- "Y"
    
    data_18$TRAIT <- rep (c(0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1),(length(datastratum18 $SUBJID)/18 ))
    
    data_18$IMPAIR <- datastratum18$IMPAIR
    
  }else {data_18 <- data.frame(NULL)}
  
 
  ################################################################################################
  # Gathering of data
  ################################################################################################
  base <- bind_rows (data_2, data_4, data_6, data_8, data_10, data_12, data_14, data_16, data_18)
  
  
  ################################################################################################
  # We randomly discard one of these outcomes for odd number
  ################################################################################################
  
  ligne <- c()
  for (i in unique(base[base$IMPAIR==1,]$SUBJID)) {
    ligne<- c(ligne, sample(as.numeric(row.names(base[base$SUBJID== i,])), 1))
  }
  base <- base[-ligne, ]
  
  
  ################################################################################################
  # Frequencies
  ################################################################################################
  
  freq1<- freq(base$Y[base$TRAIT==0], plot=F)
  P_0_freq <- freq1[2,2]
  
  freq2 <- freq(base$Y[base$TRAIT==1], plot=F)
  P_1_freq <- freq2[2,2]

  
  
  
  ################################################################################################
  # 1. ANALYSIS USING CONDITIONAL APPROACH BASED ON MIXED MODELS
  ################################################################################################
  
  ################################################################################################
  # 1.1 LOGIT LINK: ESTIMATION OF ODDS RATIO
  ################################################################################################
  
  ################################################################################################
  # 1.1.1 Correlation structure 1 (no model misspecification) 
  ################################################################################################
  
  fanalyse1 <- function(base, Y, TRAIT, SUBJID){
    
    
    # Fit regression model
    fit_analyse1 <- glmer(Y ~ TRAIT +  (1|SUBJID) + (1 |SUBJID:TRAIT),
                          data = base, family = binomial,control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000),
                                                                                check.conv.grad = .makeCC("stop", tol = 2e-3, relTol = NULL) ))
    
    # Pick up estimate
    
    #if singular fit  = T then subsitute fitS
    
    if( (isSingular (fit_analyse1, tol = 1e-4)==T)) {
      
      res_analyse1  <- data.frame(muest1="fitS", P_0_mod1="fitS", betaest1="fitS",
                              
                              P_1_mod1="fitS", betaest_se1="fitS",
                              
                              upperCI1="fitS", lowerCI1="fitS", pvalue1="fitS",
                              
                              varsubj1="fitS", varsubjttt1="fitS")
      
    } else {
    res_analyse1 <- data.frame(muest1 <- summary (fit_analyse1)$coefficients[1,1],
                               P_0_mod1 <-  ((exp(muest1))/(1+exp(muest1))),
                               
                               betaest1 <-  summary (fit_analyse1)$coefficients [2,1],
                               P_1_mod1 <-  ((exp (betaest1)*exp(muest1))/(1+(exp(betaest1)*exp(muest1)))),
                               
                               betaest_se1<-sqrt(diag(vcov(fit_analyse1)))[2],
                               
                               upperCI1<-betaest1 + 1.96*betaest_se1,
                               lowerCI1<-betaest1  - 1.96*betaest_se1,
                               
                               pvalue1<-summary (fit_analyse1)$coefficients[2,4],
                               
                               varsubj1<-as.numeric (VarCorr(fit_analyse1)[2]),
                               varsubjttt1<-as.numeric (VarCorr(fit_analyse1)[1]) )
    }
    # Return matrix containing estimate, SE, CI and p
    return(res_analyse1)
    
  }
  
  temp_analyse1 <- try(fanalyse1(base, Y, TRAIT, SUBJID), TRUE)
  

  
  
  # If converges, save output; otherwise substitute 9999s
  
  if(is.list(temp_analyse1)) {
    
    names(temp_analyse1)[1] <- "muest1"
    names(temp_analyse1)[2] <- "P_0_mod1"
    names(temp_analyse1)[3] <- "betaest1"
    names(temp_analyse1)[4] <-  "P_1_mod1"
    names(temp_analyse1)[5] <- "betaest_se1"
    names(temp_analyse1)[6] <- "upperCI1"
    names(temp_analyse1)[7] <- "lowerCI1"
    names(temp_analyse1)[8] <- "pvalue1"
    names(temp_analyse1)[9] <- "varsubj1"
    names(temp_analyse1)[10] <- "varsubjttt1"
    
    analyse1  <- as.data.frame(temp_analyse1)
    
  } else {
    
    analyse1  <- data.frame(muest1='NonC', P_0_mod1='NonC', betaest1='NonC',
                            
                            P_1_mod1='NonC', betaest_se1='NonC',
                            
                            upperCI1='NonC', lowerCI1='NonC', pvalue1='NonC',
                            
                            varsubj1='NonC', varsubjttt1='NonC')
    
  }
  
  
  
  ################################################################################################
  # 1.1.2 Correlation structure 3 (model misspecification) 
  ################################################################################################
  
  fanalyse12 <- function(base, Y, TRAIT, SUBJID){

    # Fit regression
    
    fanalyse12 <- glmer(Y ~ TRAIT +  (1|SUBJID) ,
                        data = base, family = binomial,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000),
                        check.conv.grad = .makeCC("stop", tol = 2e-3, relTol = NULL) ))
    
    # Pick up estimate
    
    #if singular fit  = T then subsitute fitS
    
    if( (isSingular (fanalyse12, tol = 1e-4)==T)) {
      
      fanalyse12  <- data.frame(muest12="fitS", P_0_mod12="fitS", betaest12="fitS",
                                   
                                   P_1_mod12="fitS", betaest_se12="fitS",
                                   
                                   upperCI12="fitS", lowerCI12="fitS", pvalue12="fitS",
                                   
                                   varsubj12="fitS", varsubjttt12="fitS")
      
    } else {
  
    res_analyse12 <- data.frame( muest12 <- summary (fit_analyse12)$coefficients[1,1],
                                 P_0_mod12 <- ((exp(muest12))/(1+exp(muest12))),
                                 
                                 betaest12 <- summary (fit_analyse12)$coefficients [2,1],
                                 betaest_se12 <- sqrt(vcov(fit_analyse12)[2,2]),
                                 
                                 P_1_mod12 <- ((exp (betaest12)*exp(muest12))/(1+(exp(betaest12)*exp(muest12)))),
                                 
                                 
                                 upperCI12 <-  betaest12 + 1.96*betaest_se12,
                                 lowerCI12 <-  betaest12  - 1.96*betaest_se12,
                                 
                                 pvalue12 <- summary (fit_analyse12)$coefficients[2,4],
                                 
                                 varsubj12 <- as.numeric (VarCorr(fit_analyse12)[1]) )
    
    }
    # Return matrix containing estimate, SE, CI and p
    return(res_analyse12)
    
  }
  
  temp_analyse12 <- try(fanalyse12(base, Y, TRAIT, SUBJID), TRUE)
  
  
  
  
  # If converges, save output; otherwise substitute 9999s
  
  if(is.list(temp_analyse12)) {
    
    names(temp_analyse12)[1] <- "muest12"
    names(temp_analyse12)[2] <- "P_0_mod12"
    names(temp_analyse12)[3] <- "betaest12"
    names(temp_analyse12)[4] <- "betaest_se12"
    names(temp_analyse12)[5] <- "P_1_mod12"
    names(temp_analyse12)[6] <- "upperCI12"
    names(temp_analyse12)[7] <- "lowerCI12"
    names(temp_analyse12)[8] <- "pvalue12"
    names(temp_analyse12)[9] <- "varsubj12"
    
    analyse12  <- as.data.frame(temp_analyse12)
    
  } else {
    
    analyse12  <- data.frame(muest12='NonC', P_0_mod12='NonC', betaest12='NonC',
                             
                             P_1_mod12='NonC', betaest_se12='NonC',
                             
                             upperCI12='NonC', lowerCI12='NonC', pvalue12='NonC',
                             
                             varsubj12='NonC')
    
  }
  

  
  
  ################################################################################################
  # 1.2 IDENTITY LINK: ESTIMATION OF RISK DIFFERENCE
  ################################################################################################
  
  ################################################################################################
  # 1.2.1 Correlation structure 1 (no model misspecification) 
  ################################################################################################
  
  
  fanalyse2 <- function(base, Y, TRAIT, SUBJID){
    
    # Fit regression
    
    fit_analyse2 <- lmer(Y ~ TRAIT +  (1|SUBJID) + (1 |SUBJID:TRAIT), data= base,
                                      control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000),
                                      check.conv.grad = .makeCC("stop", tol = 2e-3, relTol = NULL) ))
    
    # Pick up estimate
    
    #if singular fit  = T then subsitute fitS
    
    if( (isSingular (fit_analyse2, tol = 1e-4)==T)) {
      
      res_analyse2  <- data.frame(muest2="fitS", P_0_mod2="fitS", betaest2="fitS",
                                   
                                   P_1_mod2="fitS", betaest_se2="fitS",
                                   
                                   upperCI2="fitS", lowerCI2="fitS", pvalue2="fitS",
                                   
                                   varsubj2="fitS", varsubjttt2="fitS")
    } else {
      
    res_analyse2 <- data.frame(  muest2 <- summary (fit_analyse2)$coefficients[1,1],
                                 P_0_mod2 <- ((exp(muest2))/(1+exp(muest2))),
                                 
                                 betaest2 <- summary (fit_analyse2)$coefficients [2,1],
                                 betaest_se2 <- sqrt(vcov(fit_analyse2)[2,2]),
                                 
                                 P_1_mod2 <- ((exp (betaest2)*exp(muest2))/(1+(exp(betaest2)*exp(muest2)))),
                                 
                                 
                                 upperCI2 <-  betaest2 + 1.96*betaest_se2,
                                 lowerCI2 <-  betaest2  - 1.96*betaest_se2,
                                
                                 pvalue2 <- summary (fit_analyse2)$coefficients[2,5],
                                 
                                 varsubj2 <- as.numeric (VarCorr(fit_analyse2)[1]), 
                                 varsubjttt2 <- as.numeric (VarCorr(fit_analyse2)[1]) )
    
    }


    # Return matrix containing estimate, SE, CI and p
    return(res_analyse2)
    
  }
  
  temp_analyse2 <- try(fanalyse2(base, Y, TRAIT, SUBJID), TRUE)
  
  
  
  
  # If converges, save output; otherwise substitute 9999s
  
  if(is.list(temp_analyse2)) {
    
    names(temp_analyse2)[1] <- "muest2"
    names(temp_analyse2)[2] <- "P_0_mod2"
    names(temp_analyse2)[3] <- "betaest2"
    names(temp_analyse2)[4] <- "betaest_se2"
    names(temp_analyse2)[5] <- "P_1_mod2"
    names(temp_analyse2)[6] <- "upperCI2"
    names(temp_analyse2)[7] <- "lowerCI2"
    names(temp_analyse2)[8] <- "pvalue2"
    names(temp_analyse2)[9] <- "varsubj2"
    names(temp_analyse2)[10] <- "varsubjttt2"
    
    analyse2  <- as.data.frame(temp_analyse2)
    
  } else {
    
    analyse2  <- data.frame(muest2='NonC', P_0_mod2='NonC', betaest2='NonC',
                             
                             P_1_mod2='NonC', betaest_se2='NonC',
                            
                             upperCI2='NonC', lowerCI2='NonC', pvalue2='NonC',
                             
                             varsubj2='NonC', varsubjttt2='NonC')
    
  }
  

  
 
  ################################################################################################
  # 1.3.2 Correlation structure 3 (model misspecification) 
  ################################################################################################
  

  fanalyse22 <- function(base, Y, TRAIT, SUBJID){
    
    # Fit regression
    
    fit_analyse22 <- lmer(Y ~ TRAIT +  (1|SUBJID), data= base,
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000),
                                             check.conv.grad = .makeCC("stop", tol = 2e-3, relTol = NULL) ))
    
    # Pick up estimate
    
    #if singular fit  = T then subsitute fitS
    
    if( (isSingular (fit_analyse22, tol = 1e-4)==T)) {
      
      res_analyse22  <- data.frame(muest22="fitS", P_0_mod22="fitS", betaest22="fitS",
                                  
                                  P_1_mod22="fitS", betaest_se22="fitS",
                                  
                                  upperCI22="fitS", lowerCI22="fitS", pvalue22="fitS",
                                  
                                  varsubj22="fitS", varsubjttt22="fitS")
    } else {
      
    res_analyse22 <- data.frame( muest22 <- summary (fit_analyse22)$coefficients[1,1],
                                 P_0_mod22 <- ((exp(muest22))/(1+exp(muest22))),
                                 
                                 betaest22 <- summary (fit_analyse22)$coefficients [2,1],
                                 betaest_se22 <- sqrt(vcov(fit_analyse22)[2,2]),
                                 
                                 P_1_mod22 <- ((exp (betaest22)*exp(muest22))/(1+(exp(betaest22)*exp(muest22)))),
                                 
                                 
                                 upperCI22 <-  betaest22+ 1.96*betaest_se22,
                                 lowerCI22 <-  betaest22  - 1.96*betaest_se22,
                                 
                                 pvalue22 <- summary (fit_analyse22)$coefficients[2,5],
                                 
                                 varsubj22 <- as.numeric (VarCorr(fit_analyse22)[1]))
    
    
    }
    
    # Return matrix containing estimate, SE, CI and p
    return(res_analyse22)
    
  }
  
  temp_analyse22 <- try(fanalyse22(base, Y, TRAIT, SUBJID), TRUE)
  
  
  
  
  # If converges, save output; otherwise substitute 9999s
  
  if(is.list(temp_analyse22)) {
    
    names(temp_analyse22)[1] <- "muest22"
    names(temp_analyse22)[2] <- "P_0_mod22"
    names(temp_analyse22)[3] <- "betaest22"
    names(temp_analyse22)[4] <- "betaest_se22"
    names(temp_analyse22)[5] <- "P_1_mod22"
    names(temp_analyse22)[6] <- "upperCI22"
    names(temp_analyse22)[7] <- "lowerCI22"
    names(temp_analyse22)[8] <- "pvalue22"
    names(temp_analyse22)[9] <- "varsubj22"
    
    analyse22  <- as.data.frame(temp_analyse22)
    
  } else {
    
    analyse22  <- data.frame(muest22='NonC', P_0_mod22='NonC', betaest22='NonC',
                            
                            P_1_mod22='NonC', betaest_se22='NonC',
                            
                            upperCI22='NonC', lowerCI22='NonC', pvalue22='NonC',
                            
                            varsubj22='NonC')
    
  }
  

  

  ################################################################################################
  # 2. ANALYSIS USING MARGINAL APPROACH BASED ON GENERALIZED ESTIMATING EQUATIONS
  ################################################################################################
  
  ################################################################################################
  # 2.1 RELATIVE TREATMENT EFFECT ESTIMATE
  ################################################################################################
  
  #creation of a database to have all the pseudo observations
  
  # creation of database regarding arm of treatment
  baseCONTROL <- base[base$TRAIT=="0",]
  baseTTT <- base[base$TRAIT=="1",]
  baseTTT$Y1 <- baseTTT$Y
  baseCONTROL$Y0 <- baseCONTROL$Y
  
  # I delete usefulness variables 
  baseTTT <- baseTTT[,-2]
  baseTTT <- baseTTT[,-2]
  baseTTT <- baseTTT[,-2]
  baseTTT <- baseTTT[,-2]
  baseCONTROL <- baseCONTROL[,-2]
  baseCONTROL <- baseCONTROL[,-2]
  baseCONTROL <- baseCONTROL[,-2]
  baseCONTROL <- baseCONTROL[,-2]
  
  # I join the base
  baseliang <- full_join(baseTTT,baseCONTROL,by=c("SUBJID"))
  
  #estimation of the difference
  baseliang$diff <- baseliang$Y1- baseliang$Y0
  
  #I delete the variables with 0
  baseliang2 <- baseliang[-which(baseliang$diff==0),]
  
  
  #analysis
  analyse3 <-  geeglm(Y1~1, family=binomial(link = "logit"), corstr="ex",scale.fix = F,
                      data=baseliang2, id=SUBJID)
  
  estimate3 <- as.numeric (summary (analyse3)$coefficients[1])
  estimate_se3 <- as.numeric (summary (analyse3)$coefficients[2])
  
  pvalue3 <- as.numeric (summary (analyse3)$coefficients[4])
  
  #CI
  cc <- coef(summary(analyse3))
  citab <- with(as.data.frame(cc),
                cbind(lwr=Estimate-1.96*Std.err,
                      upr=Estimate+1.96*Std.err))
  rownames(citab) <- rownames(cc)
  citab
  
  lowerCI3 <- citab[1]
  upperCI3 <- citab[2]

  
  OR3 <- exp (estimate3)
  
  
  
  ################################################################################################
  # 2.2 ABSOLUTE TREATMENT EFFECT ESTIMATE
  ################################################################################################
  
  analyse4 <-  geeglm(diff~1, family=gaussian, corstr="ex",scale.fix = F,
                      data=baseliang, id=SUBJID)
  
  estimate4 <- as.numeric (summary (analyse4)$coefficients[1])
  estimate_se4 <- as.numeric (summary (analyse4)$coefficients[2])
  
  pvalue4 <- as.numeric (summary (analyse4)$coefficients[4])
  
  #CI calculation
  cc <- coef(summary(analyse4))
  citab <- with(as.data.frame(cc),
                cbind(lwr=Estimate-1.96*Std.err,
                      upr=Estimate+1.96*Std.err))
  rownames(citab) <- rownames(cc)
  citab
  
  lowerCI4 <- citab[1]
  upperCI4 <- citab[2]
  
  
  
  
  #I add the value I need for after the calcul
  #proportion difference
  diffprop_theo <- P_0 - P_1
  #odds ratio
  OR_theo <- (P_0*(1 - P_1))/ (P_1 * (1 - P_0))
  
  

  
  
  
  return(c(P_0_freq, P_1_freq,
           diffprop_theo, OR_theo,
           analyse1$muest1, analyse1$P_0_mod1, analyse1$betaest1, analyse1$P_1_mod1, analyse1$betaest_se1, analyse1$varsubj1, 
           analyse1$varsubjttt1, analyse1$pvalue1, 
           analyse12$muest12, analyse12$P_0_mod12, analyse12$betaest12, analyse12$P_1_mod12, analyse12$betaest_se12,  
           analyse12$varsubj12,  analyse12$pvalue12, 
           analyse2$muest2, analyse2$P_0_mod2, analyse2$betaest2, analyse2$P_1_mod2, analyse2$betaest_se2,  
           analyse2$varsubj2, analyse2$varsubjttt2, analyse2$pvalue2,
           analyse22$muest22, analyse22$P_0_mod22, analyse22$betaest22, analyse22$P_1_mod22, analyse22$betaest_se22,  
           analyse22$varsubj22,  analyse22$pvalue22, 
           estimate3, estimate_se3, pvalue3, OR3,
           estimate4, estimate_se4, pvalue4 ))
  
  
}


