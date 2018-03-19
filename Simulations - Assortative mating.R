#############
#############
##         ##
## SUMMARY ##
##         ##
#############
#############

#1) User-defined functions
     #1.1) Functions to do summary data Mendelian randomization
          #1.1.1) Two-stage least squares (TSLS) regression
          #1.1.2) Inverse variance weighting (IVW)
          #1.1.3) MR-Egger regression
          #1.1.4) Weighted median
          #1.1.5) MBE
     #1.2) Additional functions
          #1.2.1) NO Measurement Error (NOME) assumption
          #1.2.2) Heterogeneity test 
          #1.2.3) Check whether confidence intervals include zero and/or the true causal effect 

#2) Set up each simulation scenario
     #2.1) Cross-trait assortative mating on X and Y
     #2.2) Cross-trait assortative mating on X_line and Y_line
     #2.3) Cross-trait assortative mating on X_star and Y_star
     #2.4) Single-trait assortative mating on C_line
     #2.5) Single-trait assortative mating on C_star
     #2.6) Single-trait assortative mating on X
     #2.7) Single-trait assortative mating on Y
     #2.8) Cross-trait assortative mating on X and Y along generations
     #2.9) Power to detect an effect

#3) Run the simulations
     #3.1) Loop through each unique combination of parameters
          #3.1.1) Extract the parameters corresponding to the current scenario
          #3.1.2) Generate n_sim datasets, and compute and store results for each
               #3.1.2.1) Simulate P1 mothers
               #3.1.2.2) Simulate fathers, offspring and mothers along generations
                    #3.1.2.2.1) Introduce assortative mating
                    #3.1.2.2.2) Simulate offspring's data
               #3.1.2.3) Compute result
          #3.1.3) Finalise the analysis for the current combination of parameters

#------------------------------------------------------------------------------------------------------------------------------#

#Clean the working environment
rm(list=ls())

#Install/load required packages
if(require('MASS')==F)       {install.packages('MASS'); library('MASS')}
if(require('AER')==F)        {install.packages('AER'); library('AER')}
if(require('doParallel')==F) {install.packages('doParallel'); library('doParallel')}

registerDoParallel(cores=55) #Number of cores to use to perform the simulations

#Inform the directory and the name of the file where you want to save the simulation results
directory <- '/directory'
filename  <- paste(directory, 'filename.txt', sep='/')

#------------------------------------------------------------------------------------------------------------------------------#

###############################
###############################
##                           ##     
## 1) User-defined functions ##
##                           ##
###############################
###############################

###########################################################
#1.1) Functions to do summary data Mendelian randomization#
###########################################################

##################################################
#1.1.1) Two-stage least squares (TSLS) regression#
##################################################

TSLS <- function(model, skip_fathers=F, weights=NULL, X=NULL, Y=NULL, G_X=NULL, X_mothers=NULL, G_X_mothers=NULL, X_fathers=NULL, G_X_fathers=NULL) {
  
  if(is.null(weights)) { weights <- rep(1, length(X)) }

  if(model==1) {
    
    #Generate allele scores
    Z <- apply(t(t(G_X)*weights), 1, sum)

    fit   <- ivreg(Y~X|Z)
    coefs <- 2
  }
  
  if(model==2) {
    
    #Generate allele scores
    Z <- apply(t(t(G_X)*weights), 1, sum)
    Z_mothers <- apply(t(t(G_X_mothers)*weights), 1, sum)
    Z_fathers <- apply(t(t(G_X_fathers)*weights), 1, sum)
    
    fit <- ivreg(Y~X+Z_mothers+Z_fathers|Z+Z_mothers+Z_fathers)
    coefs <- 2:4
  }
  
  if(model==3) {
    
    #Generate allele scores
    Z_mothers <- apply(t(t(G_X_mothers)*weights), 1, sum)
    
    if(skip_fathers==F) {
      Z_fathers <- apply(t(t(G_X_fathers)*weights), 1, sum)    
      fit   <- ivreg(Y~X_mothers+X_fathers|Z_mothers+Z_fathers)
    } else {
      fit   <- ivreg(Y~X_mothers|Z_mothers)
    }
    coefs <- 2:nrow(summary(fit)$coef) 
  }
  
  if(model==4) {

    #Generate allele scores
    Z <- apply(t(t(G_X)*weights), 1, sum)
    Z_mothers <- apply(t(t(G_X_mothers)*weights), 1, sum)
    
    if(skip_fathers==F) {
      Z_fathers <- apply(t(t(G_X_fathers)*weights), 1, sum)    
      fit <- ivreg(Y~X+X_mothers+X_fathers|Z+Z_mothers+Z_fathers)
    } else {
      fit <- ivreg(Y~X+X_mothers|Z+Z_mothers)
    }

    coefs <- 2:nrow(summary(fit)$coef)
  }
  
  #Extract results
  if(length(coefs)==1) {
    TSLSresults <- c(summary(fit)$coef[coefs,1:2], confint(fit)[coefs,], summary(fit)$coef[coefs,3:4])
    names(TSLSresults) <- c('Estimate', 'SE', 'CI_low', 'CI_upp', 'T', 'P-value')
  } else {
    TSLSresults           <- cbind(summary(fit)$coef[coefs,1:2], confint(fit)[coefs,], summary(fit)$coef[coefs,3:4])    
    colnames(TSLSresults) <- c('Estimate', 'SE', 'CI_low', 'CI_upp', 'T', 'P-value')
  }

  if(model==3 & length(coefs)==1) {
    TSLSresults <- rbind(TSLSresults, rep(NA, length(TSLSresults)))
  }

  if(model==4 & length(coefs)==2) {
    TSLSresults <- rbind(TSLSresults, rep(NA, ncol(TSLSresults)))
  }
  
  return(TSLSresults)
}

#########################################
#1.1.2) Inverse variance weighting (IVW)#
#########################################

IVW <- function(BetaXG, BetaYG, seBetaYG) {
  
  IVWfit <- summary(lm(BetaYG ~ -1+BetaXG,weights=1/seBetaYG^2))
  
  # Inference with correct standard errors  
  DF <- length(BetaYG)-1
  IVWBeta <- IVWfit$coef[1,1]
  SE <- IVWfit$coef[1,2]/min(IVWfit$sigma, 1)
  IVW_p <- 2*(1-pt(abs(IVWBeta/SE),DF))
  IVW_CI <- IVWBeta + c(-1,1)*qt(df=DF, 0.975)*SE
  
  # IVWResults = (point estimate, corrected standard error, 
  # 95% Confidence interval, t-statistic, p-value)  
  IVWResults <- c(IVWBeta,SE,IVW_CI,IVWBeta/SE,IVW_p)
  names(IVWResults) <- c('Estimate', 'SE', 'CI_low', 'CI_upp', 'T', 'P-value')  
  
  return(IVWResults)
}

############################
#1.1.3) MR-Egger regression#
############################

MREgger <- function(BetaXG, BetaYG, seBetaYG) {
  
  BYG <- BetaYG*sign(BetaXG) # Pre-processing steps to ensure all  
  BXG <- abs(BetaXG)         # gene--exposure estimates are positive  
  
  ##########
  #MR-Egger#
  ##########
  MREggerFit <- summary(lm(BYG ~ BXG,weights=1/seBetaYG^2))
  
  # Inference with correct standard errors
  MREggerBeta0 <- MREggerFit$coef[1,1]
  MREggerBeta1 <- MREggerFit$coef[2,1]
  SE0 <- MREggerFit$coef[1,2]/min(MREggerFit$sigma, 1)
  SE1 <- MREggerFit$coef[2,2]/min(MREggerFit$sigma, 1)
  DF <- length(BetaYG)-2
  MRBeta0_p <- 2*(1-pt(abs(MREggerBeta0/SE0),DF))
  MRBeta1_p <- 2*(1-pt(abs(MREggerBeta1/SE1),DF))
  MRBeta0_CI <- MREggerBeta0 + c(-1,1)*qt(df=DF, 0.975)*SE0
  MRBeta1_CI <- MREggerBeta1 + c(-1,1)*qt(df=DF, 0.975)*SE1
  
  # MREggerResults = (point estimate, corrected standard error, 
  # 95% Confidence interval, t-statistic, p-value) for
  # intercept (row 1) and slope (row 2). 
  MREggerResults <- matrix(nrow = 2,ncol = 6)
  MREggerResults[1,] <- c(MREggerBeta0,SE0,MRBeta0_CI,MREggerBeta0/SE0,MRBeta0_p)  
  MREggerResults[2,] <- c(MREggerBeta1,SE1,MRBeta1_CI,MREggerBeta1/SE1,MRBeta1_p)
  colnames(MREggerResults) <- c('Estimate', 'SE', 'CI_low', 'CI_upp', 'T', 'P-value')
  rownames(MREggerResults) <- c('Intercept', 'Beta')  
  
  return(MREggerResults)
}

########################
#1.1.4) Weighted median#
########################

WeightedMedian <- function(BetaXG, BetaYG, seBetaXG, seBetaYG, alpha=0.05) {
  
  #Define a function to obtain the WM estimate
  weighted.median <- function(betaIV.in, weights.in) {
    betaIV.order <- betaIV.in[order(betaIV.in)]
    weights.order <- weights.in[order(betaIV.in)]
    weights.sum <- cumsum(weights.order)-0.5*weights.order
    weights.sum <- weights.sum/sum(weights.order)
    below <- max(which(weights.sum<0.5))
    weighted.est <- betaIV.order[below] + (betaIV.order[below+1]-betaIV.order[below])*
      (0.5-weights.sum[below])/(weights.sum[below+1]-weights.sum[below])
    return(weighted.est)
  }
  
  #Define a function to obtain the WM estimate SE through boostraping
  weighted.median.boot <- function(betaXG.in, betaYG.in, sebetaXG.in, sebetaYG.in, weights.in) {
    med <- NULL
    for(i in 1:1000) {
      betaXG.boot <- rnorm(length(betaXG.in), mean=betaXG.in, sd=sebetaXG.in)
      betaYG.boot <- rnorm(length(betaYG.in), mean=betaYG.in, sd=sebetaYG.in)
      betaIV.boot <- betaYG.boot/betaXG.boot
      med[i] <- weighted.median(betaIV.boot, weights.in)
    }
    return(sd(med))
  }
  
  betaIV <- BetaYG/BetaXG #ratio estimates
  weights <- (seBetaYG/BetaXG)^-2 #inverse-variance weights
  
  #Obtain WM estimate and SE
  betaWM <- weighted.median(betaIV, weights)
  
  sebetaWM <- weighted.median.boot(BetaXG, BetaYG,
                                   seBetaXG, seBetaYG, weights)
  
  #Calculate WM 95% CIs using the normal approximation
  betaWM_CI <- betaWM + c(-1,1)*qnorm(1-alpha/2)*sebetaWM
  
  #Calculate P-value
  pWM <- pt(abs(betaWM/sebetaWM), df=length(BetaXG)-1, lower.tail=F)*2
  
  WMResults <- c(betaWM, sebetaWM, betaWM_CI, betaWM/sebetaWM, pWM)  
  names(WMResults) <- c('Estimate', 'SE', 'CI_low', 'CI_upp', 'T', 'P-value')
  
  return(WMResults)
}

#######################
#1.1.5) Mode estimator#
#######################

MBE <- function(BetaXG, BetaYG, seBetaXG, seBetaYG, phi=1, n_boot=1e3, alpha=0.05) {
  
  #--------------------------------------#
  #Function to compute the point estimate#
  #--------------------------------------#
  #BetaIV.in: ratio estimates
  #seBetaIV.in: standard errors of ratio estimates
  beta <- function(BetaIV.in, seBetaIV.in) {
    
    #Bandwidth rule - modified Silverman's rule proposed by Bickel (2002)
    s <- 0.9*(min(sd(BetaIV.in), mad(BetaIV.in)))/length(BetaIV.in)^(1/5)
    
    #Standardised weights
    weights <- seBetaIV.in^-2/sum(seBetaIV.in^-2)
    
    beta <- NULL
    
    for(cur_phi in phi) {
      #Define the actual bandwidth
      h <- s*cur_phi
      #Compute the smoothed empirical density function
      densityIV <- density(BetaIV.in, weights=weights, bw=h)
      #Extract the point with the highest density as the point estimate 
      beta[length(beta)+1] <- densityIV$x[densityIV$y==max(densityIV$y)]
    }
    return(beta)
  }
  
  #------------------------------------------#
  #Function to estimate SEs through bootstrap#
  #------------------------------------------#
  #BetaIV.in: ratio estimates
  #seBetaIV.in: standard errors of ratio estimates
  #beta_Mode.in: point causal effect estimates
  boot <- function(BetaIV.in, seBetaIV.in, beta_Mode.in) {
    
    #Set up a matrix to store the results from each bootstrap iteration
    beta.boot <- numeric(n_boot)
    
    for(i in 1:n_boot) {
      #Re-sample each ratio estimate using SEs derived not assuming NOME
      BetaIV.boot      <- rnorm(length(BetaIV.in), mean=BetaIV.in, sd=seBetaIV.in)
      
      #Weighted mode, not assuming NOME
      beta.boot[i] <- beta(BetaIV.in=BetaIV.boot, seBetaIV.in=seBetaIV.in)
    }
    return(beta.boot)
  }
  
  #Ratio estimates
  BetaIV   <- BetaYG/BetaXG    
  
  #SEs of ratio estimates
  seBetaIV <- seBetaYG/abs(BetaXG)
  
  #Point causal effect estimate using the weighted mode
  beta_Mode <- beta(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV)
  
  #Compute SEs, confidence intervals and P-value
  beta_Mode.boot <- boot(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV, beta_Mode.in=beta_Mode)
  se_Mode <- mad(beta_Mode.boot)
  
  CIlow_Mode <- beta_Mode-qnorm(1-alpha/2)*se_Mode
  CIupp_Mode <- beta_Mode+qnorm(1-alpha/2)*se_Mode
  
  P_Mode <- pt(abs(beta_Mode/se_Mode), df=length(BetaXG)-1, lower.tail=F)*2
  
  #Return a data frame containing the results
  Results <- data.frame(beta_Mode, se_Mode, CIlow_Mode, CIupp_Mode, P_Mode)  
  colnames(Results) <- c('Estimate', 'SE', 'CI_low', 'CI_upp', 'P')
  
  return(Results)
}

############################
#1.2) Additional functions #
############################

################################################
#1.2.1) NO Measurement Error (NOME) assumption #
################################################

NOME <- function(BetaXG, seBetaXG, seBetaYG, FXG=NULL) {
  
  if(is.null(FXG)) {
    FXG <- (BetaXG/seBetaXG)^2
  }
  
  #IVW
  weights <- (seBetaYG/abs(BetaXG))^-2 
  FXG     <- (weighted.mean(FXG, weights)-1)/weighted.mean(FXG, weights)
  
  #MR-Egger
  Isq <- function(y,s) {
    k <- length(y)
    w <- 1/s^2; sum.w <- sum(w)
    mu.hat <- sum(y*w)/sum.w
    Q <- sum(w*(y-mu.hat)^2)
    Isq <- (Q-(k-1))/Q
    Isq <- max(0, Isq)
    return(Isq)
  }
  Isq_res <- Isq(abs(BetaXG)/seBetaYG, seBetaXG/seBetaYG)
  
  return(c(FXG=FXG, Isq_res=Isq_res)) 
}

############################
#1.2.2) Heterogeneity test #
############################

HeterogeneityTest <- function(BetaXG, BetaYG, seBetaXG, seBetaYG) {
  
  #Cocharn's Q
  BetaIV      <- BetaYG/BetaXG                   #Ratio estimates
  weights1    <- (seBetaYG/abs(BetaXG))^-2       #1st order weights
  IVW_beta_w1 <- weighted.mean(BetaIV, weights1) #IVW point estimate using 1st order weights
  
  weights2 <- 1/(seBetaYG^2/BetaXG^2 + (IVW_beta_w1^2)*seBetaXG^2/BetaXG^2) #Modified 2nd order weights
  BIVw2    <- BetaIV*sqrt(weights2)
  sW2      <- sqrt(weights2)
  IVWfit   <- summary(lm(BIVw2 ~ -1+sW2)) #IVW model with modified 2nd order weights
  
  DF         <- length(BetaIV)-1
  phi_IVWfit <- IVWfit$sigma^2
  Q_stat     <- DF*phi_IVWfit                  #Q statistic
  I2         <- max(0, (Q_stat-DF)/Q_stat)*100 #I2
  Q_p        <- 1-pchisq(Q_stat, DF)           #P-value
  
  return(c(Q_stat, I2, Q_p))
}

#######################################################################################
#1.2.3) Check whether confidence intervals include zero and/or the true causal effect #
#######################################################################################

coverage_power <- function(CI_low, CI_upp, beta) {
  
  coverage <- CI_low<=beta & beta<=CI_upp
  power    <- CI_low>=0 | CI_upp<=0
  return(c(coverage, power))
}

########################################
########################################
##                                    ##
## 2) Set up each simulation scenario ##
##                                    ##
########################################
########################################

#Set a dataframe for varying parameters
parameters <- data.frame(matrix(nrow=2e3, ncol=22)) #Set up an empty dataframe to store the parameters

colnames(parameters) <- c('scenario',       #Simulations scenario
                          'omega',          #Type of assortative mating
                          'zeta',           #Number of generations
                          'n_X',            #Sample size of the dataset used to estimate instrument-exposure associations
                          'n_Y',            #Sample size of the dataset used to estimate instrument-outcome associations
                          'L_X',            #Number of independent genetic variants that cause the exposure
                          'L_Y',            #Number of independent genetic variants that cause the outcome,
                          'P',              #Strength of cross-trait assortative mating
                          'H',              #Number of genetic instruments
                          'delta_X',        #Collective additive effect of the L_X genetic variants associated with X on X
                          'delta_Y',        #Collective additive effect of the L_Y genetic variants associated with Y on Y
                          'delta_X_line',   #Collective additive effect of the L_X genetic variants associated with X on X_line
                          'delta_Y_line',   #Collective additive effect of the L_Y genetic variants associated with Y on Y_line
                          'delta_X_C_line', #Collective additive effect of the L_X genetic variants associated with X on C_line
                          'delta_Y_C_line', #Collective additive effect of the L_Y genetic variants associated with Y on C_line
                          'kappa_X_star',   #Effect of the exposure X on X_star
                          'kappa_C_star',   #Effect of the exposure X on C_star
                          'lambda_Y_star',  #Effect of the outcome Y on Y_star
                          'lambda_C_star',  #Effect of the outcome Y on C_star
                          'theta_X',        #Effect of the confounder U on the exposure X
                          'theta_Y',        #Effect of the confounder U on the outcome Y
                          'beta'            #Effect of the exposure X on the outcome Y
)

count <- 0 #A counter to facilitate adding each combination of parameters to the parameters dataframe

#Function to facilitate filling in the parameters of each simulation scenario
set_simulation_scenario <- function(parameters, count,
                                    scenario, omega, zeta, n_X, n_Y, L_X, L_Y, P, H,
                                    delta_X, delta_Y, delta_X_line, delta_Y_line, delta_X_C_line, delta_Y_C_line,
                                    kappa_X_star, kappa_C_star, lambda_Y_star, lambda_C_star,
                                    theta_X, theta_Y, beta) {
  
  for(cur_scenario in scenario) {
    for(cur_omega in omega) {
      for(cur_zeta in zeta) {
        for(cur_n_X in n_X) {
          for(cur_n_Y in n_Y) {
            for(cur_L_X in L_X) {
              for(cur_L_Y in L_Y) {
                for(cur_P in P) {
                  for(cur_H in H) {
                    for(cur_delta_X in delta_X) {
                      for(cur_delta_Y in delta_Y) {
                        for(cur_delta_X_line in delta_X_line) {
                          for(cur_delta_Y_line in delta_Y_line) {
                            for(cur_delta_X_C_line in delta_X_C_line) {
                              for(cur_delta_Y_C_line in delta_Y_C_line) {
                                for(cur_kappa_X_star in kappa_X_star) {
                                  for(cur_kappa_C_star in kappa_C_star) {
                                    for(cur_lambda_Y_star in lambda_Y_star) {
                                      for(cur_lambda_C_star in lambda_C_star) {
                                        for(cur_theta_X in theta_X) {
                                          for(cur_theta_Y in theta_Y) {
                                            for(cur_beta in beta) {
                                              count <- count+1
                                              
                                              parameters$scenario[count]       <- cur_scenario
                                              parameters$omega[count]          <- cur_omega
                                              parameters$zeta[count]           <- cur_zeta
                                              parameters$n_X[count]            <- cur_n_X
                                              parameters$n_Y[count]            <- cur_n_Y
                                              parameters$L_X[count]            <- cur_L_X
                                              parameters$L_Y[count]            <- cur_L_Y
                                              parameters$P[count]              <- cur_P
                                              parameters$H[count]              <- cur_H
                                              parameters$delta_X[count]        <- cur_delta_X
                                              parameters$delta_Y[count]        <- cur_delta_Y
                                              parameters$delta_X_line[count]   <- cur_delta_X_line
                                              parameters$delta_Y_line[count]   <- cur_delta_Y_line
                                              parameters$delta_X_C_line[count] <- cur_delta_X_C_line
                                              parameters$delta_Y_C_line[count] <- cur_delta_Y_C_line
                                              parameters$kappa_X_star[count]   <- cur_kappa_X_star
                                              parameters$kappa_C_star[count]   <- cur_kappa_C_star
                                              parameters$lambda_Y_star[count]  <- cur_lambda_Y_star
                                              parameters$lambda_C_star[count]  <- cur_lambda_C_star
                                              parameters$theta_X[count]        <- cur_theta_X
                                              parameters$theta_Y[count]        <- cur_theta_Y
                                              parameters$beta[count]           <- cur_beta
                                            }
                                          } 
                                        }
                                      }
                                    }
                                  }
                                }   
                              } 
                            } 
                          } 
                        } 
                      }
                    }
                  }
                }
              } 
            }
          }
        }
      }
    }
  }
  return(list(parameters,count))
}

################################################
#2.1) Cross-trait assortative mating on X and Y#
################################################

simulation_1 <- set_simulation_scenario(parameters=parameters, count=count,
                                        scenario       = 1,
                                        omega          = 1,
                                        zeta           = 1,
                                        n_X            = 2e4,
                                        n_Y            = 2e4,
                                        L_X            = 50,
                                        L_Y            = 50,
                                        P              = seq(0,1,0.2),
                                        H              = 50,
                                        delta_X        = c(sqrt(0.1), sqrt(0.5)),
                                        delta_Y        = c(sqrt(0.1), sqrt(0.5)),
                                        delta_X_line   = 0,
                                        delta_Y_line   = 0,
                                        delta_X_C_line = 0,
                                        delta_Y_C_line = 0,
                                        kappa_X_star   = 0,
                                        kappa_C_star   = 0,
                                        lambda_Y_star  = 0,
                                        lambda_C_star  = 0,
                                        theta_X        = sqrt(0.4),
                                        theta_Y        = sqrt(0.4),
                                        beta           = 0
)
parameters <- simulation_1[[1]]
count      <- simulation_1[[2]]

##########################################################
#2.2) Cross-trait assortative mating on X_line and Y_line#
##########################################################

simulation_2 <- set_simulation_scenario(parameters=parameters, count=count,
                                        scenario       = 2,
                                        omega          = 2,
                                        zeta           = 1,
                                        n_X            = 2e4,
                                        n_Y            = 2e4,
                                        L_X            = 50,
                                        L_Y            = 50,
                                        P              = 1,
                                        H              = 50,
                                        delta_X        = sqrt(0.5),
                                        delta_Y        = sqrt(0.5),
                                        delta_X_line   = sqrt(c(0, 0.5)),
                                        delta_Y_line   = sqrt(c(0, 0.5)),
                                        delta_X_C_line = 0,
                                        delta_Y_C_line = 0,
                                        kappa_X_star   = 0,
                                        kappa_C_star   = 0,
                                        lambda_Y_star  = 0,
                                        lambda_C_star  = 0,
                                        theta_X        = sqrt(0.4),
                                        theta_Y        = sqrt(0.4),
                                        beta           = 0
)
parameters <- simulation_2[[1]]
count      <- simulation_2[[2]]

##########################################################
#2.3) Cross-trait assortative mating on X_star and Y_star#
##########################################################

simulation_3 <- set_simulation_scenario(parameters=parameters, count=count,
                                        scenario       = 3,
                                        omega          = 3,
                                        zeta           = 1,
                                        n_X            = 2e4,
                                        n_Y            = 2e4,
                                        L_X            = 50,
                                        L_Y            = 50,
                                        P              = 1,
                                        H              = 50,
                                        delta_X        = sqrt(0.5),
                                        delta_Y        = sqrt(0.5),
                                        delta_X_line   = 0,
                                        delta_Y_line   = 0,
                                        delta_X_C_line = 0,
                                        delta_Y_C_line = 0,
                                        kappa_X_star   = sqrt(c(0, 0.5)),
                                        kappa_C_star   = 0,
                                        lambda_Y_star  = sqrt(c(0, 0.5)),
                                        lambda_C_star  = 0,
                                        theta_X        = sqrt(0.4),
                                        theta_Y        = sqrt(0.4),
                                        beta           = 0
)
parameters <- simulation_3[[1]]
count      <- simulation_3[[2]]

################################################
#2.4) Single-trait assortative mating on C_line#
################################################

simulation_4 <- set_simulation_scenario(parameters=parameters, count=count,
                                        scenario       = 4,
                                        omega          = 4,
                                        zeta           = 1,
                                        n_X            = 2e4,
                                        n_Y            = 2e4,
                                        L_X            = 50,
                                        L_Y            = 50,
                                        P              = 1,
                                        H              = 50,
                                        delta_X        = sqrt(0.5),
                                        delta_Y        = sqrt(0.5),
                                        delta_X_line   = 0,
                                        delta_Y_line   = 0,
                                        delta_X_C_line = sqrt(c(0, 0.4)),
                                        delta_Y_C_line = sqrt(c(0, 0.4)),
                                        kappa_X_star   = 0,
                                        kappa_C_star   = 0,
                                        lambda_Y_star  = 0,
                                        lambda_C_star  = 0,
                                        theta_X        = sqrt(0.4),
                                        theta_Y        = sqrt(0.4),
                                        beta           = 0
)
parameters <- simulation_4[[1]]
count      <- simulation_4[[2]]

################################################
#2.5) Single-trait assortative mating on C_star#
################################################

simulation_5 <- set_simulation_scenario(parameters=parameters, count=count,
                                        scenario       = 5,
                                        omega          = 5,
                                        zeta           = 1,
                                        n_X            = 2e4,
                                        n_Y            = 2e4,
                                        L_X            = 50,
                                        L_Y            = 50,
                                        P              = 1,
                                        H              = 50,
                                        delta_X        = sqrt(0.5),
                                        delta_Y        = sqrt(0.5),
                                        delta_X_line   = 0,
                                        delta_Y_line   = 0,
                                        delta_X_C_line = 0,
                                        delta_Y_C_line = 0,
                                        kappa_X_star   = 0,
                                        kappa_C_star   = sqrt(c(0, 0.3)),
                                        lambda_Y_star  = 0,
                                        lambda_C_star  = sqrt(c(0, 0.3)),
                                        theta_X        = sqrt(0.4),
                                        theta_Y        = sqrt(0.4),
                                        beta           = 0
)
parameters <- simulation_5[[1]]
count      <- simulation_5[[2]]

###########################################
#2.6) Single-trait assortative mating on X#
###########################################

simulation_6 <- set_simulation_scenario(parameters=parameters, count=count,
                                        scenario       = 6,
                                        omega          = 6,
                                        zeta           = 1,
                                        n_X            = 2e4,
                                        n_Y            = 2e4,
                                        L_X            = 50,
                                        L_Y            = 50,
                                        P              = 1,
                                        H              = 50,
                                        delta_X        = sqrt(0.5),
                                        delta_Y        = sqrt(0.5),
                                        delta_X_line   = 0,
                                        delta_Y_line   = 0,
                                        delta_X_C_line = 0,
                                        delta_Y_C_line = 0,
                                        kappa_X_star   = 0,
                                        kappa_C_star   = 0,
                                        lambda_Y_star  = 0,
                                        lambda_C_star  = 0,
                                        theta_X        = sqrt(0.4),
                                        theta_Y        = sqrt(0.4),
                                        beta           = 0
)
parameters <- simulation_6[[1]]
count      <- simulation_6[[2]]

###########################################
#2.7) Single-trait assortative mating on Y#
###########################################

simulation_7 <- set_simulation_scenario(parameters=parameters, count=count,
                                        scenario       = 7,
                                        omega          = 7,
                                        zeta           = 1,
                                        n_X            = 2e4,
                                        n_Y            = 2e4,
                                        L_X            = 50,
                                        L_Y            = 50,
                                        P              = 1,
                                        H              = 50,
                                        delta_X        = sqrt(0.5),
                                        delta_Y        = sqrt(0.5),
                                        delta_X_line   = 0,
                                        delta_Y_line   = 0,
                                        delta_X_C_line = 0,
                                        delta_Y_C_line = 0,
                                        kappa_X_star   = 0,
                                        kappa_C_star   = 0,
                                        lambda_Y_star  = 0,
                                        lambda_C_star  = 0,
                                        theta_X        = sqrt(0.4),
                                        theta_Y        = sqrt(0.4),
                                        beta           = 0
)
parameters <- simulation_7[[1]]
count      <- simulation_7[[2]]

##################################################################
#2.8) Cross-trait assortative mating on X and Y along generations#
##################################################################

simulation_8 <- set_simulation_scenario(parameters=parameters, count=count,
                                        scenario       = 8,
                                        omega          = 1,
                                        zeta           = seq(1,10,2),
                                        n_X            = 2e4,
                                        n_Y            = 2e4,
                                        L_X            = 50,
                                        L_Y            = 50,
                                        P              = c(0.5,1),
                                        H              = 50,
                                        delta_X        = c(sqrt(0.1), sqrt(0.5)),
                                        delta_Y        = c(sqrt(0.1), sqrt(0.5)),
                                        delta_X_line   = 0,
                                        delta_Y_line   = 0,
                                        delta_X_C_line = 0,
                                        delta_Y_C_line = 0,
                                        kappa_X_star   = 0,
                                        kappa_C_star   = 0,
                                        lambda_Y_star  = 0,
                                        lambda_C_star  = 0,
                                        theta_X        = sqrt(0.4),
                                        theta_Y        = sqrt(0.4),
                                        beta           = 0
)
parameters <- simulation_8[[1]]
count      <- simulation_8[[2]]

################################
#2.9) Power to detect an effect#
################################

simulation_9 <- set_simulation_scenario(parameters=parameters, count=count,
                                        scenario       = 9,
                                        omega          = 1,
                                        zeta           = 1,
                                        n_X            = 2e4,
                                        n_Y            = 2e4,
                                        L_X            = 50,
                                        L_Y            = 50,
                                        P              = 0,
                                        H              = 50,
                                        delta_X        = c(sqrt(0.1), sqrt(0.2)),
                                        delta_Y        = c(sqrt(0.1), sqrt(0.2)),
                                        delta_X_line   = 0,
                                        delta_Y_line   = 0,
                                        delta_X_C_line = 0,
                                        delta_Y_C_line = 0,
                                        kappa_X_star   = 0,
                                        kappa_C_star   = 0,
                                        lambda_Y_star  = 0,
                                        lambda_C_star  = 0,
                                        theta_X        = sqrt(0.4),
                                        theta_Y        = sqrt(0.4),
                                        beta           = 0.05
)
parameters <- simulation_9[[1]]
count      <- simulation_9[[2]]

parameters <- parameters[!is.na(parameters$n_X),] #Remove empty rows
rm(list=c('count', ls(pattern='simulation_')))

############################
############################
##                        ##
## 3) Run the simulations ##
##                        ##
############################
############################

n_sim <- 5e3 #Number of simulated datasets per combination of parameters

#########################################################
#3.1) Loop through each unique combination of parameters#
#########################################################

names_simulation_results <- c(colnames(parameters),
                              'r2X', 'r2Y',
                              'varX', 'varX_line', 'varX_star', 'varY', 'varY_line', 'varY_star',
                              'varC_line', 'varC_star', 'varZ_X', 'varZ_Y',
                              'corPhenoParents_XX', 'corPhenoParents_YY', 'corPhenoParents_XY', 'corPhenoParents_YX', 'corPheno_XY', 'corGeno_XY', 'corGeno_X', 'corGeno_Y',
                              paste(rep(c(paste('TSLS', 1:9, sep=''), 'IVW', 'MREgger', 'WeightedMedian', 'ModeEstimator'), each=4), c('beta', 'se', 'coverage', 'power'), sep='_'),
                              'FXG', 'IGX2', 'MREggerInt_beta', 'MREggerInt_power', 'cochran_q', 'cochran_power')

if(filename%in%list.files(directory, full.names=T)) {
  results      <- read.table(filename, header=T, sep='\t')
  nrow_results <- nrow(results)
} else {
  results      <- NULL
  nrow_results <- 0
}

for(a in (nrow_results+1):nrow(parameters)) {
  
  print(a) 
  
  #####################################################################
  #3.1.1) Extract the parameters corresponding to the current scenario#
  #####################################################################
  
  scenario       <- parameters$scenario[a]
  omega          <- parameters$omega[a]
  zeta           <- parameters$zeta[a]
  n_X            <- parameters$n_X[a]
  n_Y            <- parameters$n_Y[a]
  n              <- n_X+n_Y
  L_X            <- parameters$L_X[a]
  L_Y            <- parameters$L_Y[a]
  sum_L          <- L_X+L_Y
  P              <- parameters$P[a]
  H              <- parameters$H[a]
  delta_X        <- parameters$delta_X[a]
  delta_Y        <- parameters$delta_Y[a]
  delta_X_line   <- parameters$delta_X_line[a]
  delta_Y_line   <- parameters$delta_Y_line[a]
  delta_X_C_line <- parameters$delta_X_C_line[a]
  delta_Y_C_line <- parameters$delta_Y_C_line[a]
  kappa_X_star   <- parameters$kappa_X_star[a]
  kappa_C_star   <- parameters$kappa_C_star[a]
  lambda_Y_star  <- parameters$lambda_Y_star[a]
  lambda_C_star  <- parameters$lambda_C_star[a]
  theta_X        <- parameters$theta_X[a]
  theta_Y        <- parameters$theta_Y[a]
  beta           <- parameters$beta[a]
  
  ########################################################################
  #3.1.2) Generate n_sim datasets, and compute and store results for each#
  ########################################################################
  
  simulation_results <- foreach(cur_sim=1:n_sim, .combine='rbind', .inorder=F, .packages=c('MASS', 'AER')) %dopar% {
    
    SimulateData <- function(G_X, G_Y) {
      
      Z_X     <- apply(G_X*gamma_X, 1, sum) #Weighted additive allele score combining all exposure-associated genetic variants
      Z_Y     <- apply(G_Y*gamma_Y, 1, sum) #Weighted additive allele score combining all outcome-associated genetic variants

      U <- rnorm(n) #Confounder of the exposure-outcome association
      
      X      <- theta_X*U + delta_X*Z_X + rnorm(n, 0, sqrt(1-theta_X^2-delta_X^2)) #Exposure phenotype
      X_line <- delta_X_line*Z_X + rnorm(n, 0, sqrt(1-delta_X_line^2))             #Phenotype also influenced directly by Z_X
      X_star <- kappa_X_star*X + rnorm(n, 0, sqrt(1-kappa_X_star^2))               #Phenotype caused by X

      Y      <- beta*X + theta_Y*U + delta_Y*Z_Y + rnorm(n, 0, sqrt(1-beta^2-theta_Y^2-delta_Y^2-2*theta_X*beta*theta_Y)) #Outcome phenotype
      Y_line <- delta_Y_line*Z_Y + rnorm(n, 0, sqrt(1-delta_Y_line^2))                                            #Phenotype also influenced directly by Z_Y
      Y_star <- lambda_Y_star*Y + rnorm(n, 0 , sqrt(1-lambda_Y_star^2))                                           #Phenotype caused by Y

      C_line <- delta_X_C_line*Z_X + delta_Y_C_line*Z_Y + rnorm(n, 0, sqrt(1-delta_X_C_line^2-delta_Y_C_line^2)) #Phenotype influenced directly by both Z_X and Z_Y
      C_star <- kappa_C_star*X + lambda_C_star*Y + rnorm(n, 0, sqrt(1-(kappa_C_star^2+lambda_C_star^2+2*((theta_X)*(theta_Y)+2*beta*theta_X*theta_Y)*kappa_C_star*kappa_C_star)))        #Phenotype caused by both X and Y

      return(data.frame(X, X_line, X_star, Y, Y_line, Y_star, C_line, C_star, U, Z_X, Z_Y, G_X, G_Y))
    }
    
    ##############################
    #3.1.2.1) Simulate P1 mothers#
    ##############################
    
    #Genetic variants associated with exposure
    G_X   <- matrix(nrow=n, ncol=L_X)
    phi_X <- runif(L_X, 0.1, 0.9) #Effect allele frequencies
    for(k in 1:L_X) {
      G_X[,k] <- rbinom(n, 2, phi_X[k])
    }
    G_X_names     <- paste('G_X_', 1:L_X, sep='')
    colnames(G_X) <- G_X_names 
    
    #Calculate weights of the L_X X-associated genetic variants, so that a weighted allele score naturally has an expected variance of 1
    sigma_X <- 2*phi_X*(1-phi_X)    #Expected variance of each genetic variant
    a_X     <- 1/L_X/10             #Minimum individual r2
    b_X     <- 2/L_X-a_X            #Maximum individual r2
    r2_k    <- runif(L_X, a_X, b_X) #r2 of each X-associated genetic variant
    gamma_X <- sqrt(r2_k/sigma_X)   #Direct effect on X of each X-associated genetic variant
    gamma_X <- matrix(rep(gamma_X, nrow=n), byrow=T, ncol=L_X, nrow=n)
    
    #Genetic variants associated with outcome
    G_Y   <- matrix(nrow=n, ncol=L_Y)
    phi_Y <- runif(L_Y, 0.1, 0.9) #Effect allele frequencies
    for(k in 1:L_Y) {
      G_Y[,k] <- rbinom(n, 2, phi_Y[k])
    }
    G_Y_names     <- paste('G_Y_', 1:L_Y, sep='')
    colnames(G_Y) <- G_Y_names 
    
    #Calculate weights of the L_Y Y-associated genetic variants, so that a weighted allele score naturally has an expected variance of 1
    sigma_Y <- 2*phi_Y*(1-phi_Y)    #Expected variance of each genetic variant
    a_Y     <- 1/L_Y/10             #Minimum individual r2
    b_Y     <- 2/L_Y-a_Y            #Maximum individual r2
    r2_k    <- runif(L_Y, a_Y, b_Y) #r2 of each X-associated genetic variant
    gamma_Y <- sqrt(r2_k/sigma_Y)   #Direct effect on Y of each Y-associated genetic variant
    gamma_Y <- matrix(rep(gamma_Y, nrow=n), byrow=T, ncol=L_Y, nrow=n)
    
    phi_X_Y <- c(phi_X, phi_Y)
  
    #Combine data in a single data frame
    data_mothers <- SimulateData(G_X, G_Y); rm(list=c('G_X', 'G_Y'))
    
    ####################################################################
    #3.1.2.2) Simulate fathers, offspring and mothers along generations#
    ####################################################################
    
    cur_zeta <- 1
    while(cur_zeta<=zeta) {
      
      #Compute correlation coefficients between all genetic variants
      G_matrix <- data_mothers[,substr(colnames(data_mothers), 1, 3)%in%c('G_X', 'G_Y')]  
      G_cor_matrix <- matrix(nrow=sum_L, ncol=sum_L)
      for(b in 1:sum_L) {
        for(c in b:sum_L) {
          G_cor <- cor(G_matrix[,b], G_matrix[,c])
          G_cor_matrix[b,c] <- G_cor
          G_cor_matrix[c,b] <- G_cor
        } 
      }; rm(G_matrix)
      
      #Simulate paternal genetic variants with a similar genetic correlation pattern than the mothers
      #Using a multivariate normal rather than a multivariate binomial distribution to reduce computational time
      temp_G  <- mvrnorm(n=n, mu=rep(0, sum_L), Sigma=G_cor_matrix); rm(G_cor_matrix)
      G       <- matrix(nrow=n, ncol=sum_L)  
      for(d in 1:sum_L) {
        cur_phi <- phi_X_Y[d]
        cur_G   <- temp_G[,d]
        cur_G   <- cut(cur_G, include.lowest=T, quantile(cur_G, probs=c(0,
                                                                        (1-cur_phi)^2,
                                                                        (1-cur_phi)^2+2*cur_phi*(1-cur_phi),
                                                                        1)), labels=0:2)
        cur_G <- as.numeric(as.character(cur_G))
        G[,d] <- cur_G 
      }; rm(temp_G)
      
      #Genetic variants associated with exposure
      G_X <- G[,1:L_X]; 
      colnames(G_X) <- G_X_names
      
      #Genetic variants associated with outcome
      G_Y <- G[,(L_X+1):sum_L]; rm(G)
      colnames(G_Y) <- G_Y_names
      
      data_fathers <- SimulateData(G_X, G_Y); rm(list=c('G_X', 'G_Y'))
      
      #########################################
      #3.1.2.2.1) Introduce assortative mating#
      #########################################
      
      #Sort maternal and paternal datasets to induce a perfect positive correlation between mother's and father's phenotypes
      if(omega==1) {
        phenotypes <- c('X', 'Y')

      } else if(omega==2) {
        phenotypes <- c('X_line', 'Y_line')

      } else if(omega==3) {
        phenotypes <- c('X_star', 'Y_star')
        
      } else if(omega==4) {
        phenotypes <- c('C_line', 'C_line')

      } else if(omega==5) {
        phenotypes <- c('C_star', 'C_star')

      } else if(omega==6) {
        phenotypes <- c('X', 'X')
        
      } else if(omega==7) {
        phenotypes <- c('Y', 'Y')
      }

      ind_1 <- sample(1:n, n/2, replace=F) #Subset 1
      ind_2 <- (1:n)[!(1:n)%in%ind_1]      #Subset 2
      
      #Generate imperfect (unless P=1) versions of the phenotypes under assortment
      data_mothers$phenotype1 <- P*data_mothers[,phenotypes[1]]
      data_mothers$phenotype2 <- P*data_mothers[,phenotypes[2]]

      data_fathers$phenotype1 <- P*data_fathers[,phenotypes[1]]
      data_fathers$phenotype2 <- P*data_fathers[,phenotypes[2]]

      if(P<1) {
        data_mothers$phenotype1 <- data_mothers$phenotype1 + rnorm(n, 0, sqrt(1-var(data_mothers$phenotype1)))
        data_mothers$phenotype2 <- data_mothers$phenotype2 + rnorm(n, 0, sqrt(1-var(data_mothers$phenotype2)))

        data_fathers$phenotype1 <- data_fathers$phenotype1 + rnorm(n, 0, sqrt(1-var(data_fathers$phenotype1)))
        data_fathers$phenotype2 <- data_fathers$phenotype2 + rnorm(n, 0, sqrt(1-var(data_fathers$phenotype2)))
      }

      #Subset 1: women with high phenotype 1 tend to partner with men with high phenotype 2
      data_mothers1 <- data_mothers[ind_1,]
      data_fathers1 <- data_fathers[ind_1,]
      data_mothers1 <- data_mothers1[order(data_mothers1$phenotype1),]
      data_fathers1 <- data_fathers1[order(data_fathers1$phenotype2),]

      #Subset 2: women with high phenotype 2 tend to partner with men with high phenotype 1
      data_mothers2 <- data_mothers[ind_2,]
      data_fathers2 <- data_fathers[ind_2,]
      data_mothers2 <- data_mothers2[order(data_mothers2$phenotype2),]
      data_fathers2 <- data_fathers2[order(data_fathers2$phenotype1),]

      #Combine subsets
      data_mothers <- rbind(data_mothers1, data_mothers2)
      data_fathers <- rbind(data_fathers1, data_fathers2)

      ind          <- sample(1:n, replace=F)
      data_mothers <- data_mothers[ind,]  
      data_fathers <- data_fathers[ind,]  
      
      ######################################
      #3.1.2.2.2) Simulate offspring's data#
      ######################################
      
      #Simulate offspring's exposure-associated genetic variants based on their parents
      G_mothers <- data_mothers[,substr(colnames(data_mothers), 1, 3)%in%c('G_X', 'G_Y')]
      G_fathers <- data_fathers[,substr(colnames(data_fathers), 1, 3)%in%c('G_X', 'G_Y')]
      
      G_mothers_nontransmitted <- G_fathers_nontransmitted <- G <- matrix(nrow=n, ncol=sum_L)
      for(j in 1:sum_L) {
        G_mothers_j <- G_mothers[,j]/2
        indHet <- G_mothers_j==0.5
        G_mothers_j[indHet] <- rbinom(sum(indHet), 1, 0.5)
        G_mothers_nontransmitted_j <- G_mothers_j
        G_mothers_nontransmitted_j[indHet] <- 1-G_mothers_j[indHet]
        G_mothers_nontransmitted[,j] <- G_mothers_nontransmitted_j
        
        G_fathers_j <- G_fathers[,j]/2
        indHet <- G_fathers_j==0.5
        G_fathers_j[indHet] <- rbinom(sum(indHet), 1, 0.5)
        G_fathers_nontransmitted_j <- G_fathers_j
        G_fathers_nontransmitted_j[indHet] <- 1-G_fathers_j[indHet]
        G_fathers_nontransmitted[,j] <- G_fathers_nontransmitted_j
        
        G[,j] <- G_mothers_j + G_fathers_j
      }
      
      #Genetic variants associated with exposure
      G_X <- G[,1:L_X]
      colnames(G_X) <- G_X_names
      
      #Genetic variants associated with outcome
      G_Y <- G[,(L_X+1):sum_L]; rm(G)
      colnames(G_Y) <- G_Y_names
      
      #Combine the data in a single data frame
      data_offspring <- SimulateData(G_X, G_Y); rm(list=c('G_X', 'G_Y'))
      
      cur_zeta <- cur_zeta+1
      if(cur_zeta<=zeta) { data_mothers <- data_offspring }
    }
    
    #Calculate r2
    r2_res <- c(summary(lm(X~Z_X, data=data_offspring))$r.squared,
                summary(lm(Y~Z_Y, data=data_offspring))$r.squared)
    
    #Calculate variances
    G_X_matrix <- data_offspring[,G_X_names]
    Z_X        <- apply(G_X_matrix, 1, sum)
    
    G_Y_matrix <- data_offspring[,G_Y_names]
    Z_Y        <- apply(G_Y_matrix, 1, sum)
    
    var_res <- c(var(data_offspring$X), var(data_offspring$X_line), var(data_offspring$X_star),
                 var(data_offspring$Y), var(data_offspring$Y_line), var(data_offspring$Y_star),
                 var(data_offspring$C_line), var(data_offspring$C_star),
                 var(data_offspring$Z_X), var(data_offspring$Z_Y))
    
    #Calculate correlations
    G_X_cor       <- cor(G_X_matrix); rm(G_X_matrix)
    diag(G_X_cor) <- NA
    G_X_cor       <- mean(G_X_cor, na.rm=T)

    G_Y_cor       <- cor(G_Y_matrix); rm(G_Y_matrix)
    diag(G_Y_cor) <- NA
    G_Y_cor       <- mean(G_Y_cor, na.rm=T)
    
    cor_res <- c(cor(data_mothers$X,   data_fathers$X),
                 cor(data_mothers$Y,   data_fathers$Y),
                 cor(data_mothers$X,   data_fathers$Y),
                 cor(data_mothers$Y,   data_fathers$X),
                 cor(data_offspring$X, data_offspring$Y),
                 cor(Z_X, Z_Y), G_X_cor, G_Y_cor)
    
    data <- data_offspring; rm(data_offspring) 
    
    #########################
    #3.1.2.3) Compute result#
    #########################
    
    #Prepare the data for TSLS regression
    ind <- sample(1:n, n_X, replace=F)
    
    data_TSLS                       <- data[ind,]
    G_X_TSLS                        <- data_TSLS[,G_X_names] 
    phenotypes_TSLS                 <- data_TSLS[,c('X','Y')]
    G_X_mothers_TSLS                <- G_mothers[ind, G_X_names]; rm(G_mothers)
    G_X_fathers_TSLS                <- G_fathers[ind, G_X_names]; rm(G_fathers)
    G_X_mothers_nontransmitted_TSLS <- G_mothers_nontransmitted[ind, 1:L_X]; rm(G_mothers_nontransmitted)
    G_X_fathers_nontransmitted_TSLS <- G_fathers_nontransmitted[ind, 1:L_X]; rm(G_fathers_nontransmitted); rm(data_TSLS)
    X_parents_TSLS                  <- data.frame(X_mothers=data_mothers$X, X_fathers=data_fathers$X); rm(list=c('data_mothers', 'data_fathers'))
    X_parents_TSLS                  <- X_parents_TSLS[ind,]
    
    #Split the data into two non-overlapping datasets    
    data_X    <- data[1:n_X    , c(G_X_names, 'X')]
    data_Y    <- data[(n_X+1):n, c(G_X_names, 'Y')]; rm(data)

    #Sample H genetic variants from the L_X exposure-associated variants to use as genetic instruments
    IV_index <- sample(1:L_X, H, F)
    
    G_X_TSLS                        <- G_X_TSLS[,IV_index]
    G_X_mothers_TSLS                <- G_X_mothers_TSLS[,IV_index]
    G_X_fathers_TSLS                <- G_X_fathers_TSLS[,IV_index]
    G_X_mothers_nontransmitted_TSLS <- G_X_mothers_nontransmitted_TSLS[,IV_index]
    G_X_fathers_nontransmitted_TSLS <- G_X_fathers_nontransmitted_TSLS[,IV_index]
    
    data_X <- data_X[,c(IV_index, ncol(data_X) )]
    data_Y <- data_Y[,c(IV_index, ncol(data_Y) )]
    
    #Compute summary associations
    BetaXG <- BetaYG <- seBetaXG <- seBetaYG <- numeric(H)
    
    for(h in 1:H) {
      #Instrument-exposure associations
      fit_X <- lm(data_X$X~data_X[,h])
      BetaXG[h]   <- coef(fit_X)[2]
      seBetaXG[h] <- summary(fit_X)$coef[2,2]
      
      #Instrument-outcome associations      
      fit_Y <- lm(data_Y$Y~data_Y[,h])
      BetaYG[h]   <- coef(fit_Y)[2]
      seBetaYG[h] <- summary(fit_Y)$coef[2,2]
    }
    
    #Mendelian randomization
    TSLS_M1_res         <- TSLS(model=1, weights=BetaXG,
                                X=phenotypes_TSLS$X, Y=phenotypes_TSLS$Y, G_X=G_X_TSLS)
    TSLS_M2_res         <- TSLS(model=2, weights=BetaXG,
                                X=phenotypes_TSLS$X, Y=phenotypes_TSLS$Y, G_X=G_X_TSLS,
                                G_X_mothers=G_X_mothers_TSLS, G_X_fathers=G_X_fathers_TSLS)
    TSLS_M3_res         <- TSLS(model=3, skip_fathers=omega==6 & P==1, weights=BetaXG,
                                Y=phenotypes_TSLS$Y,
                                X_mothers=X_parents_TSLS$X_mothers, G_X_mothers=G_X_mothers_nontransmitted_TSLS,
                                X_fathers=X_parents_TSLS$X_fathers, G_X_fathers=G_X_fathers_nontransmitted_TSLS)
    TSLS_M4_res         <- TSLS(model=4, skip_fathers=omega==6 & P==1, weights=BetaXG,
                                X=phenotypes_TSLS$X, Y=phenotypes_TSLS$Y, G_X=G_X_TSLS,
                                X_mothers=X_parents_TSLS$X_mothers, G_X_mothers=G_X_mothers_nontransmitted_TSLS,
                                X_fathers=X_parents_TSLS$X_fathers, G_X_fathers=G_X_fathers_nontransmitted_TSLS)
    IVW_res            <- IVW(BetaXG, BetaYG, seBetaYG)
    MREgger_res        <- MREgger(BetaXG, BetaYG, seBetaYG)
    WeightedMedian_res <- WeightedMedian(BetaXG, BetaYG, seBetaXG, seBetaYG)
    ModeEstimator_res  <- as.numeric(MBE(BetaXG, BetaYG, seBetaXG, seBetaYG))
    
    MR_res_temp <- rbind(TSLS_M1_res[-5], TSLS_M2_res[,-5], TSLS_M3_res[,-5], TSLS_M4_res[,-5],
                         IVW_res[-5], MREgger_res[2,-5], WeightedMedian_res[-5], ModeEstimator_res)
    MR_res <- NULL
    for(i in 1:nrow(MR_res_temp)) {
      MR_res <- c(MR_res,
                  c(MR_res_temp[i,1:2], coverage_power(MR_res_temp[i,3], MR_res_temp[i,4], beta)))
    }; rm(MR_res_temp)

    NOME_res <- NOME(BetaXG, seBetaXG, seBetaYG) #Computing instrument strength statistics
    
    #Horizontal pleiotropy tests
    MR_Egger_Int  <- c(MREgger_res[1,1], MREgger_res[1,3]>=0 | MREgger_res[1,4]<=0) 
    Cochran_test  <- HeterogeneityTest(BetaXG=BetaXG, BetaYG=BetaYG, seBetaXG=seBetaXG, seBetaYG=seBetaYG)
    Cochran_Q     <- Cochran_test[1] 
    Cochran_power <- Cochran_test[3]<0.05
    
    c(r2_res, var_res, cor_res, MR_res, NOME_res, MR_Egger_Int, Cochran_Q, Cochran_power)
  }
  
  ########################################################################
  #3.1.3) Finalise the analysis for the current combination of parameters#
  ########################################################################
  
  simulation_results <- apply(simulation_results, 2, mean)                #Average the results across all simulated datasets
  simulation_results <- c(as.numeric(parameters[a,]), simulation_results) #Add the parameters corresponding to the current simulations
  
  results           <- rbind(results, simulation_results) #Add the current results to the previous results
  colnames(results) <- names_simulation_results 
  write.table(results, filename, row.names=F, sep='\t')
}
