################################################################################
# summarize posterior samples from three chains
# gelman-rubin to assess convergence
# summarize posterior means and credible intervals:
#   model parameters   
#   alarm function
# posterior prediction mean and credible intervals
################################################################################

# Sourced and adapted from: Ward, C., Deardon, R., and Schmidt, A. M. (2023).
#   Bayesian modeling of dynamic behavioral change during an epidemic.
#   Infectious Disease Modelling, 8(4):947–963.

################################################################################

library(coda)
library(nimble)

# source relevant scripts
source('./scripts/modelCodes.R')
source('./scripts/getModelInputs.R')
source('./scripts/postPredFit.R')

summarizePost <- function(resThree, incData, smoothI, smoothWindow,
                          N, I0, R0, Istar0, Rstar0, lengthI, 
                          alarmFit, prior, disease, year, dates) {
  
  if (alarmFit == 'power') {
    paramSamples1 <- resThree[[1]][,-grep('alarm|yAlarm|Rstar|R0',
                                          colnames(resThree[[1]]))]
    paramSamples2 <- resThree[[2]][,-grep('alarm|yAlarm|Rstar|R0', 
                                          colnames(resThree[[2]]))]
    paramSamples3 <- resThree[[3]][,-grep('alarm|yAlarm|Rstar|R0', 
                                          colnames(resThree[[3]]))]
  } else if (alarmFit == 'basic') {
    paramSamples1 <- resThree[[1]][,-grep('Rstar|R0', 
                                          colnames(resThree[[1]])), drop = F]
    paramSamples2 <- resThree[[2]][,-grep('Rstar|R0', 
                                          colnames(resThree[[2]])), drop = F]
    paramSamples3 <- resThree[[3]][,-grep('Rstar|R0', 
                                          colnames(resThree[[3]])), drop = F]
  }
  
  # Rstar Posterior
  RstarSamples1 <- resThree[[1]][,grep('Rstar', colnames(resThree[[1]]))]
  RstarSamples2 <- resThree[[2]][,grep('Rstar', colnames(resThree[[2]]))]
  RstarSamples3 <- resThree[[3]][,grep('Rstar', colnames(resThree[[3]]))]
  RstarPost <- rbind(RstarSamples1, RstarSamples2, RstarSamples3)
  
  ##############################################################################
  ### gelman-rubin
  res_mcmc <- mcmc.list(mcmc(paramSamples1),
                        mcmc(paramSamples2),
                        mcmc(paramSamples3))
  gdiag <- data.frame(gelman.diag(res_mcmc, multivariate = F)$psrf)
  colnames(gdiag) <- c('gr', 'grUpper')
  gdiag$param <- rownames(gdiag)
  rownames(gdiag) <- NULL
  
  ##############################################################################
  ### posterior mean and 95% CI for parameters
  paramsPost <- rbind(paramSamples1, paramSamples2, paramSamples3)
  postMeans <- colMeans(paramsPost)
  postCI <- apply(paramsPost, 2, quantile, probs = c(0.025, 0.975))
  postParams <- data.frame(param = names(postMeans),
                           mean = postMeans,
                           lower = postCI[1,],
                           upper = postCI[2,])
  rownames(postParams) <- NULL
  
  ##############################################################################
  ### posterior distribution of alarm function 
  
  if (alarmFit == 'power') {
    alarmSamples1 <- t(resThree[[1]][,grep('yAlarm', colnames(resThree[[1]]))])
    alarmSamples2 <- t(resThree[[2]][,grep('yAlarm', colnames(resThree[[2]]))])
    alarmSamples3 <- t(resThree[[3]][,grep('yAlarm', colnames(resThree[[3]]))])
    alarmSamples <- cbind(alarmSamples1, alarmSamples2, alarmSamples3)
    
    postMeans <- rowMeans(alarmSamples)
    postCI <- apply(alarmSamples, 1, quantile, probs = c(0.025, 0.975))
    
    n <- 50
    maxI <- ceiling(max(smoothI))
    xAlarm <- seq(0, maxI, length.out = n)
    
    postAlarm <- data.frame(xAlarm = xAlarm, 
                            mean = postMeans,
                            lower = postCI[1,],
                            upper = postCI[2,])
    
    rownames(postAlarm) <- NULL
    
  } else {
    postAlarm <- data.frame(xAlarm = NA, 
                            mean = NA,
                            lower = NA,
                            upper = NA)
  }
  
  
  ##############################################################################
  ### posterior predictive model fit
  
  postPredObs <- postPredFit(incData = incData, smoothI = smoothI,
                             N = N, I0 = I0, R0 = R0, Istar0 = Istar0,
                             Rstar0 = Rstar0, lengthI = lengthI, 
                             alarmFit = alarmFit, prior = prior, 
                             smoothWindow = smoothWindow, 
                             paramsPost = paramsPost, alarmSamples = alarmSamples)
  
  postMean <- rowMeans(postPredObs)
  postCI <- apply(postPredObs, 1, quantile, probs = c(0.025, 0.975))
  postPredFit <- data.frame(time = 1:length(incData),
                            mean = postMean,
                            lower = postCI[1,],
                            upper = postCI[2,])
  
  ##############################################################################
  
  ### output
  out <- list(gdiag = gdiag,
              postParams = postParams,
              postAlarm = postAlarm,
              postPredFit = postPredFit,
              N = N, 
              resThree = resThree, 
              dates = dates, 
              prior = prior, 
              incData = incData, 
              disease = disease,
              year = year,
              removals = R0) 
  
  # Save the output to results folder
  filename <- paste0("results/output_", alarmFit, "_", disease, "_", year, "_",
                     R0,".rds")
  saveRDS(out, file = filename)
  
  return(out)
  
}








