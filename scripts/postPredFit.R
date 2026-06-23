################################################################################
# function to do posterior predictive distribution of observed epidemic curve 
################################################################################

# Sourced and adapted from: Ward, C., Deardon, R., and Schmidt, A. M. (2023).
#   Bayesian modeling of dynamic behavioral change during an epidemic.
#   Infectious Disease Modelling, 8(4):947–963.

################################################################################

postPredFit <- function(incData, smoothI, N, I0, R0, Istar0, Rstar0, lengthI,
                        alarmFit, prior, smoothWindow, 
                        paramsPost, alarmSamples) {
  
  # model-specific constants, data, and inits
  modelInputs <- getModelInput(alarmFit = alarmFit, 
                               incData = incData, smoothI = smoothI, 
                               prior = prior, 
                               smoothWindow = smoothWindow,
                               N = N, I0 = I0, R0 = R0)
  
  # model code
  if (alarmFit == 'power') {
    modelCode <- get(paste0('SIR_', alarmFit, '_sim'))
    modelInputs$constantsList$bw <- smoothWindow
  } else {
    modelCode <- get(paste0('SIR_', alarmFit))
  }
  
  # compile model and simulator
  myModelPred <- nimbleModel(modelCode, 
                             constants = modelInputs$constantsList)
  
  compiledPred <- compileNimble(myModelPred) 
  
  tau <- modelInputs$constantsList$tau
  dataNodes <- paste0('Istar[', 1:tau, ']')
  dataNodes <- c(dataNodes, paste0('Rstar[', 1:tau, ']'))
  
  sim_R <- simulator(myModelPred, dataNodes)
  sim_C <- compileNimble(sim_R)
  
  # get order of parameters
  parentNodes <- myModelPred$getParents(dataNodes, stochOnly = TRUE)
  parentNodes <- parentNodes[-which(parentNodes %in% dataNodes)]
  parentNodes <- myModelPred$expandNodeNames(parentNodes, returnScalarComponents = TRUE)
  
  nPost <- 10000
  postPredInc <- matrix(NA, nrow = tau, ncol = nPost)
  set.seed(1)
  for (j in 1:nPost) {
    
    postIdx <- sample(1:nrow(paramsPost), 1)
    
    betaPost <- paramsPost[postIdx,'beta']
    
    # model specific parameters
    if (alarmFit == 'power') {
      alarmParamPost <- paramsPost[postIdx, 'k']
      trueVals <- c(betaPost, alarmParamPost)
    } 
    
    # for exponential infectious period
    rateIPost <- paramsPost[postIdx, 'rateI']
    trueVals <- c(trueVals, rateIPost)
  
    trueVals <- trueVals[parentNodes]
    
    postPredInc[,j] <- apply(sim_C$run(trueVals, 10), 2, median)[grep('Istar', dataNodes)]
  }
  
  postPredInc
}
