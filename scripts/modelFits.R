################################################################################
# function to fit models 
# inputs: 
#   incData - observed incidence data
#   alarmFit - type of function to use to describe the alarm
################################################################################

# Sourced and adapted from: Ward, C., Deardon, R., and Schmidt, A. M. (2023).
#   Bayesian modeling of dynamic behavioral change during an epidemic.
#   Infectious Disease Modelling, 8(4):947–963.

################################################################################

fitAlarmModel <- function(incData, smoothI, N, I0, R0, prior, smoothWindow, 
                          alarmFit, seed) {
  
  source('./scripts/modelCodes.R')
  source('./scripts/getModelInputs.R')
  
  # Load appropriate model code (SIR_power)
  modelCode <- get(paste0('SIR_', alarmFit))
  
  # Set seed for reproducibility
  set.seed(seed + 3)
  
  # Prepare model inputs
  modelInputs <- getModelInput(
    alarmFit = alarmFit,
    incData = incData,
    smoothI = smoothI,
    prior = prior,
    smoothWindow = smoothWindow,
    N = N,
    I0 = I0,
    R0 = R0
  )
  
  niter <- modelInputs$niter
  nburn <- modelInputs$nburn
  nthin <- modelInputs$nthin
  
  # Create and configure model
  myModel <- nimbleModel(modelCode, 
                         data = modelInputs$dataList, 
                         constants = modelInputs$constantsList,
                         inits = modelInputs$initsList)
  
  myConfig <- configureMCMC(myModel)
  
  # Monitor alarm components
  myConfig$addMonitors(c('yAlarm', 'alarm'))
  
  # Replace sampler for Rstar with custom update
  myConfig$removeSamplers('Rstar')
  myConfig$addSampler(target = 'Rstar', type = "RstarUpdate")
  
  # Model-specific samplers
  if (alarmFit == 'power') {
    myConfig$removeSampler(c('beta', 'rateI'))
    myConfig$addSampler(target = c('beta', 'rateI'), type = "AF_slice")
    
    myConfig$removeSampler('k')
    myConfig$addSampler(target = 'k', type = "slice")
    
  } 
  
  myConfig$addMonitors(c('Rstar'))
  
  print(myConfig)
  
  nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
  myMCMC <- buildMCMC(myConfig)
  compiled <- compileNimble(myModel, myMCMC)
  
  runMCMC(compiled$myMCMC, 
          niter = niter, 
          nburnin = nburn, 
          thin = nthin,
          setSeed = seed)
}



