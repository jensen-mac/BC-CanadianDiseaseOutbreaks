################################################################################
# Get model inputs for nimble for power alarm model
# used for posterior prediction
################################################################################

# Sourced and adapted from: Ward, C., Deardon, R., and Schmidt, A. M. (2023).
#   Bayesian modeling of dynamic behavioral change during an epidemic.
#   Infectious Disease Modelling, 8(4):947–963.

################################################################################

getModelInput <- function(alarmFit, incData, smoothI, prior, smoothWindow,
                          N, I0, R0) {
  
  # Initial number of susceptibles
  S0 <- N - I0 - R0
  tau <- length(incData)
  
  # Shared data for most models
  dataList <- list(Istar = incData, smoothI = smoothI)
  Rstar <- c(I0, head(incData, -1))
  
  # Prior parameters for infectious period rate (gamma distribution)
  if (prior == 1) {
    # Centered on 8 days, 90% interval: 7-9 days
    # Used for measles
    # Data is in weekly format, so priors reflect this weekly format
    aa <- 172
    bb <- 194
  } else if (prior == 2) {
    # Centered on 5 days, 90% interval: 4-6 days
    # Used for influenza
    # Data is in weekly format, so priors reflect this weekly format
    aa <- 45.8
    bb <- 66.3
  } 
  
  # Handle different alarm models
  if (alarmFit == 'power') {
    maxI <- ceiling(max(smoothI))
    n <- 50
    xAlarm <- seq(0, maxI, length.out = n)
    
    constantsList <- list(tau = tau, N = N, I0 = I0, R0 = R0,
                          n = n, xAlarm = xAlarm, aa = aa, bb = bb)
    
    initsList <- list(beta = runif(1, 1/7, 1),
                      k = runif(1, 0, 1),
                      rateI = rgamma(1, aa, bb),
                      Rstar = Rstar)
    
  } else if (alarmFit == 'basic') {
    # No alarm-related constants needed
    constantsList <- list(tau = tau, N = N, I0 = I0, R0 = R0,
                          aa = aa, bb = bb)
    
    dataList <- list(Istar = incData)  # no smoothI needed
    
    initsList <- list(beta = runif(1, 1/7, 1),
                      rateI = rgamma(1, aa, bb),
                      Rstar = Rstar)
    
    xAlarm <- NULL
    
  }
  
  niter <- 300000
  nburn <- 100000
  nthin <- 10
  
  return(list(constantsList = constantsList,
              dataList = dataList,
              initsList = initsList,
              niter = niter,
              nburn = nburn,
              nthin = nthin,
              xAlarm = xAlarm))
}

