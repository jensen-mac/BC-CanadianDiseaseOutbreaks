################################################################################
# MCMC Fitting
# Run all years for both diseases
# Run only a selected year (e.g., new fitting, changing R0)
# Manual action required is identified in comments with TODO
################################################################################

# Required Packages
library(nimble)
library(parallel)

nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)

# Source scripts
source('./scripts/modelCodes.R')
source('./scripts/getModelInputs.R')
source('./scripts/modelFits.R')
source('./scripts/summarizePost.R')
source('./scripts/postPredFit.R')

# Load data
influenza <- read.csv("data/influenza_CA_ON_weekly.csv")
measles <- read.csv("data/measles_CA_ON_weekly.csv")

############################# Run all models ###################################
# Runs the model for all outbreak years and both diseases

diseases <- c("influenza", "measles")

disease_params <- list(
  influenza = list(smoothWindow = 4,  prior = 2),
  measles = list(smoothWindow = 16, prior = 1)
)

lengthI <- 1
idxStart <- 2
alarmFit_i <- "power"

for (disease in diseases) {
  
  outbreakDates <- read.csv(paste0("data/OutbreakDates_", disease, ".csv"))
  outbreakDates$start <- as.Date(outbreakDates$start)
  outbreakDates$end <- as.Date(outbreakDates$end)
  
  smoothWindow_i <- disease_params[[disease]]$smoothWindow
  prior_i <- disease_params[[disease]]$prior
  
  dat <- get(disease)
  dat$period_start_date <- as.Date(dat$period_start_date)
  dat$period_end_date <- as.Date(dat$period_end_date)
  
  for (i in seq_len(nrow(outbreakDates))) {
    
    start_yr <- outbreakDates$start[i]
    end_yr <- outbreakDates$end[i]
    year <- outbreakDates$year[i]
    
    # Define outbreak window according to pre-defined dates
    disease_yr <- dat[
      dat$period_end_date >= start_yr &
        dat$period_end_date <= end_yr, ]
    
    N <- disease_yr$population[1]
    disease_yr$smoothedCases <- round(movingAverage(disease_yr$cases_this_period, 2))
    
    # Initial states, 0 removals assumption
    I0 <- sum(disease_yr$smoothedCases[max(1, idxStart - lengthI + 1):idxStart])
    R0 <- 0
    
    # smoothed incidence to inform alarm function 
    # (shifted so alarm is informed only by data up to time t-1)
    disease_yr$smoothI <- head(
      movingAverage(c(I0, disease_yr$smoothedCases), smoothWindow_i),
      -1
    )
    
    # Initialize current number of infectious and removed individuals
    incData <- disease_yr$smoothedCases[-1]
    smoothI <- disease_yr$smoothI[-1]
    
    # Initialize removal vector and previously observed incidence
    Rstar0 <- disease_yr$smoothedCases[max(1, idxStart - lengthI + 1):idxStart]
    Istar0 <- disease_yr$smoothedCases[max(1, idxStart - smoothWindow_i + 1):idxStart]
    
    # Run three chains in parallel
    cl <- makeCluster(3)
    clusterExport(cl, list('incData', 'smoothI', 'alarmFit_i', 'prior_i',
                           'smoothWindow_i', 'N', 'I0', 'R0', 'Rstar0',
                           'lengthI'))
    
    resThree <- parLapplyLB(cl, 1:3, function(x) {
      library(nimble)
      nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
      source('./scripts/modelFits.R')
      
      fitAlarmModel(incData = incData,
                    smoothI = smoothI,
                    N = N,
                    I0 = I0,
                    R0 = R0,
                    prior = prior_i,
                    smoothWindow = smoothWindow_i,
                    alarmFit = alarmFit_i,
                    seed = x)
    })
    
    stopCluster(cl)
    
    postSummaries <- summarizePost(
      resThree = resThree,
      incData = incData,
      smoothI = smoothI,
      smoothWindow = smoothWindow_i,
      N = N,
      I0 = I0,
      R0 = R0,
      Istar0 = Istar0,
      Rstar0 = Rstar0,
      lengthI = lengthI,
      alarmFit = alarmFit_i,
      prior = prior_i,
      disease = disease,
      year = year,
      dates = disease_yr$period_end_date
    )
  }
}

############################### Run one year ###################################
# Used to run model after one-off adjustments of start/end dates for outbreak years
# Used for running models on selected years (1945, 1960, 1979) for Section 3.2

# TODO: set for desired <disease> and <year> 
disease <- "measles"
year <- 1945
R0 <- 0 # TODO: change removal assumptions

# Parameters
disease_params <- list(
  influenza = list(smoothWindow = 4,  prior = 2),
  measles = list(smoothWindow = 16, prior = 1)
)

lengthI <- 1
idxStart <- 2
alarmFit_i <- "power"

smoothWindow_i <- disease_params[[disease]]$smoothWindow
prior_i <- disease_params[[disease]]$prior

outbreakDates <- read.csv(paste0("data/OutbreakDates_", disease, ".csv"))
outbreakDates$start <- as.Date(outbreakDates$start)
outbreakDates$end <- as.Date(outbreakDates$end)

start_yr <- outbreakDates$start[outbreakDates$year == year]
end_yr <- outbreakDates$end[outbreakDates$year == year]

dat <- get(disease)
dat$period_start_date <- as.Date(dat$period_start_date)
dat$period_end_date <- as.Date(dat$period_end_date)

# Define outbreak window according to pre-defined dates
disease_yr <- dat[
  dat$period_end_date >= start_yr &
    dat$period_end_date <= end_yr, ]

N <- disease_yr$population[1]
disease_yr$smoothedCases <- round(movingAverage(disease_yr$cases_this_period, 2))

# Initial states,
I0 <- sum(disease_yr$smoothedCases[max(1, idxStart - lengthI + 1):idxStart])

# Smoothed incidence to inform alarm function
# (shifted so alarm is informed only by data up to time t-1)
disease_yr$smoothI <- head(
  movingAverage(c(I0, disease_yr$smoothedCases), smoothWindow_i),
  -1
)

# Initialize current number of infectious and removed individuals
incData <- disease_yr$smoothedCases[-1]
smoothI <- disease_yr$smoothI[-1]

# Initialize removal vector and previously observed incidence
Rstar0 <- disease_yr$smoothedCases[max(1, idxStart - lengthI + 1):idxStart]
Istar0 <- disease_yr$smoothedCases[max(1, idxStart - smoothWindow_i + 1):idxStart]

# Run three chains in parallel
cl <- makeCluster(3)
clusterExport(cl, list('incData', 'smoothI', 'alarmFit_i', 'prior_i',
                       'smoothWindow_i', 'N', 'I0', 'R0', 'Rstar0', 
                       'lengthI'))

resThree <- parLapplyLB(cl, 1:3, function(x) {
  library(nimble)
  nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
  source('./scripts/modelFits.R')
  
  fitAlarmModel(incData = incData,
                smoothI = smoothI,
                N = N,
                I0 = I0,
                R0 = R0,
                prior = prior_i,
                smoothWindow = smoothWindow_i,
                alarmFit = alarmFit_i,
                seed = x)
})

stopCluster(cl)

postSummaries <- summarizePost(
  resThree = resThree,
  incData = incData,
  smoothI = smoothI,
  smoothWindow = smoothWindow_i,
  N = N,
  I0 = I0,
  R0 = R0,
  Istar0 = Istar0,
  Rstar0 = Rstar0,
  lengthI = lengthI,
  alarmFit = alarmFit_i,
  prior = prior_i,
  disease = disease,
  year = year,
  dates = disease_yr$period_end_date
)







