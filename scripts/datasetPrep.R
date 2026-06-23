################################################################################
# Dataset Preparation
# Extracts desired disease datasets, following initial data exploration
################################################################################

# Required Packages
library(dplyr)

dir.create("data", showWarnings = FALSE)

# CANDID Normalized Dataset
# Obtained from https://github.com/canmod/iidda
canmod <- read.csv("canmod-cdi-normalized.csv")
head(canmod)

# Weekly influenza cases, in Ontario, Canada
influenza <- canmod %>%
  filter(
    iso_3166_2 == "CA-ON",
    time_scale == "wk",
    disease == "influenza"
  )
write.csv(influenza, file.path("data", "influenza_CA_ON_weekly.csv"), row.names = FALSE)

# Weekly measles cases, in Ontario, Canada
measles <- canmod %>%
  filter(
    iso_3166_2 == "CA-ON",
    time_scale == "wk",
    disease == "measles"
  )
write.csv(influenza, file.path("data", "measles_CA_ON_weekly.csv"), row.names = FALSE)













