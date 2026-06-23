################################################################################
# Obtain removal estimates for sensitivity scenarios 2 & 3
# 2: individuals >15 years have previously contracted measles and are immune
# 3: 95% of the population is removed

# Download population data from Statistics Canada:
# https://open.canada.ca/data/en/dataset/04e5b2c5-219b-408c-b9ea-4da828ccab81
# https://open.canada.ca/data/en/dataset/ecdee020-5919-4996-8d3d-c3df75f50ca0 
################################################################################

library(dplyr)

# Load data for N estimates
influenza <- read.csv("data/influenza_CA_ON_weekly.csv")
measles <- read.csv("data/measles_CA_ON_weekly.csv")

# Population estimates by age (Statistics Canada, 2000, 2024)
ON_pop_1921_1971 <- read.csv("data/17100029.csv")
ON_pop_1971_2025 <- read.csv("data/17100005.csv") 

head(ON_pop_1921_1971)
head(ON_pop_1971_2025)

############################## Scenario 2 ######################################

# population units are in thousands
age_0_15_1945and1960 <- ON_pop_1921_1971 %>% mutate(Population = VALUE * 1000) %>%
  filter(GEO == "Ontario",
         REF_DATE %in% c(1945, 1960),
         Age.group == "0-15 years",
         Sex == "Both sexes") %>%
  select(REF_DATE, GEO, Age.group, Population)

# total population estimates are stored in results from 0 removal scenario
results_1945 <- readRDS("results/output_power_measles_1945_0.rds")
population_1945 <- results_1945$N
removals_1945_scenario2 <- population_1945 - age_0_15_1945and1960$Population[1]
removals_1945_scenario2

results_1960 <- readRDS("results/output_power_measles_1960_0.rds")
population_1960 <- results_1960$N
removals_1960_scenario2 <- population_1960 - age_0_15_1945and1960$Population[2]
removals_1960_scenario2

age_0_15_1979 <- ON_pop_1971_2025  %>%
  filter(GEO == "Ontario",
         REF_DATE %in% c(1979),
         Age.group == "0 to 15 years",
         Gender %in% c("Total - gender")) %>%
  select(REF_DATE, GEO, Age.group, Population = VALUE)

results_1979 <- readRDS("results/output_power_measles_1979_0.rds")
population_1979 <- results_1979$N
removals_1979_scenario2 <- population_1979 - age_0_15_1979$Population[1]
removals_1979_scenario2

############################## Scenario 3 ######################################

removals_1945_scenario3 <- results_1945$N*0.95
removals_1945_scenario3

removals_1960_scenario3 <- results_1960$N*0.95
removals_1960_scenario3

removals_1979_scenario3 <- results_1979$N*0.95
removals_1979_scenario3










