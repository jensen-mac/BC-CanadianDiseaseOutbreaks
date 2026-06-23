################################################################################
# Date Selection
# Use to identify the start/end dates for outbreaks to be analyzed, label years
# Manual action required is identified in comments with TODO
################################################################################

# Required packages
library(ggplot2)

# Loading disease and population data
influenza <- read.csv("./data/influenza_CA_ON_weekly.csv")
measles   <- read.csv("./data/measles_CA_ON_weekly.csv")

influenza$period_end_date   <- as.Date(influenza$period_end_date)
influenza$period_start_date <- as.Date(influenza$period_start_date)
measles$period_end_date     <- as.Date(measles$period_end_date)
measles$period_start_date   <- as.Date(measles$period_start_date)

# Initialize tables for start and end dates of each outbreak for diseases
if (!exists("OutbreakDates_influenza")) {
  OutbreakDates_influenza <- data.frame(
    year  = character(),
    start = character(),
    end   = character(),
    stringsAsFactors = FALSE
  )
}

if (!exists("OutbreakDates_measles")) {
  OutbreakDates_measles <- data.frame(
    year  = character(),
    start = character(),
    end   = character(),
    stringsAsFactors = FALSE
  )
}

########################### 1: Analyze each season #############################
disease <- "influenza" # TODO: change <disease> to either "influenza" or "measles"
df <- get(disease)
yr <- "1976"   # TODO: change <yr> for each year analyzed

# most outbreaks start in the fall
# classifies the year of the outbreak according to the fall, 
# even if the outbreak starts in the winter, to avoid duplicate analysis years
yr_data <- df[df$period_start_date >= as.Date(paste0(yr,"-09-01")) &
                df$period_end_date   <= as.Date(paste0(as.integer(yr) + 1, "-08-31")), ]

ggplot(yr_data, aes(x = period_end_date, y = cases_this_period)) +
  geom_line() +
  labs(title   = paste(yr, disease, "full season"),
       x = "Date", y = "Cases") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA))

####################### 2: Manually adjust for outbreak ########################
# TODO: adjust <year_start> and <year_end> for that year's outbreak
# Start: first week of sustained increase toward peak, avoid leading plateaus
# End: when the overall epidemic curve has declined from peak,
#      does not have to go down to zero cases

year_start <- "1977-01-22"   
year_end   <- "1977-04-02"   

outbreak_data <- df[df$period_end_date >= as.Date(year_start) &
                      df$period_end_date <= as.Date(year_end), ]

ggplot(outbreak_data, aes(x = period_end_date, y = cases_this_period)) +
  geom_line() +
  labs(title    = paste(yr, disease, "outbreak"),
       subtitle  = paste("Start:", year_start, "  End:", year_end),
       x = "Date", y = "Cases") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA))

outbreak_data$cases_this_period

###################### 3: Add that defined season to table #####################

tableName <- paste0("OutbreakDates_", disease)
OutbreakDates <- get(tableName)

if (!yr %in% OutbreakDates$year) {
  OutbreakDates <- rbind(OutbreakDates, data.frame(
    year  = yr,
    start = year_start,
    end   = year_end,
    stringsAsFactors = FALSE
  ))
  assign(tableName, OutbreakDates)
  cat("Added", disease, "season", yr)
} else {
  cat(disease, "season", yr, "already in table")
}

# Save after every year
write.csv(OutbreakDates,
          file.path("data", paste0("OutbreakDates_", disease, ".csv")),
          row.names = FALSE)




























