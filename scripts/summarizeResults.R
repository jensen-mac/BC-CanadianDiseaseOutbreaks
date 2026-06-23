################################################################################
# Take raw results files for each outbreak and summarize the key results into
# smaller files to be used for figures/tables.
################################################################################

# sourcing all result file paths
results_files <- list.files("results", pattern = "^output_power_.*\\.rds$", full.names = TRUE)

dir.create("summarized_results", showWarnings = FALSE)

results <- map(results_files, function(f) {
  post <- readRDS(f)
  combined_samples <- do.call(rbind, post$resThree)
  
  # Basic reproductive number = beta/gamma
  R0 <- combined_samples[, "beta"] / combined_samples[, "rateI"]
  
  alarm_cols <- grep("^alarm\\[", colnames(combined_samples))
  # maximum alarm
  max_alarms <- apply(combined_samples[, alarm_cols], 1, max)
  
  # posterior predictive fit
  n <- nrow(post$postPredFit)
  ppf <- post$postPredFit %>%
    mutate(disease  = post$disease,
           year = post$year,
           removals = as.numeric(post$removals),
           date = tail(post$dates, n),
           obs = tail(post$incData, n))
  
  list(
    params = data.frame(disease = post$disease, year = post$year,
                        removals = as.numeric(post$removals), post$postParams),
    alarm = data.frame(disease = post$disease, year = post$year,
                        removals = as.numeric(post$removals),
                        mean = mean(max_alarms),
                        lower = quantile(max_alarms, 0.025),
                        upper = quantile(max_alarms, 0.975)),
    R0 = data.frame(disease = post$disease, year = post$year,
                        removals = as.numeric(post$removals), param = "R0",
                        mean = mean(R0),
                        lower = quantile(R0, 0.025),
                        upper = quantile(R0, 0.975)),
    ppf = ppf
  )
})

saveRDS(map_dfr(results, "params"), "summarized_results/params.rds")
saveRDS(map_dfr(results, "alarm"), "summarized_results/alarm.rds")
saveRDS(map_dfr(results, "R0"), "summarized_results/R0.rds")
saveRDS(map_dfr(results, "ppf"), "summarized_results/ppf.rds")

