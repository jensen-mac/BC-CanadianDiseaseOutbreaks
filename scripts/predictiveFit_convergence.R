################################################################################
# Save posterior predictive fit and convergence diagnostics to fits folder
# Used to assess model fit for all years
################################################################################

library(ggplot2)
library(coda)

dir.create("fits", showWarnings = FALSE)

# sourcing all result file paths
results_files <- list.files("results", pattern = "^output_power_.*\\.rds$", full.names = TRUE)

all_gdiag <- list()

for (f in results_files) {
  
  file_label <- tools::file_path_sans_ext(basename(f))
  
  postSummaries <- readRDS(f)
  
  postPredFit <- postSummaries$postPredFit
  gdiag <- postSummaries$gdiag
  dates <- postSummaries$dates
  incData <- postSummaries$incData
  resThree <- postSummaries$resThree
  
  # Posterior predictive fit
  postPredFit$date <- tail(dates,   nrow(postPredFit))
  postPredFit$obs <- tail(incData, nrow(postPredFit))
  
  predictive_fit <- ggplot(postPredFit, aes(x = date)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = adjustcolor("steelblue", alpha.f = 0.3)) +
    geom_line(aes(y = mean), color = "steelblue", linetype = "dashed", linewidth = 0.7) +
    geom_line(aes(y = obs),  color = "black", linewidth = 0.7) +
    scale_x_date(date_labels = "%b", date_breaks = "2 months") +
    labs(x = "Date", y = "Incidence", title = file_label) +
    theme_minimal(base_size = 15) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 15)
    )
  
  ggsave(
    filename = file.path("fits", paste0(file_label, "_postpredfit.png")),
    plot = predictive_fit, width = 8, height = 4, dpi = 300
  )
  
  # Gelman-Rubin diagnostic values
  ggdiag_df <- gdiag
  gdiag_df$file <- file_label
  all_gdiag[[file_label]] <- gdiag_df
  
  # Traceplots
  power_params <- c("beta", "rateI", "k")
  available_params <- colnames(resThree[[1]])
  power_params <- power_params[power_params %in% available_params]
  
  mcmc_list <- lapply(resThree, function(res) as.mcmc(res[, power_params, drop = FALSE]))
  combined_mcmc <- mcmc.list(mcmc_list)
  
  png(
    filename = file.path("fits", paste0(file_label, "_traceplots.png")),
    width = 12, height = 10, units = "in", res = 300       
  )
  par(mfrow = c(2, 2))
  for (param in power_params) {
    traceplot(combined_mcmc[, param], main = paste(param, ",", file_label))
  }
  par(mfrow = c(1, 1))
  dev.off()
  
}

gdiag_all <- do.call(rbind, all_gdiag)
rownames(gdiag_all) <- NULL
gdiag_all <- gdiag_all[, c("file", "param", "gr", "grUpper")]

write.csv(gdiag_all, file.path("fits", "gdiag_summary.csv"), row.names = FALSE)

