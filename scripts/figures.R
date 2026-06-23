################################################################################
# Creating Figures for paper
################################################################################

library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(grid)
library(scales)

################################################################################
# Loading in data from summarized results folder used for figures
################################################################################

# Summarized estimates obtained through summarizeResults.R
parameter_estimates <- readRDS("summarized_results/params.rds")
alarm <- readRDS("summarized_results/alarm.rds")
reproductive_number <- readRDS("summarized_results/R0.rds")
predicted_fit <- readRDS("summarized_results/ppf.rds")

# Figure 1 and 2 
no_removals_parameters <- parameter_estimates %>% filter(removals == 0)
no_removals_R0 <- reproductive_number %>% filter(removals == 0)

# Figure 3
no_removals_alarm <- alarm %>% filter(removals == 0)
influenza_alarm <- filter(no_removals_alarm, disease == "influenza")
measles_alarm <- filter(no_removals_alarm, disease == "measles")

# Figure 4 & Table 1
selected_years <- c(1945, 1960, 1979)
sensitivity_posterior_fit <- predicted_fit %>% 
  filter(year %in% selected_years, disease == "measles")
sensitivity_results <- parameter_estimates %>% 
  filter(year %in% selected_years, disease == "measles")

# to map removal numbers from file names into the scenario considered
# sourced from sensitivityAnalysis.R
removals_legend <- c(
  "0"       = "No removals",
  # 1945
  "3028204" = "Age >15 removed",
  "3865839" = "95% removals",
  # 1960
  "4187137" = "Age >15 removed",
  "5912930" = "95% removals",
  # 1979
  "6608153" = "Age >15 removed",
  "8309572" = "95% removals"
)

################################################################################
# Figure 1 & 2: Posterior mean estimates and 95% credible intervals for k and R0

# Figure S1 & S2: Posterior mean estimates and 95% credible intervals for beta and gamma
################################################################################

plot_estimates <- function(data, parameter, quantity, y_format, y_lab = NULL) {
  df <- filter(data, param == parameter)
  ggplot(df, aes(x = year, y = mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.5) +
    scale_y_continuous(labels = y_format) +
    scale_x_continuous(breaks = pretty(df$year, n = 6)) +
    labs(title = quantity, x = "Year", y = y_lab) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      plot.title = element_text(face = "italic", hjust = 0.5, size = 22),
      axis.title.y = if (is.null(y_lab)) element_blank() else element_text(size = 14),
      axis.title.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      strip.background = element_blank(),
      strip.text = element_blank()
    )
}

# Combines k and R0 plots
plot_k_R0 <- function(disease_name, filename) {
  df_params <- filter(no_removals_parameters, disease == disease_name)
  df_r0 <- filter(no_removals_R0, disease == disease_name)
  
  p_k <- plot_estimates(df_params, "k",  expression(italic(k)),  
                        label_number(accuracy = 0.0001), y_lab = "Posterior Mean")
  p_R0 <- plot_estimates(df_r0, "R0", bquote(italic(.("\u211B"))[0]), 
                         label_number(accuracy = 0.1))
  ggsave(filename, p_k | p_R0, width = 10.5, height = 4.25, dpi = 600)
}

plot_k_R0("influenza", "influenza_k_R0.png") # Figure 1
plot_k_R0("measles", "measles_k_R0.png") # Figure 2

# Combines beta and gamma plots
plot_beta_gamma <- function(disease_name, filename) {
  df <- filter(no_removals_parameters, disease == disease_name)
  
  p_beta <- plot_estimates(df, "beta",  expression(italic(beta)),  
                           label_number(accuracy = 0.1), y_lab = "Posterior Mean")
  p_gamma <- plot_estimates(df, "rateI", expression(italic(gamma)), 
                            label_number(accuracy = 0.1))
  
  ggsave(filename, p_beta | p_gamma, width = 10.5, height = 4.25, dpi = 600)
}

plot_beta_gamma("influenza", "influenza_beta_gamma.png") # Figure S1
plot_beta_gamma("measles", "measles_beta_gamma.png") # Figure S2

################################################################################
# Figure 3: Maximum alarm value for each outbreak
################################################################################

plot_alarm <- function(data, title, x_breaks, x_limits, y_lab = NULL) {
  ggplot(data, aes(x = year, y = mean)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, linewidth = 0.5) +
    geom_point(size = 2, color = "black") +
    scale_x_continuous(breaks = x_breaks, limits = x_limits) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Year", y = y_lab, title = title) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.title.y = if (is.null(y_lab)) element_blank() else element_text(size = 14),
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 13),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      strip.background = element_blank(),
      strip.text = element_blank()
    )
}

p_influenza <- plot_alarm(influenza_alarm, "Influenza", seq(1950, 2020, 5), 
                          c(1950, 2020), y_lab = "Maximum Alarm")
p_measles <- plot_alarm(measles_alarm, "Measles", seq(1940, 1980, 5), c(1940, 1980))

ggsave("max_alarm.png", p_influenza + p_measles + plot_layout(ncol = 2), 
       width = 10.5, height = 4.25, dpi = 600)

################################################################################
# Table 1: Posterior estimates for k and R0 for selected years

# Table S1: Posterior estimates for beta and gamma for selected years
################################################################################

# Table 1: k and R0
k_values <- sensitivity_results %>%
  filter(param == "k") %>%
  mutate(simulation = recode(as.character(removals), !!!removals_legend),
         simulation = factor(simulation, levels = unique(removals_legend))) %>%
  arrange(year, simulation) %>%
  select(year, simulation, param, mean, lower, upper) %>%
  mutate(across(c(mean, lower, upper), ~ formatC(., format = "f", digits = 7)))

R0_values <- map_dfr(results, "R0") %>%
  filter(year %in% selected_years, disease == "measles") %>%
  mutate(simulation = recode(as.character(removals), !!!removals_legend),
         simulation = factor(simulation, levels = unique(removals_legend)),
         param = "R0") %>%
  arrange(year, simulation) %>%
  select(year, simulation, param, mean, lower, upper) %>%
  mutate(across(c(mean, lower, upper), ~ formatC(., format = "f", digits = 4)))

bind_rows(k_values, R0_values) %>%
  arrange(year, simulation) %>%
  as.data.frame() %>%
  print()

# Table S1: beta and gamma
sensitivity_results %>%
  filter(param %in% c("beta", "rateI")) %>%
  mutate(simulation = recode(as.character(removals), !!!removals_legend),
         simulation = factor(simulation, levels = unique(removals_legend))) %>%
  arrange(year, simulation, param) %>%
  select(year, simulation, param, mean, lower, upper) %>%
  mutate(across(c(mean, lower, upper), ~ as.character(round(., 4)))) %>%
  as.data.frame() %>%
  print()

################################################################################
# Figure 4: Posterior predictive fits for selected measles years
################################################################################

predictive_fit <- function(data, assumption) {
  ggplot(data, aes(x = date)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = adjustcolor("steelblue", alpha.f = 0.3)) +
    geom_line(aes(y = mean), color = "steelblue", linetype = "dashed", linewidth = 0.7) +
    geom_line(aes(y = obs),  color = "black", linewidth = 0.7) +
    scale_x_date(date_labels = "%b", date_breaks = "2 months") +
    labs(x = "Date", y = "Incidence", title = assumption) +
    theme_minimal(base_size = 15) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 15),
      plot.margin = margin(t = 5, r = 5, b = 20, l = 5)
    )
}

row_plots <- map(selected_years, function(yr) {
  year_data <- filter(all_fits, year == yr)
  sims <- levels(year_data$simulation)
  
  plots <- map(sims, ~ predictive_fit(filter(year_data, simulation == .x), .x))
  
  year_label <- wrap_elements(textGrob(as.character(yr),
                                       gp = gpar(fontsize = 20, fontface = "bold")))
  
  year_label / wrap_plots(plots, nrow = 1) + plot_layout(heights = c(0.08, 1))
})

final_plot <- wrap_plots(row_plots, ncol = 1) +
  plot_annotation(theme = theme(plot.margin = margin(20, 20, 20, 20)))

ggsave("fits_grid.png", final_plot, width = 11.5, height = 13.5, dpi = 600)






