# BC-CanadianDiseaseOutbreaks
This repository contains code associated with the analysis presented in "Behavioural Change in Canadian Disease Outbreaks" by MacLean & Deardon.

## Data
Incidence data for this study were obtained from the Canadian Disease Incidence Dataset (CANDID) by Earn et al. (2024). Population data required for Section 3.2 were obtained from Statistics Canada (2000, 2024). Filtered datasets for measles and influenza in Ontario, Canada, are provided in the `data/` folder.

## File Descriptions

| File | Description |
|------|-------------|
| `datasetPrep.R` | Filters CANDID data for measles and influenza in Ontario. |
| `dateSelection.R` | Produces `OutbreakDates_influenza.csv` and `OutbreakDates_measles.csv`. |
| `runModels.R` | Runs the full analysis across all years or individual select years. |
| `sensitivityAnalysis.R` | Computes initial removals under the three immunity assumptions (Section 3.2). |
| `getModelInputs.R` | Adapted from Ward et al. (2023) |
| `modelCodes.R` | Adapted from Ward et al. (2023) |
| `modelFits.R` | Adapted from Ward et al. (2023) |
| `postPredFit.R` | Adapted from Ward et al. (2023) |
| `summarizePost.R` | Adapted from Ward et al. (2023) |
| `predictiveFit_convergence.R` | Creates plots of posterior predictive fit, parameter traceplots, and exports Gelman-Rubin convergence diagnostics for all years. |
| `summarizeResults.R` | Summarizes estimates and posterior fit from the `results/` folder (not included in this repository) into smaller files to be used for figures and tables. |
| `figures.R` | Produces figures and output for tables in the manuscript. |

---

## Reproducing Results
Parameter estimates and summaries needed for figure reproduction are provided in `summarized_results/`. Raw output files for all years can be reproduced by running `runModels.R`, but are too large to include in this repository.
