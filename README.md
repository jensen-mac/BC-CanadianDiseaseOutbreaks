# BC-CanadianDiseaseOutbreaks
This repository contains code associated with the analysis presented in "Behavioural Change in Canadian Disease Outbreaks" by MacLean & Deardon.

## Data
Incidence data for this study were obtained from the Canadian Disease Incidence Dataset (CANDID) by Earn et al. (2024). Filtered datasets for measles and influenza in Ontario, Canada, are provided in the `Data/` folder.

## File Descriptions

| File | Description |
|------|-------------|
| `datasetPrep.R` | Filters CANDID data for measles and influenza in Ontario. |
| `dateSelection.R` | Produces `OutbreakDates_influenza.csv` and `OutbreakDates_measles.csv`. |
| `sensitivityAnalysis.R` | Computes initial removals under the three immunity assumptions (Section 3.2). |
| `runModels.R` | Runs the full analysis across all years or individual select years. |
| `getModelInputs.R` | Adapted from Ward et al. (2023) |
| `modelCodes.R` | Adapted from Ward et al. (2023) |
| `modelFits.R` | Adapted from Ward et al. (2023) |
| `postPredFit.R` | Adapted from Ward et al. (2023) |
| `summarizePost.R` | Adapted from Ward et al. (2023) |

---

## Reproducing Results
Parameter estimates and summaries needed for figure reproduction are provided in `summarized_results/`. Raw output files for all years can be reproduced by running `runModels.R`, but are too large to include in this repository.
