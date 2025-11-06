# cnma: Comprehensive Network Meta-Analysis

<!-- badges: start -->
[![R-CMD-check](https://github.com/mahmood726-cyber/HFN786/workflows/R-CMD-check/badge.svg)](https://github.com/mahmood726-cyber/HFN786/actions)
[![test-coverage](https://github.com/mahmood726-cyber/HFN786/workflows/test-coverage/badge.svg)](https://github.com/mahmood726-cyber/HFN786/actions)
<!-- badges: end -->

## Overview

**cnma** is a comprehensive toolkit for conducting network meta-analysis with both frequentist and Bayesian approaches. It provides a complete workflow for:

- Network meta-analysis with frequentist and Bayesian methods
- Transportability analysis and weighting
- GRADE evidence quality assessment
- Risk of bias (RoB2) integration
- Extensive diagnostics and heterogeneity assessment
- Advanced visualization and reporting

## Installation

You can install the development version from GitHub:

```r
# Install remotes if needed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install cnma
remotes::install_github("mahmood726-cyber/HFN786")
```

## Quick Start

```r
library(cnma)

# Generate example data
data <- simulate_cnma_data(n_studies = 30)

# Configure analysis
config <- setup_cnma(
  sm = "HR",              # Summary measure: Hazard Ratio
  use_bayesian = FALSE,   # Use frequentist approach
  export_plots = TRUE     # Save plots
)

# Run analysis
results <- run_cnma_analysis(data, config = config)

# View results
print(results)
summary(results)
```

## Features

### Core Analysis
- **Frequentist NMA**: Based on netmeta package with random-effects models
- **Bayesian NMA**: MCMC-based Bayesian network meta-analysis
- **Transportability**: Weight studies by similarity to target population
- **GRADE Weighting**: Incorporate evidence quality into analysis

### Diagnostics
- Heterogeneity assessment (τ², I²)
- Inconsistency detection (node-splitting, UME models)
- Publication bias tests (Egger's test, funnel plots)
- Leave-one-out analysis
- Net heat plots

### Meta-Regression
- Covariate adjustment
- Spline transformations with cross-validation
- Cluster-robust variance estimation
- Dose-response modeling

### Visualization
- Network plots
- Forest plots
- League tables
- SUCRA plots
- Contribution plots

## Example Workflow

```r
# Step 1: Prepare your data
data <- simulate_cnma_data(n_studies = 40)

# Step 2: Set up configuration
config <- setup_cnma(
  sm = "OR",
  use_transport = TRUE,
  use_grade_weighting = TRUE,
  use_bayesian = FALSE,
  run_metareg = TRUE,
  export_plots = TRUE,
  plot_dir = "my_plots"
)

# Step 3: Run analysis
results <- run_cnma_analysis(
  data = data,
  ref_treatment = "Placebo",
  config = config
)

# Step 4: Examine results
summary(results)
```

## Data Format

Your input data should have the following columns:

- `studlab`: Study identifier
- `treat1`: First treatment in comparison
- `treat2`: Second treatment in comparison
- `TE`: Treatment effect estimate (log scale for ratios)
- `seTE`: Standard error of treatment effect

Optional columns for advanced features:
- `grade`: GRADE quality rating (High/Moderate/Low/Very low)
- `is_rct`: Study design indicator (1 = RCT, 0 = observational)
- `rob2`: Risk of bias assessment (low/some concerns/high)
- Covariate columns: `age_mean`, `female_pct`, `bmi_mean`, etc.

## Documentation

For detailed documentation, see:

```r
?run_cnma_analysis  # Main analysis function
?setup_cnma         # Configuration options
?simulate_cnma_data # Data generation
```

## Citation

If you use this package in your research, please cite:

```
[Citation information to be added]
```

## License

MIT License - see LICENSE file for details

## Issues and Contributions

Report issues at: https://github.com/mahmood726-cyber/HFN786/issues

Contributions are welcome via pull requests.
