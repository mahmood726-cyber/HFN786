# CNMA: Comprehensive Network Meta-Analysis

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![R](https://img.shields.io/badge/R-%E2%89%A5%204.0.0-blue.svg)](https://www.r-project.org/)

A comprehensive toolkit for conducting network meta-analysis with both frequentist and Bayesian approaches. CNMA provides extensive features for systematic reviews and meta-analyses in healthcare and other research fields.

## Features

### Core Functionality
- **Frequentist Network Meta-Analysis**: Based on the proven `netmeta` package
- **Bayesian Analysis**: Support for MCMC-based Bayesian NMA (planned)
- **Automatic Reference Selection**: Intelligently selects the most common treatment as reference
- **Parallel Processing**: Built-in support for multi-core processing with `future` framework

### Advanced Capabilities (Planned)
- **Transportability Analysis**: Weight studies based on similarity to target population
- **GRADE Quality Weighting**: Incorporate evidence quality into analysis
- **Risk of Bias (RoB2) Assessment**: Weight studies by bias risk
- **Meta-Regression**: Explore treatment effect modifiers with spline support
- **Extensive Diagnostics**:
  - Leave-one-out analysis
  - Leave-one-treatment-out analysis
  - Node-splitting for inconsistency detection
  - Publication bias assessment (PET-PEESE, trim-and-fill, selection models)
- **Interactive Visualizations**: Network plots, forest plots, and rankings

## Installation

### From GitHub (Development Version)

```r
# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install cnma
devtools::install_github("mahmood726-cyber/HFN786")
```

### Dependencies

CNMA requires the following packages:
- **Core**: `netmeta`, `ggplot2`, `dplyr`, `tidyr`, `purrr`, `tibble`
- **Optional**: `digest` (for caching), `future` and `future.apply` (for parallel processing)
- **Suggested**: `testthat`, `knitr`, `rmarkdown` (for development)

## Quick Start

### Basic Example

```r
library(cnma)

# Generate simulated data for demonstration
data <- simulate_cnma_data(n_studies = 30, seed = 123)

# Configure analysis (using basic settings)
config <- setup_cnma(
  sm = "HR",                    # Hazard ratio
  use_bayesian = FALSE,         # Use frequentist approach
  export_plots = FALSE,         # Don't export plots to files
  report_html = FALSE          # Don't generate HTML report
)

# Run analysis
results <- run_cnma_analysis(
  data = data,
  config = config
)

# View results
print(results)
summary(results)
```

### One-Line Quick Start

```r
# Quickstart with all defaults
results <- cnma_quickstart(n_studies = 20)
```

## Detailed Usage

### 1. Data Preparation

CNMA expects data in pairwise format with the following columns:
- `studlab`: Study identifier
- `treat1`: First treatment
- `treat2`: Second treatment
- `TE`: Treatment effect estimate (e.g., log hazard ratio)
- `seTE`: Standard error of treatment effect

```r
# Your data should look like this:
data <- data.frame(
  studlab = c("Study1", "Study1", "Study2"),
  treat1 = c("Placebo", "Placebo", "Placebo"),
  treat2 = c("DrugA", "DrugB", "DrugA"),
  TE = c(log(0.85), log(0.75), log(0.80)),
  seTE = c(0.15, 0.20, 0.18)
)

# Clean and validate
data_clean <- cnma_clean_data(data)
validate_cnma_input(data_clean)
```

### 2. Configuration

Create a configuration object to control analysis parameters:

```r
config <- setup_cnma(
  # Summary measure
  sm = "HR",                    # "HR", "OR", "RR", "MD", "SMD"

  # Analysis options
  use_bayesian = FALSE,         # Use Bayesian methods (requires JAGS)
  run_metareg = FALSE,          # Run meta-regression

  # Parallel processing
  parallel_strategy = "sequential",  # or "multisession", "multicore"
  n_cores = 2,

  # Reproducibility
  seed = 42,
  enable_cache = FALSE,

  # Output
  export_results = FALSE,
  export_plots = FALSE,
  report_html = FALSE
)

# View configuration
print(config)
```

### 3. Run Analysis

```r
# Run with automatic reference selection
results <- run_cnma_analysis(
  data = data,
  config = config
)

# Or specify reference treatment
results <- run_cnma_analysis(
  data = data,
  ref_treatment = "Placebo",
  config = config
)
```

### 4. Access Results

```r
# Main NMA results (netmeta object)
nma <- results$results$main_nma

# View network structure
if (requireNamespace("netmeta", quietly = TRUE)) {
  netmeta::netgraph(nma, plastic = FALSE)
}

# Treatment effects vs reference
print(nma$TE.random)

# Heterogeneity
cat("Tau:", nma$tau, "\n")
cat("I²:", nma$I2.random * 100, "%\n")
```

## Parallel Processing

CNMA supports parallel processing for computationally intensive operations:

```r
# Enable parallel processing
cnma_parallel_on(strategy = "multisession", workers = 4)

# Run analysis (will use parallel processing where applicable)
results <- run_cnma_analysis(data, config = config)

# Disable when done
cnma_parallel_off()
```

## Advanced Features

### Meta-Regression (Planned)

Explore treatment effect modifiers:

```r
config <- setup_cnma(
  run_metareg = TRUE,
  metareg_covariates = c("age_mean", "female_pct", "bmi_mean"),
  metareg_spline_covars = c("age_mean"),  # Use splines for age
  metareg_cr2 = TRUE                       # Use CR2 variance correction
)
```

### Transportability Analysis (Planned)

Weight studies based on similarity to target population:

```r
# Define target population
target_pop <- data.frame(
  age_mean = 70,
  female_pct = 0.52,
  bmi_mean = 27
)

config <- setup_cnma(
  use_transport = TRUE,
  transport_metric = "mahalanobis"
)

results <- run_cnma_analysis(
  data = data,
  target_population = target_pop,
  config = config
)
```

### Bayesian Analysis (Planned)

```r
config <- setup_cnma(
  use_bayesian = TRUE,
  bayes_chains = 3,
  bayes_iter = 10000,
  bayes_warmup = 5000,
  bayes_nodesplit = TRUE        # Run node-splitting for inconsistency
)
```

## Data Simulation

Generate realistic network meta-analysis data for testing:

```r
# Generate data with different network structures
data_small <- simulate_cnma_data(n_studies = 20, seed = 1)
data_large <- simulate_cnma_data(n_studies = 100, seed = 2)

# Includes covariates for meta-regression
head(data_small)
# Columns: studlab, treat1, treat2, TE, seTE, year, is_rct,
#          study_design, grade, rob2, age_mean, female_pct,
#          bmi_mean, charlson, baseline_risk
```

## Package Options

Control package behavior with global options:

```r
# Suppress messages
options(cnma.quiet = TRUE)

# Restore messages
options(cnma.quiet = FALSE)
```

## Validation and Data Quality

```r
# Validate input data
validate_cnma_input(data)  # Throws error if invalid

# Clean data automatically
data_clean <- cnma_clean_data(data)  # Removes invalid rows

# Check for:
# - Missing required columns
# - Non-finite values
# - Non-positive standard errors
# - Proper data types
```

## Performance Tips

1. **Caching**: Enable caching for repeated analyses
   ```r
   config <- setup_cnma(enable_cache = TRUE)
   ```

2. **Parallel Processing**: Use for large networks (>50 studies)
   ```r
   cnma_parallel_on("multisession", workers = parallel::detectCores() - 1)
   ```

3. **Selective Features**: Disable unused features to speed up analysis
   ```r
   config <- setup_cnma(
     use_bayesian = FALSE,
     run_metareg = FALSE,
     run_loo = FALSE
   )
   ```

## Current Limitations

**Note**: This is version 1.0.0. Many advanced features are planned but not yet implemented:

- ❌ Bayesian analysis (JAGS/Stan backends)
- ❌ Transportability weighting
- ❌ GRADE weighting
- ❌ Meta-regression
- ❌ Leave-one-out / Leave-one-treatment-out
- ❌ Publication bias tests (PET-PEESE, Copas, trim-and-fill)
- ❌ Interactive visualizations
- ❌ HTML report generation

Currently implemented:
- ✅ Core frequentist NMA (via netmeta)
- ✅ Data validation and cleaning
- ✅ Automatic reference selection
- ✅ Parallel processing infrastructure
- ✅ Data simulation
- ✅ Comprehensive configuration system

## Examples

### Complete Workflow

```r
library(cnma)

# 1. Load or simulate data
data <- simulate_cnma_data(n_studies = 40, seed = 42)

# 2. Explore data
table(data$treat1, data$treat2)
summary(data$TE)

# 3. Configure analysis
config <- setup_cnma(
  sm = "HR",
  use_bayesian = FALSE,
  seed = 42
)

# 4. Run analysis
results <- run_cnma_analysis(
  data = data,
  config = config
)

# 5. Examine results
print(results)
summary(results)

# 6. Access underlying netmeta object
nma <- results$results$main_nma
print(nma)

# 7. Create visualizations (requires netmeta)
if (requireNamespace("netmeta", quietly = TRUE)) {
  # Network graph
  netmeta::netgraph(nma, plastic = FALSE)

  # Forest plot
  netmeta::forest(nma, reference.group = results$ref_treatment)

  # Net heat plot (inconsistency)
  netmeta::netheat(nma)

  # Treatment ranking
  netmeta::netrank(nma)
}
```

## Getting Help

- **Documentation**: See function help with `?run_cnma_analysis`, `?setup_cnma`, etc.
- **Issues**: Report bugs at https://github.com/mahmood726-cyber/HFN786/issues
- **Questions**: Open a discussion on GitHub

## Citation

If you use CNMA in your research, please cite:

```
CNMA Development Team (2025). cnma: Comprehensive Network Meta-Analysis.
R package version 1.0.0. https://github.com/mahmood726-cyber/HFN786
```

## License

Apache License 2.0 - see [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please see CONTRIBUTING.md for guidelines.

## Acknowledgments

CNMA builds on the excellent work of:
- [`netmeta`](https://cran.r-project.org/package=netmeta) by Guido Schwarzer et al.
- [`dplyr`](https://dplyr.tidyverse.org/) and the tidyverse ecosystem
- [`future`](https://future.futureverse.org/) for parallel processing

## Development Status

🚧 **Active Development**: This package is under active development. The API may change. Please check the [CHANGELOG](CHANGELOG.md) for updates.

### Roadmap

**Version 1.1** (Q2 2025):
- Meta-regression implementation
- Publication bias diagnostics
- Enhanced visualization options

**Version 1.2** (Q3 2025):
- Bayesian analysis backends (JAGS, Stan)
- Node-splitting for inconsistency detection
- Interactive network plots

**Version 2.0** (Q4 2025):
- Transportability analysis
- GRADE/RoB2 weighting
- Machine learning features
- HTML report generation

---

**Built with ❤️ for evidence synthesis**
