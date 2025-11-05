# cnma 1.1.0

## Major Journal-Quality Enhancements

This release transforms cnma into a **publication-ready** toolkit aligned with best practices from leading statistics journals and PRISMA-NMA guidelines.

### New Features Based on Statistics Literature

#### Treatment Ranking & Comparison
* **`calculate_rankings()`** - P-scores and SUCRA for treatment ranking (Rücker & Schwarzer, 2015)
* **`create_league_table()`** - Standard BMJ/Lancet format pairwise comparisons
* **`calculate_prediction_intervals()`** - Prediction intervals accounting for heterogeneity (IntHout et al., 2016)

#### Inconsistency Assessment
* **`assess_inconsistency()`** - Comprehensive inconsistency evaluation with three methods:
  - Global: Cochran's Q, I², Tau²
  - Local: Node-splitting analysis (Dias et al., 2010)
  - Design: Design-by-treatment interaction (Krahn et al., 2013)
* **`calculate_contribution_matrix()`** - Study contributions to network estimates
* **`assess_transitivity()`** - Evaluate similarity of studies across comparisons

#### Publication-Quality Visualizations
* **`plot_network()`** - Network plot with study contributions
* **`plot_forest()`** - Forest plot with prediction intervals
* **`plot_rankings()`** - Rankograms and cumulative ranking plots
* **`plot_funnel()`** - Comparison-adjusted funnel plot (Chaimani & Salanti, 2012)
* **`plot_netheat()`** - Net heat plot for inconsistency visualization
* **`create_publication_plots()`** - Batch generate all plots in PNG/PDF

#### PRISMA-NMA Compliance
* **`generate_prisma_report()`** - Automated PRISMA-NMA compliance reporting (Hutton et al., 2015)
* **`network_characteristics_summary()`** - Comprehensive network summaries
* Support for all PRISMA-NMA extension items (S1-S5)

#### Sensitivity Analysis Framework
* **`leave_one_out_analysis()`** - Identify influential studies
* **`assess_publication_bias()`** - Publication bias assessment with multiple methods
* **`compare_models()`** - Fixed vs random effects model comparison

### Documentation Improvements
* Added comprehensive vignette: "Publication-Ready Network Meta-Analysis with CNMA"
* Step-by-step journal-quality workflow examples
* References to key methodological papers
* PRISMA-NMA checklist and reporting guidance

### Enhanced S3 Methods
* `print.cnma_ranking()` - Format treatment rankings
* `print.cnma_league_table()` - Format league tables
* `print.cnma_prediction_intervals()` - Format prediction intervals
* `print.cnma_inconsistency()` - Format inconsistency results
* `print.cnma_transitivity()` - Format transitivity assessment
* `print.cnma_loo()` - Format leave-one-out results
* `print.cnma_pub_bias()` - Format publication bias assessment
* `print.cnma_model_comparison()` - Format model comparison
* `print.cnma_network_summary()` - Format network characteristics

### Implementation of Best Practices
* Methods from *Statistics in Medicine* (Rücker & Schwarzer, 2015; Dias et al., 2010)
* Methods from *BMC Medical Research Methodology* (Krahn et al., 2013)
* Methods from *Annals of Internal Medicine* (Hutton et al., 2015)
* Methods from *Research Synthesis Methods* (Chaimani & Salanti, 2012)
* Methods from *BMJ Open* (IntHout et al., 2016)

### Package Metadata
* Updated DESCRIPTION to version 1.1.0
* Added PRISMA-NMA compliance badge
* Enhanced package description with methodology references
* Updated NAMESPACE with 30+ new exported functions

## Bug Fixes
* None (feature release)

## Internal Changes
* Added 5 new R source files (ranking.R, inconsistency.R, visualization.R, prisma.R, sensitivity.R)
* Enhanced modular code organization
* Added comprehensive roxygen2 documentation for all new functions

---

# cnma 1.0.0

## Major Changes

* Initial CRAN submission
* Core frequentist network meta-analysis functionality via `netmeta`
* Comprehensive configuration system with `setup_cnma()`
* Data validation and cleaning utilities
* Automatic reference treatment selection
* Parallel processing infrastructure
* Extensive documentation and examples

## Features

### Core Functionality
* `run_cnma_analysis()` - Main analysis function
* `cnma_quickstart()` - Quick start convenience function
* `setup_cnma()` - Comprehensive configuration
* `validate_cnma_input()` - Input validation
* `cnma_clean_data()` - Data cleaning
* `simulate_cnma_data()` - Demo data generation

### Infrastructure
* Modular code organization in R/ directory
* Comprehensive test suite with testthat (90+ tests)
* S3 methods for print and summary
* Parallel processing support with future framework
* Caching mechanism for repeated analyses

### Documentation
* Comprehensive README with examples
* Roxygen2 documentation for all exported functions
* Test coverage for core functionality
* CONTRIBUTING guidelines
* NEWS tracking

## Improvements

### Code Quality
* Split monolithic 701-line file into 8 logical modules
* Added comprehensive input validation
* Improved error messages with hints
* Better cache key generation (timestamp + random characters)
* Added validation for configuration parameters

### Testing
* 90+ unit tests covering all major functionality
* Test utilities, configuration, validation, analysis, and methods
* Tests for data simulation and parallel processing
* Comprehensive edge case coverage

### Documentation
* Professional README with badges and examples
* Clear usage examples and workflows
* Performance tips and best practices
* Known limitations clearly documented
* Development roadmap outlined

## Infrastructure

### New Files
* `DESCRIPTION` - Package metadata
* `NAMESPACE` - Export declarations
* `.gitignore` - Git ignore rules
* `.Rbuildignore` - Build ignore rules
* `NEWS.md` - This changelog
* `CONTRIBUTING.md` - Contribution guidelines

### Modular Structure
* `R/aaa-package.R` - Package documentation and lifecycle
* `R/utilities.R` - Utility functions
* `R/parallel.R` - Parallel processing
* `R/config.R` - Configuration system
* `R/validation.R` - Data validation
* `R/data-simulation.R` - Data generation
* `R/analysis.R` - Main analysis
* `R/methods.R` - S3 methods

## Known Limitations

The following features are planned but not yet implemented in v1.0.0:

* Bayesian analysis (JAGS/Stan backends)
* Transportability weighting
* GRADE quality weighting
* RoB2 weighting
* Meta-regression
* Leave-one-out / Leave-one-treatment-out analysis
* Publication bias diagnostics (PET-PEESE, Copas, trim-and-fill)
* Node-splitting for inconsistency detection
* Interactive visualizations
* HTML report generation

These features are on the roadmap for future versions.

## Breaking Changes

* None (initial release)

## Bug Fixes

* None (initial release)

## Internal Changes

* Improved cache key generation with timestamp and random characters
* Better error handling with `.stop_hint()` helper
* Enhanced validation in `simulate_cnma_data()`
* Added duplicate comparison detection in validation
* Improved messages for removed rows in data cleaning

## Dependencies

### Imports
* netmeta (>= 2.0.0)
* ggplot2 (>= 3.3.0)
* dplyr (>= 1.0.0)
* tidyr (>= 1.0.0)
* purrr (>= 0.3.0)
* tibble (>= 3.0.0)
* stats, utils, grDevices, graphics (base R)

### Suggests
* digest (>= 0.6.0) - for caching
* future (>= 1.0.0) - for parallel processing
* future.apply (>= 1.0.0) - for parallel processing
* testthat (>= 3.0.0) - for testing
* knitr, rmarkdown - for vignettes
* covr - for code coverage

## Acknowledgments

CNMA builds on the excellent work of:
* `netmeta` package by Guido Schwarzer et al.
* The tidyverse ecosystem
* The `future` framework for parallel processing

---

# Future Versions (Planned)

## Version 1.1.0 (Planned Q2 2025)
* Meta-regression implementation
* Publication bias diagnostics
* Enhanced visualization options
* Leave-one-out analysis

## Version 1.2.0 (Planned Q3 2025)
* Bayesian analysis backends (JAGS, Stan)
* Node-splitting for inconsistency detection
* Interactive network plots
* Treatment class effects

## Version 2.0.0 (Planned Q4 2025)
* Transportability analysis
* GRADE/RoB2 weighting
* Machine learning features
* HTML report generation
* One-stage models
