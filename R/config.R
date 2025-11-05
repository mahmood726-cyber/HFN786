# =========================================================
# Configuration Functions
# =========================================================

#' Setup CNMA Configuration
#'
#' Creates a configuration object for CNMA analysis with all necessary parameters.
#'
#' @param sm Summary measure (default: "HR" for hazard ratio)
#' @param use_transport Use transportability weighting
#' @param transport_metric Distance metric for transportability
#' @param transport_kernel Kernel function for transportability
#' @param transport_truncation Truncation for transport weights
#' @param min_weight Minimum weight value
#' @param use_grade_weighting Use GRADE quality weighting
#' @param use_design_weighting Use study design weighting
#' @param design_var Variable for study design
#' @param design_weight_map Mapping of design values to weights
#' @param use_rob2_weighting Use RoB2 weighting
#' @param rob2_var Variable for risk of bias
#' @param use_bayesian Run Bayesian analysis
#' @param bayes_link Link function for Bayesian model
#' @param bayes_likelihood Likelihood for Bayesian model
#' @param bayes_priors Prior specifications
#' @param bayes_nodesplit Run Bayesian node-splitting
#' @param bayes_metareg_covariates Covariates for Bayesian meta-regression
#' @param bayes_chains Number of MCMC chains
#' @param bayes_iter MCMC iterations
#' @param bayes_warmup MCMC warmup iterations
#' @param one_stage_glmm Run one-stage GLMM
#' @param glmm_outcome Outcome type for GLMM
#' @param glmm_link Link function for GLMM
#' @param run_loo Run leave-one-out analysis
#' @param run_loto Run leave-one-treatment-out analysis
#' @param pet_peese_min_k Minimum studies for PET-PEESE
#' @param run_selection_models Run selection models
#' @param run_copas Run Copas analysis
#' @param run_trimfill Run trim-and-fill
#' @param run_netheat Generate net heat plot
#' @param nodesplit_min_k Minimum studies for node-splitting
#' @param run_metareg Run meta-regression
#' @param metareg_covariates Covariates for meta-regression
#' @param metareg_spline_covars Covariates for spline transformation
#' @param metareg_spline_df Degrees of freedom for splines
#' @param metareg_spline_cv Use cross-validation for spline df
#' @param metareg_cv_folds Number of CV folds
#' @param metareg_cr2 Use CR2 variance correction
#' @param baseline_risk_var Variable for baseline risk
#' @param class_map Treatment class mapping
#' @param dose_response_var Variable for dose-response
#' @param use_ml Use machine learning for heterogeneity
#' @param use_interactive Generate interactive plots
#' @param ggtheme ggplot2 theme
#' @param prediction_profiles Profiles for prediction
#' @param ref_rotation Perform reference rotation
#' @param sucra_boot Bootstrap SUCRA
#' @param sucra_boot_iter SUCRA bootstrap iterations
#' @param enable_ume Enable UME inconsistency model
#' @param enable_egger Enable network-level Egger test
#' @param export_results Export results to files
#' @param export_plots Export plots
#' @param plot_dir Directory for plots
#' @param output_dir Directory for output
#' @param report_html Generate HTML report
#' @param report_file HTML report filename
#' @param write_manifest Write run manifest
#' @param auto_install Auto-install missing packages
#' @param parallel_strategy Parallel processing strategy
#' @param n_cores Number of cores for parallel processing
#' @param seed Random seed
#' @param enable_cache Enable caching
#' @param log_to_file Log to file
#' @param log_file Log filename
#' @return Configuration list
#' @export
#' @examples
#' config <- setup_cnma(sm = "OR", use_bayesian = FALSE)
setup_cnma <- function(
  sm = "HR",
  use_transport = TRUE,
  transport_metric = c("mahalanobis", "euclidean"),
  transport_kernel = c("gaussian", "tricube"),
  transport_truncation = 0.02,
  min_weight = 1e-6,
  use_grade_weighting = TRUE,
  use_design_weighting = TRUE,
  design_var = "is_rct",
  design_weight_map = c(`1` = 1.0, `0` = 0.6),
  use_rob2_weighting = FALSE,
  rob2_var = "rob2",
  use_bayesian = TRUE,
  bayes_link = "identity",
  bayes_likelihood = "normal",
  bayes_priors = list(),
  bayes_nodesplit = TRUE,
  bayes_metareg_covariates = NULL,
  bayes_chains = 3,
  bayes_iter = 10000,
  bayes_warmup = 5000,
  one_stage_glmm = FALSE,
  glmm_outcome = c("binary", "continuous"),
  glmm_link = "logit",
  run_loo = TRUE,
  run_loto = TRUE,
  pet_peese_min_k = 10,
  run_selection_models = TRUE,
  run_copas = TRUE,
  run_trimfill = TRUE,
  run_netheat = TRUE,
  nodesplit_min_k = 10,
  run_metareg = TRUE,
  metareg_covariates = c("age_mean", "female_pct", "bmi_mean", "charlson"),
  metareg_spline_covars = c("age_mean", "bmi_mean"),
  metareg_spline_df = 3,
  metareg_spline_cv = TRUE,
  metareg_cv_folds = 5,
  metareg_cr2 = TRUE,
  baseline_risk_var = NULL,
  class_map = NULL,
  dose_response_var = NULL,
  use_ml = TRUE,
  use_interactive = FALSE,
  ggtheme = ggplot2::theme_minimal(),
  prediction_profiles = list(),
  ref_rotation = FALSE,
  sucra_boot = FALSE,
  sucra_boot_iter = 300,
  enable_ume = TRUE,
  enable_egger = TRUE,
  export_results = FALSE,
  export_plots = TRUE,
  plot_dir = "cnma_plots",
  output_dir = "cnma_results",
  report_html = TRUE,
  report_file = "cnma_report.html",
  write_manifest = TRUE,
  auto_install = FALSE,
  parallel_strategy = "sequential",
  n_cores = max(1, min(4, parallel::detectCores() - 1)),
  seed = 42,
  enable_cache = FALSE,
  log_to_file = FALSE,
  log_file = "cnma_log.txt"
) {
  # Validate inputs
  stopifnot(
    "sm must be character" = is.character(sm) && length(sm) == 1,
    "bayes_chains must be positive" = bayes_chains > 0,
    "bayes_iter must be positive" = bayes_iter > 0,
    "bayes_warmup must be positive" = bayes_warmup > 0,
    "bayes_warmup must be less than bayes_iter" = bayes_warmup < bayes_iter,
    "n_cores must be positive" = n_cores > 0,
    "seed must be numeric" = is.numeric(seed)
  )

  transport_metric <- match.arg(transport_metric)
  transport_kernel <- match.arg(transport_kernel)
  glmm_outcome <- match.arg(glmm_outcome)

  config <- as.list(environment())
  class(config) <- "cnma_config"
  config
}

#' Print CNMA configuration
#' @param x Configuration object
#' @param ... Additional arguments
#' @export
print.cnma_config <- function(x, ...) {
  cat("CNMA Configuration:\n")
  cat("  Summary measure:", x$sm, "\n")
  cat("  Bayesian analysis:", x$use_bayesian, "\n")
  cat("  Transportability:", x$use_transport, "\n")
  cat("  Meta-regression:", x$run_metareg, "\n")
  cat("  Parallel strategy:", x$parallel_strategy, "\n")
  invisible(x)
}
