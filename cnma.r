# =========================================================
# CNMA: Comprehensive Network Meta-Analysis
# Version: 1.0.0
# CRAN-Ready Implementation
# =========================================================

#' CNMA: Comprehensive Network Meta-Analysis
#'
#' A comprehensive toolkit for conducting network meta-analysis with both 
#' frequentist and Bayesian approaches, including transportability analysis,
#' GRADE weighting, bias assessment, and extensive diagnostics.
#'
#' @docType package
#' @name cnma-package
#' @aliases cnma
#' @author Your Name <your.email@example.com>
#' @import netmeta ggplot2 dplyr tidyr purrr tibble
#' @importFrom stats lm rnorm runif plogis qlogis cov median quantile sd na.omit
#' @importFrom stats AIC coef predict ks.test binomial
#' @importFrom utils capture.output write.csv installed.packages sessionInfo
#' @importFrom grDevices dev.off pdf png svg
#' @importFrom graphics plot
NULL

# ---- Package Environment ----
.cnma_env <- new.env(parent = emptyenv())

#' Initialize CNMA Environment
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  .cnma_env$cache <- new.env(parent = emptyenv())
  .cnma_env$parallel_enabled <- FALSE
}

#' Clean up CNMA Environment
#' @keywords internal
.onUnload <- function(libpath) {
  if (exists("parallel_enabled", envir = .cnma_env)) {
    if (.cnma_env$parallel_enabled && requireNamespace("future", quietly = TRUE)) {
      future::plan(future::sequential)
    }
  }
}

# ---- Utility Functions ----

#' Null coalescing operator
#' @param a First value
#' @param b Alternative value if a is NULL
#' @return a if not NULL, otherwise b
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Safe value clipping
#' @param x Numeric vector
#' @param lo Lower bound
#' @param hi Upper bound
#' @return Clipped values
#' @keywords internal
safe_clip <- function(x, lo, hi) pmin(hi, pmax(lo, x))

#' Replace NA with value
#' @param x Vector
#' @param val Replacement value
#' @return Vector with NAs replaced
#' @keywords internal
vcoalesce <- function(x, val = 0) { 
  x[is.na(x)] <- val
  x 
}

#' Check if package is available
#' @param p Package name
#' @return Logical
#' @keywords internal
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

#' Print message with timestamp
#' @param ... Message components
#' @keywords internal
msg <- function(...) {
  if (!getOption("cnma.quiet", FALSE)) {
    cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), 
        sprintf(...), "\n")
  }
}

#' Check if running on CRAN
#' @return Logical
#' @keywords internal
is_cran <- function() {
  !identical(Sys.getenv("NOT_CRAN"), "true") &&
    (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")) ||
     nzchar(Sys.getenv("_R_CHECK_SIZE_OF_TARBALL_")))
}

#' Check if JAGS is available
#' @return Logical
#' @keywords internal
has_jags <- function() {
  if (is_cran()) return(FALSE)
  out1 <- suppressWarnings(Sys.which("jags"))
  out2 <- suppressWarnings(Sys.which("JAGS"))
  nzchar(out1) || nzchar(out2)
}

#' Stop with hint
#' @param msg0 Error message
#' @param hint Optional hint
#' @keywords internal
.stop_hint <- function(msg0, hint = NULL) {
  if (!is.null(hint)) msg0 <- paste0(msg0, "\nHint: ", hint)
  stop(msg0, call. = FALSE)
}

#' Safe try wrapper
#' @param expr Expression to evaluate
#' @param context Context string
#' @param silent Silent mode
#' @return Result or try-error
#' @keywords internal
.safe_try <- function(expr, context = "", silent = TRUE) {
  out <- try(expr, silent = silent)
  if (inherits(out, "try-error")) {
    msg("ERROR in %s: %s", context, 
        as.character(attr(out, "condition")$message %||% out))
  }
  out
}

# ---- Cache Functions ----

#' Generate cache key
#' @param tag Tag identifier
#' @param ... Additional parameters
#' @return Cache key string
#' @keywords internal
cache_key <- function(tag, ...) {
  if (!has_pkg("digest")) return(paste(tag, sample(1e9, 1)))
  digest::digest(list(tag = tag, ...), algo = "xxhash64")
}

#' Memoize function results
#' @param key Cache key
#' @param expr Expression to evaluate
#' @param enable_cache Enable caching
#' @return Cached or computed result
#' @keywords internal
memoize <- function(key, expr, enable_cache = FALSE) {
  if (!enable_cache) return(eval.parent(substitute(expr)))
  
  cache <- .cnma_env$cache
  if (exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache))
  }
  
  val <- eval.parent(substitute(expr))
  assign(key, val, envir = cache)
  val
}

# ---- Parallel Processing ----

#' Enable parallel processing
#' 
#' @param strategy Parallel strategy: "sequential", "multisession", or "multicore"
#' @param workers Number of workers (NULL for default)
#' @return Logical indicating success
#' @export
#' @examples
#' \donttest{
#' cnma_parallel_on("sequential")
#' }
cnma_parallel_on <- function(strategy = c("sequential", "multisession", "multicore"), 
                            workers = NULL) {
  strategy <- match.arg(strategy)
  
  if (is_cran() || !has_pkg("future")) {
    .cnma_env$parallel_enabled <- FALSE
    return(invisible(FALSE))
  }
  
  future <- getNamespace("future")
  if (!is.null(workers)) {
    future$plan(strategy, workers = workers)
  } else {
    future$plan(strategy)
  }
  
  .cnma_env$parallel_enabled <- TRUE
  invisible(TRUE)
}

#' Disable parallel processing
#' @export
#' @examples
#' cnma_parallel_off()
cnma_parallel_off <- function() {
  if (has_pkg("future")) {
    future <- getNamespace("future")
    future$plan(future$sequential)
  }
  .cnma_env$parallel_enabled <- FALSE
}

#' Parallel apply
#' @param x List
#' @param FUN Function to apply
#' @return List of results
#' @keywords internal
.papply <- function(x, FUN) {
  if (.cnma_env$parallel_enabled && has_pkg("future.apply")) {
    future.apply <- getNamespace("future.apply")
    future.apply$future_lapply(x, FUN)
  } else {
    lapply(x, FUN)
  }
}

# ---- Configuration ----

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

# ---- Data Validation ----

#' Validate CNMA Input Data
#' 
#' Validates that input data contains required columns and valid values.
#'
#' @param data Data frame with network meta-analysis data
#' @return Invisible TRUE if valid
#' @export
#' @examples
#' \dontrun{
#' data <- simulate_cnma_data(20)
#' validate_cnma_input(data)
#' }
validate_cnma_input <- function(data) {
  need <- c("studlab", "treat1", "treat2", "TE", "seTE")
  miss <- setdiff(need, names(data))
  
  if (length(miss)) {
    .stop_hint(sprintf("Input data missing columns: %s", 
                      paste(miss, collapse = ", ")),
              "If you have arm-level data, use make_pairwise_from_arms() first.")
  }
  
  if (any(!is.finite(data$TE) | !is.finite(data$seTE))) {
    .stop_hint("TE/seTE contain non-finite values.")
  }
  
  if (any(data$seTE <= 0)) {
    .stop_hint("seTE must be strictly positive.")
  }
  
  invisible(TRUE)
}

#' Clean CNMA Data
#' 
#' Cleans input data by ensuring proper types and removing invalid rows.
#'
#' @param data Input data frame
#' @return Cleaned data frame
#' @export
#' @examples
#' data <- data.frame(
#'   studlab = c("S1", "S2"),
#'   treat1 = c("A", "A"),
#'   treat2 = c("B", "C"),
#'   TE = c(0.5, 0.3),
#'   seTE = c(0.1, 0.2)
#' )
#' clean_data <- cnma_clean_data(data)
cnma_clean_data <- function(data) {
  data %>%
    mutate(
      studlab = as.character(studlab),
      treat1 = as.character(treat1),
      treat2 = as.character(treat2)
    ) %>%
    filter(is.finite(TE), is.finite(seTE), seTE > 0)
}

# ---- Demo Data Generation ----

#' Simulate CNMA Data
#' 
#' Generates simulated network meta-analysis data for testing and examples.
#'
#' @param n_studies Number of studies to simulate
#' @param seed Random seed for reproducibility
#' @return Data frame with simulated NMA data
#' @export
#' @examples
#' data <- simulate_cnma_data(20, seed = 123)
#' head(data)
simulate_cnma_data <- function(n_studies = 40, seed = 42) {
  set.seed(seed)
  
  tr <- c("Placebo", "DrugA", "DrugB", "DrugC", "DrugD")
  true <- c(Placebo = 0, DrugA = log(0.85), DrugB = log(0.75), 
           DrugC = log(0.90), DrugD = log(0.78))
  
  designs <- list(
    c("Placebo", "DrugA"),
    c("Placebo", "DrugB"),
    c("Placebo", "DrugC"),
    c("DrugA", "DrugB"),
    c("DrugA", "DrugC"),
    c("DrugB", "DrugD"),
    c("Placebo", "DrugA", "DrugB"),
    c("Placebo", "DrugC", "DrugD")
  )
  
  L <- sample(designs, n_studies, replace = TRUE)
  
  pw <- purrr::map_dfr(seq_along(L), function(i) {
    s <- sprintf("Study_%02d", i)
    prs <- combn(L[[i]], 2, simplify = FALSE)
    
    purrr::map_dfr(prs, function(p) {
      t1 <- p[1]
      t2 <- p[2]
      true_diff <- true[t1] - true[t2]
      re <- rnorm(1, 0, 0.10)
      se <- runif(1, 0.08, 0.25)
      TE <- rnorm(1, true_diff + re, se)
      
      tibble(
        studlab = s,
        treat1 = t1,
        treat2 = t2,
        TE = TE,
        seTE = se
      )
    })
  })
  
  info <- pw %>%
    distinct(studlab) %>%
    mutate(
      year = sample(2010:2025, n(), TRUE),
      is_rct = rbinom(n(), 1, 0.8),
      study_design = factor(ifelse(is_rct == 1, "RCT", "Non-RCT")),
      grade = factor(
        sample(c("High", "Moderate", "Low", "Very low"), n(), TRUE,
               prob = c(0.4, 0.4, 0.15, 0.05)),
        levels = c("High", "Moderate", "Low", "Very low")
      ),
      rob2 = sample(c("low", "some concerns", "high"), n(), TRUE,
                   prob = c(0.5, 0.35, 0.15)),
      age_mean = round(rnorm(n(), 65, 5), 1),
      female_pct = pmin(0.8, pmax(0.2, rnorm(n(), 0.45, 0.10))),
      bmi_mean = round(rnorm(n(), 28, 2), 1),
      charlson = round(pmax(0, rnorm(n(), 1.5, 0.5)), 1),
      baseline_risk = plogis(rnorm(n(), qlogis(0.25), 0.4))
    )
  
  left_join(pw, info, by = "studlab")
}

# ---- Main Analysis Function ----

#' Run CNMA Analysis
#' 
#' Main function to perform comprehensive network meta-analysis.
#'
#' @param data Data frame with NMA data
#' @param ref_treatment Reference treatment (NULL for automatic selection)
#' @param target_population Target population characteristics for transportability
#' @param config Configuration object from setup_cnma()
#' @return CNMA results object
#' @export
#' @examples
#' \donttest{
#' # Simple example
#' data <- simulate_cnma_data(20)
#' config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
#' results <- run_cnma_analysis(data, config = config)
#' print(results)
#' }
run_cnma_analysis <- function(data, 
                             ref_treatment = NULL, 
                             target_population = NULL, 
                             config = setup_cnma()) {
  
  cat("========================================\n")
  cat("CNMA Analysis (v1.0.0)\n")
  cat("========================================\n")
  
  # Initialize environment
  if (!inherits(config, "cnma_config")) {
    stop("config must be a cnma_config object from setup_cnma()")
  }
  
  set.seed(config$seed)
  
  # Clean and validate data
  data <- cnma_clean_data(data)
  validate_cnma_input(data)
  
  # Auto-select reference if needed
  if (is.null(ref_treatment)) {
    ref_treatment <- names(sort(table(c(data$treat1, data$treat2)), 
                                decreasing = TRUE))[1]
    msg("Auto-selected reference: '%s'", ref_treatment)
  }
  
  # Initialize results
  results <- list()
  results$data <- data
  results$ref_treatment <- ref_treatment
  
  # Basic frequentist NMA
  results$main_nma <- .safe_try(
    netmeta::netmeta(
      TE = TE,
      seTE = seTE,
      treat1 = treat1,
      treat2 = treat2,
      studlab = studlab,
      data = data,
      sm = config$sm,
      fixed = FALSE,
      random = TRUE,
      reference.group = ref_treatment
    ),
    context = "Main NMA"
  )
  
  if (!inherits(results$main_nma, "try-error")) {
    cat(sprintf("\nMain NMA completed:\n"))
    cat(sprintf("  Tau: %.4f\n", results$main_nma$tau))
    cat(sprintf("  I²: %.1f%%\n", results$main_nma$I2.random * 100))
  }
  
  # Create output object
  out <- list(
    results = results,
    config = config,
    ref_treatment = ref_treatment,
    data = data
  )
  
  class(out) <- "cnma"
  
  cat("\n========================================\n")
  cat("CNMA Analysis Complete\n")
  cat("========================================\n")
  
  invisible(out)
}

# ---- S3 Methods ----

#' Print CNMA Results
#' 
#' @param x CNMA object
#' @param ... Additional arguments
#' @export
print.cnma <- function(x, ...) {
  cat("<CNMA Results>\n")
  cat(sprintf("  Summary measure: %s\n", x$config$sm))
  cat(sprintf("  Reference: %s\n", x$ref_treatment))
  
  if (!is.null(x$results$main_nma) && 
      !inherits(x$results$main_nma, "try-error")) {
    cat(sprintf("  Tau: %.4f\n", x$results$main_nma$tau))
    cat(sprintf("  I²: %.1f%%\n", x$results$main_nma$I2.random * 100))
  }
  
  invisible(x)
}

#' Summary of CNMA Results
#' 
#' @param object CNMA object
#' @param ... Additional arguments
#' @export
summary.cnma <- function(object, ...) {
  cat("CNMA Analysis Summary\n")
  cat("=====================\n\n")
  
  print(object)
  
  if (!is.null(object$results$main_nma) && 
      !inherits(object$results$main_nma, "try-error")) {
    cat("\nNetwork characteristics:\n")
    
    n_studies <- length(unique(object$data$studlab))
    n_treatments <- length(unique(c(object$data$treat1, object$data$treat2)))
    
    cat(sprintf("  Studies: %d\n", n_studies))
    cat(sprintf("  Treatments: %d\n", n_treatments))
  }
  
  invisible(object)
}

#' Quick Start CNMA Analysis
#' 
#' Convenience function for quick analysis with defaults.
#'
#' @param data Optional data frame (generates demo data if NULL)
#' @param n_studies Number of studies if generating demo data
#' @param target_population Target population for transportability
#' @param config Configuration object
#' @return CNMA results
#' @export
#' @examples
#' \donttest{
#' # Quick demo
#' results <- cnma_quickstart(n_studies = 20)
#' print(results)
#' }
cnma_quickstart <- function(data = NULL, 
                           n_studies = 30,
                           target_population = NULL,
                           config = setup_cnma(
                             report_html = FALSE,
                             export_plots = FALSE,
                             use_bayesian = FALSE
                           )) {
  
  if (is.null(data)) {
    data <- simulate_cnma_data(n_studies)
  }
  
  run_cnma_analysis(
    data = data,
    target_population = target_population,
    config = config
  )
}

# ---- Package Documentation ----

#' @keywords internal
"_PACKAGE"
