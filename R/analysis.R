# =========================================================
# Main Analysis Function
# =========================================================

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
    treat_freq <- table(c(data$treat1, data$treat2))
    if (length(treat_freq) == 0) {
      .stop_hint("No treatments found in data.",
                 "Check that treat1 and treat2 columns contain valid treatment names.")
    }
    ref_treatment <- names(sort(treat_freq, decreasing = TRUE))[1]
    msg("Auto-selected reference: '%s'", ref_treatment)
  }

  # Validate reference treatment exists in data
  all_treatments <- unique(c(data$treat1, data$treat2))
  if (!ref_treatment %in% all_treatments) {
    .stop_hint(sprintf("Reference treatment '%s' not found in data.", ref_treatment),
               sprintf("Available treatments: %s", paste(all_treatments, collapse = ", ")))
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
