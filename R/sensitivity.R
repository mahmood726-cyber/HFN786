# =========================================================
# Sensitivity Analysis Functions
# Based on Cochrane Handbook and journal best practices
# =========================================================

#' Perform Leave-One-Out Sensitivity Analysis
#'
#' Evaluates the influence of each study by excluding it and re-running
#' the analysis. Recommended by Cochrane Handbook and BMJ guidelines.
#'
#' @param data Original data frame
#' @param config Configuration object
#' @param progress Show progress messages (default: TRUE)
#' @return Data frame with leave-one-out results
#' @export
#' @references
#' Higgins JPT, Thomas J, Chandler J, et al. (2023). Cochrane Handbook for
#' Systematic Reviews of Interventions version 6.4. Cochrane.
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(30)
#' config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
#' loo_results <- leave_one_out_analysis(data, config)
#' print(loo_results)
#' }
leave_one_out_analysis <- function(data, config = setup_cnma(), progress = TRUE) {
  if (!is.data.frame(data)) {
    .stop_hint("data must be a data.frame")
  }

  # Get unique studies
  studies <- unique(data$studlab)
  n_studies <- length(studies)

  if (progress) {
    msg("Running leave-one-out analysis for %d studies...", n_studies)
  }

  # Store results
  loo_results <- list()

  # Full model first
  full_model <- .safe_try(
    run_cnma_analysis(data, config = config),
    context = "Full model",
    silent = TRUE
  )

  if (inherits(full_model, "try-error") ||
      inherits(full_model$results$main_nma, "try-error")) {
    .stop_hint("Full model failed. Cannot proceed with sensitivity analysis.")
  }

  full_nma <- full_model$results$main_nma
  full_tau <- full_nma$tau
  full_I2 <- full_nma$I2.random

  # Leave-one-out loop
  for (i in seq_along(studies)) {
    study_name <- studies[i]

    if (progress && i %% 5 == 0) {
      msg("  Progress: %d/%d studies", i, n_studies)
    }

    # Exclude one study
    data_loo <- data[data$studlab != study_name, ]

    # Run analysis
    loo_model <- .safe_try(
      run_cnma_analysis(data_loo, config = config),
      context = sprintf("LOO: %s", study_name),
      silent = TRUE
    )

    if (!inherits(loo_model, "try-error") &&
        !inherits(loo_model$results$main_nma, "try-error")) {
      loo_nma <- loo_model$results$main_nma

      loo_results[[study_name]] <- data.frame(
        excluded_study = study_name,
        tau = loo_nma$tau,
        tau_change = loo_nma$tau - full_tau,
        tau_pct_change = ((loo_nma$tau - full_tau) / full_tau) * 100,
        I2 = loo_nma$I2.random,
        I2_change = loo_nma$I2.random - full_I2,
        influential = abs((loo_nma$tau - full_tau) / full_tau) > 0.10,
        stringsAsFactors = FALSE
      )
    } else {
      loo_results[[study_name]] <- data.frame(
        excluded_study = study_name,
        tau = NA,
        tau_change = NA,
        tau_pct_change = NA,
        I2 = NA,
        I2_change = NA,
        influential = NA,
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine results
  loo_df <- do.call(rbind, loo_results)
  rownames(loo_df) <- NULL

  # Add full model reference
  loo_df <- rbind(
    data.frame(
      excluded_study = "Full model (reference)",
      tau = full_tau,
      tau_change = 0,
      tau_pct_change = 0,
      I2 = full_I2,
      I2_change = 0,
      influential = FALSE,
      stringsAsFactors = FALSE
    ),
    loo_df
  )

  if (progress) {
    msg("Leave-one-out analysis completed!")
    influential_studies <- sum(loo_df$influential, na.rm = TRUE)
    if (influential_studies > 0) {
      msg("  Found %d influential study(ies)", influential_studies)
    } else {
      msg("  No influential studies detected")
    }
  }

  class(loo_df) <- c("cnma_loo", "data.frame")
  return(loo_df)
}

#' Print Leave-One-Out Results
#' @param x cnma_loo object
#' @param digits Number of decimal places (default: 4)
#' @param ... Additional arguments
#' @export
print.cnma_loo <- function(x, digits = 4, ...) {
  cat("Leave-One-Out Sensitivity Analysis\n")
  cat("====================================\n\n")

  # Format numeric columns
  output <- x
  output$tau <- round(output$tau, digits)
  output$tau_change <- round(output$tau_change, digits)
  output$tau_pct_change <- round(output$tau_pct_change, 2)
  output$I2 <- round(output$I2, digits)
  output$I2_change <- round(output$I2_change, digits)

  print(as.data.frame(output), row.names = FALSE)

  cat("\n")
  cat("Influential studies (>10% change in Tau):\n")
  influential <- output[output$influential == TRUE & !is.na(output$influential), ]
  if (nrow(influential) > 0) {
    print(influential[, c("excluded_study", "tau_pct_change")], row.names = FALSE)
  } else {
    cat("  None detected\n")
  }

  cat("\n")
  cat("Note: Influential studies substantially affect heterogeneity estimates\n")
  cat("      when excluded from the analysis.\n")

  invisible(x)
}

#' Assess Publication Bias
#'
#' Comprehensive assessment of publication bias using multiple methods
#' appropriate for network meta-analysis.
#'
#' @param netmeta_obj A netmeta object
#' @param methods Character vector of methods: "visual" (funnel plot),
#'   "comparison" (comparison-adjusted), "egger" (regression test)
#' @return List with publication bias assessment results
#' @export
#' @references
#' Chaimani A, Salanti G (2012). Using network meta-analysis to evaluate
#' the existence of small-study effects in a network of interventions.
#' Research Synthesis Methods, 3(2):161-76.
#'
#' Egger M, Davey Smith G, Schneider M, Minder C (1997). Bias in meta-analysis
#' detected by a simple, graphical test. BMJ, 315(7109):629-34.
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(40)
#' config <- setup_cnma(use_bayesian = FALSE)
#' results <- run_cnma_analysis(data, config = config)
#' pub_bias <- assess_publication_bias(results$results$main_nma)
#' print(pub_bias)
#' }
assess_publication_bias <- function(netmeta_obj,
                                    methods = c("visual", "comparison", "egger")) {
  if (!requireNamespace("netmeta", quietly = TRUE)) {
    .stop_hint("Package 'netmeta' is required.",
               "Install it with: install.packages('netmeta')")
  }

  if (inherits(netmeta_obj, "try-error") || is.null(netmeta_obj)) {
    .stop_hint("Invalid netmeta object provided.")
  }

  results <- list()

  # Visual inspection (funnel plot)
  if ("visual" %in% methods) {
    msg("Generating funnel plot for visual inspection...")
    plot_funnel(netmeta_obj)
    results$visual <- "Funnel plot generated (visual inspection required)"
  }

  # Comparison-adjusted funnel plot
  if ("comparison" %in% methods) {
    results$comparison <- "Comparison-adjusted funnel plot available via plot_funnel()"
  }

  # Egger's test (if available)
  if ("egger" %in% methods) {
    # Note: Egger's test for NMA is complex
    # The standard Egger test is not directly applicable
    results$egger <- list(
      note = "Egger's test for NMA requires specialized implementation",
      recommendation = "Use visual inspection of comparison-adjusted funnel plot",
      reference = "Chaimani & Salanti (2012) Research Synthesis Methods"
    )
  }

  # General assessment
  k <- length(unique(netmeta_obj$studlab))
  results$summary <- list(
    n_studies = k,
    recommendation = ifelse(
      k < 10,
      "Small number of studies - publication bias assessment unreliable",
      "Sufficient studies for publication bias assessment"
    ),
    methods_applied = methods
  )

  class(results) <- "cnma_pub_bias"
  return(results)
}

#' Print Publication Bias Assessment
#' @param x cnma_pub_bias object
#' @param ... Additional arguments
#' @export
print.cnma_pub_bias <- function(x, ...) {
  cat("Publication Bias Assessment\n")
  cat("============================\n\n")

  cat("Summary:\n")
  cat(sprintf("  Number of studies: %d\n", x$summary$n_studies))
  cat(sprintf("  Recommendation: %s\n\n", x$summary$recommendation))

  cat("Methods applied:\n")
  for (method in x$summary$methods_applied) {
    cat(sprintf("  - %s\n", method))
  }
  cat("\n")

  if (!is.null(x$visual)) {
    cat("Visual inspection:\n")
    cat(sprintf("  %s\n\n", x$visual))
  }

  if (!is.null(x$egger)) {
    cat("Egger's test:\n")
    cat(sprintf("  Note: %s\n", x$egger$note))
    cat(sprintf("  Recommendation: %s\n", x$egger$recommendation))
    cat(sprintf("  Reference: %s\n\n", x$egger$reference))
  }

  cat("Interpretation:\n")
  cat("  - Asymmetry in funnel plot may indicate publication bias\n")
  cat("  - However, asymmetry can also result from heterogeneity\n")
  cat("  - Use comparison-adjusted funnel plot for NMA\n")
  cat("  - Consider contacting authors for unpublished data\n")

  invisible(x)
}

#' Compare Model Assumptions
#'
#' Compares fixed-effect vs random-effects models and evaluates
#' model assumptions.
#'
#' @param data Original data frame
#' @param config Configuration object
#' @return List comparing models
#' @export
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(30)
#' config <- setup_cnma(use_bayesian = FALSE)
#' comparison <- compare_models(data, config)
#' print(comparison)
#' }
compare_models <- function(data, config = setup_cnma()) {
  if (!requireNamespace("netmeta", quietly = TRUE)) {
    .stop_hint("Package 'netmeta' is required.",
               "Install it with: install.packages('netmeta')")
  }

  # Clean data
  data <- cnma_clean_data(data)
  validate_cnma_input(data)

  # Auto-select reference
  ref_treatment <- names(sort(table(c(data$treat1, data$treat2)),
                              decreasing = TRUE))[1]

  # Fixed-effect model
  fixed_model <- .safe_try(
    netmeta::netmeta(
      TE = TE,
      seTE = seTE,
      treat1 = treat1,
      treat2 = treat2,
      studlab = studlab,
      data = data,
      sm = config$sm,
      fixed = TRUE,
      random = FALSE,
      reference.group = ref_treatment
    ),
    context = "Fixed-effect model"
  )

  # Random-effects model
  random_model <- .safe_try(
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
    context = "Random-effects model"
  )

  results <- list(
    fixed_model = fixed_model,
    random_model = random_model
  )

  if (!inherits(random_model, "try-error")) {
    results$heterogeneity <- list(
      tau = random_model$tau,
      tau2 = random_model$tau^2,
      I2 = random_model$I2.random,
      Q = random_model$Q,
      p_value = random_model$pval.Q,
      recommendation = ifelse(
        random_model$I2.random > 0.40 || random_model$pval.Q < 0.10,
        "Substantial heterogeneity - use random-effects model",
        "Low heterogeneity - either model acceptable"
      )
    )
  }

  class(results) <- "cnma_model_comparison"
  return(results)
}

#' Print Model Comparison
#' @param x cnma_model_comparison object
#' @param ... Additional arguments
#' @export
print.cnma_model_comparison <- function(x, ...) {
  cat("Model Comparison: Fixed vs Random Effects\n")
  cat("===========================================\n\n")

  if (!is.null(x$heterogeneity)) {
    cat("Heterogeneity Assessment:\n")
    cat(sprintf("  Tau²: %.4f\n", x$heterogeneity$tau2))
    cat(sprintf("  Tau: %.4f\n", x$heterogeneity$tau))
    cat(sprintf("  I²: %.1f%%\n", x$heterogeneity$I2 * 100))
    cat(sprintf("  Q: %.2f (p = %.4f)\n",
                x$heterogeneity$Q, x$heterogeneity$p_value))
    cat(sprintf("\n  Recommendation: %s\n\n", x$heterogeneity$recommendation))
  }

  cat("Interpretation of I² values:\n")
  cat("  0% - 40%:   Low heterogeneity\n")
  cat("  30% - 60%:  Moderate heterogeneity\n")
  cat("  50% - 90%:  Substantial heterogeneity\n")
  cat("  75% - 100%: Considerable heterogeneity\n\n")

  cat("Reference: Higgins JPT, Thompson SG (2002). Statistics in Medicine, 21:1539-58.\n")

  invisible(x)
}
