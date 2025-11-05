# =========================================================
# Comprehensive Analysis Wrapper
# Complete publication-ready analysis in one function
# =========================================================

#' Run Comprehensive Publication-Ready NMA
#'
#' Performs a complete network meta-analysis with all journal-quality features
#' including treatment rankings, inconsistency assessment, transitivity evaluation,
#' sensitivity analyses, and publication-quality visualizations. Perfect for
#' journal submissions following PRISMA-NMA guidelines.
#'
#' @param data Data frame with network meta-analysis data
#' @param ref_treatment Reference treatment (NULL for automatic selection)
#' @param config Configuration object from setup_cnma()
#' @param output_dir Directory to save all outputs (NULL for no saving)
#' @param transitivity_vars Variables for transitivity assessment
#' @param run_sensitivity Run sensitivity analyses (LOO, publication bias)
#' @param generate_plots Generate all publication plots
#' @param generate_report Generate PRISMA-NMA compliance report
#' @return List with comprehensive results
#' @export
#' @examples
#' \donttest{
#' # Complete analysis in one call
#' data <- simulate_cnma_data(40, seed = 123)
#'
#' results <- run_comprehensive_nma(
#'   data = data,
#'   ref_treatment = "Placebo",
#'   output_dir = "nma_outputs",
#'   generate_plots = TRUE,
#'   generate_report = TRUE
#' )
#'
#' # Access all components
#' print(results$main_results)
#' print(results$rankings)
#' print(results$inconsistency)
#' }
run_comprehensive_nma <- function(data,
                                  ref_treatment = NULL,
                                  config = setup_cnma(use_bayesian = FALSE),
                                  output_dir = NULL,
                                  transitivity_vars = NULL,
                                  run_sensitivity = TRUE,
                                  generate_plots = TRUE,
                                  generate_report = TRUE) {

  cat("\n")
  cat("========================================\n")
  cat("COMPREHENSIVE NMA ANALYSIS\n")
  cat("Publication-Ready with PRISMA-NMA\n")
  cat("========================================\n\n")

  # Create output directory if needed
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      msg("Created output directory: %s", output_dir)
    }
  }

  # Initialize results list
  comp_results <- list()

  # 1. Data preparation and validation
  msg("Step 1/10: Data preparation and validation...")
  data_clean <- cnma_clean_data(data)
  validate_cnma_input(data_clean)

  comp_results$data_clean <- data_clean
  comp_results$network_summary <- network_characteristics_summary(data_clean)
  print(comp_results$network_summary)
  cat("\n")

  # 2. Main network meta-analysis
  msg("Step 2/10: Running main network meta-analysis...")
  main_results <- run_cnma_analysis(
    data = data_clean,
    ref_treatment = ref_treatment,
    config = config
  )

  comp_results$main_results <- main_results
  nma <- main_results$results$main_nma

  if (inherits(nma, "try-error")) {
    .stop_hint("Main NMA failed. Cannot proceed with comprehensive analysis.")
  }

  # 3. Treatment rankings
  msg("Step 3/10: Calculating treatment rankings (P-scores)...")
  comp_results$rankings <- calculate_rankings(nma)
  print(comp_results$rankings)
  cat("\n")

  # 4. League table
  msg("Step 4/10: Generating league table...")
  comp_results$league_table <- create_league_table(nma, digits = 2)
  cat("League table created (use print() to view)\n\n")

  # 5. Prediction intervals
  msg("Step 5/10: Calculating prediction intervals...")
  comp_results$prediction_intervals <- calculate_prediction_intervals(nma)
  cat(sprintf("Calculated prediction intervals for %d comparisons\n\n",
              nrow(comp_results$prediction_intervals)))

  # 6. Inconsistency assessment
  msg("Step 6/10: Assessing network inconsistency...")
  comp_results$inconsistency <- assess_inconsistency(
    nma,
    methods = c("global", "local", "design")
  )
  print(comp_results$inconsistency)
  cat("\n")

  # 7. Transitivity assessment
  msg("Step 7/10: Evaluating transitivity assumption...")
  if (is.null(transitivity_vars)) {
    # Auto-detect numeric variables
    numeric_vars <- names(data_clean)[sapply(data_clean, is.numeric)]
    transitivity_vars <- setdiff(numeric_vars, c("TE", "seTE"))
    if (length(transitivity_vars) > 0) {
      transitivity_vars <- transitivity_vars[1:min(3, length(transitivity_vars))]
    }
  }

  if (length(transitivity_vars) > 0) {
    comp_results$transitivity <- assess_transitivity(
      data_clean,
      variables = transitivity_vars,
      by_comparison = TRUE
    )
    cat(sprintf("Transitivity assessed for: %s\n\n",
                paste(transitivity_vars, collapse = ", ")))
  } else {
    msg("No covariates available for transitivity assessment\n")
    comp_results$transitivity <- NULL
  }

  # 8. Sensitivity analyses
  if (run_sensitivity) {
    msg("Step 8/10: Running sensitivity analyses...")

    # Leave-one-out
    msg("  - Leave-one-out analysis...")
    comp_results$leave_one_out <- leave_one_out_analysis(
      data_clean,
      config = config,
      progress = FALSE
    )

    # Publication bias
    msg("  - Publication bias assessment...")
    comp_results$publication_bias <- assess_publication_bias(
      nma,
      methods = c("comparison", "egger")
    )

    # Model comparison
    msg("  - Model comparison (fixed vs random)...")
    comp_results$model_comparison <- compare_models(data_clean, config)

    cat("Sensitivity analyses completed\n\n")
  } else {
    msg("Step 8/10: Skipping sensitivity analyses...\n")
  }

  # 9. Publication-quality plots
  if (generate_plots) {
    msg("Step 9/10: Generating publication-quality plots...")

    if (!is.null(output_dir)) {
      plots_dir <- file.path(output_dir, "figures")
      create_publication_plots(
        nma,
        output_dir = plots_dir,
        reference = main_results$ref_treatment,
        device = "both"
      )
      comp_results$plots_dir <- plots_dir
    } else {
      msg("  Output directory not specified - displaying plots only\n")
      # Display key plots
      plot_network(nma)
      plot_forest(nma, reference = main_results$ref_treatment)
      plot_rankings(nma)
    }
    cat("\n")
  } else {
    msg("Step 9/10: Skipping plot generation...\n")
  }

  # 10. PRISMA-NMA compliance report
  if (generate_report) {
    msg("Step 10/10: Generating PRISMA-NMA compliance report...")

    report_format <- ifelse(!is.null(output_dir), "markdown", "text")
    report_file <- if (!is.null(output_dir)) {
      file.path(output_dir, "prisma_nma_report.md")
    } else {
      NULL
    }

    comp_results$prisma_report <- generate_prisma_report(
      main_results,
      output_format = report_format,
      file = report_file
    )

    if (!is.null(report_file)) {
      msg("  Report saved to: %s", report_file)
    }
    cat("\n")
  } else {
    msg("Step 10/10: Skipping PRISMA-NMA report generation...\n")
  }

  # Summary
  cat("\n")
  cat("========================================\n")
  cat("COMPREHENSIVE ANALYSIS COMPLETE!\n")
  cat("========================================\n\n")

  cat("Results Summary:\n")
  cat(sprintf("  - %d studies analyzed\n",
              length(unique(data_clean$studlab))))
  cat(sprintf("  - %d treatments compared\n",
              length(unique(c(data_clean$treat1, data_clean$treat2)))))
  cat(sprintf("  - Reference treatment: %s\n", main_results$ref_treatment))
  cat(sprintf("  - Heterogeneity (I²): %.1f%%\n", nma$I2.random * 100))
  cat(sprintf("  - Between-study SD (Tau): %.4f\n", nma$tau))

  if (!is.null(output_dir)) {
    cat(sprintf("\nAll outputs saved to: %s\n", output_dir))
    cat("  - Figures: figures/\n")
    if (generate_report) {
      cat("  - PRISMA report: prisma_nma_report.md\n")
    }
  }

  cat("\nAccess results:\n")
  cat("  - results$main_results - Main NMA object\n")
  cat("  - results$rankings - Treatment rankings (P-scores)\n")
  cat("  - results$league_table - Pairwise comparisons\n")
  cat("  - results$inconsistency - Inconsistency assessment\n")
  cat("  - results$prediction_intervals - Prediction intervals\n")
  if (run_sensitivity) {
    cat("  - results$leave_one_out - LOO sensitivity analysis\n")
    cat("  - results$publication_bias - Publication bias assessment\n")
  }
  cat("\n✓ Ready for publication!\n\n")

  # Set class
  class(comp_results) <- "cnma_comprehensive"

  return(comp_results)
}

#' Print Comprehensive CNMA Results
#' @param x cnma_comprehensive object
#' @param ... Additional arguments
#' @export
print.cnma_comprehensive <- function(x, ...) {
  cat("<Comprehensive CNMA Results>\n\n")

  cat("Main Analysis:\n")
  print(x$main_results)
  cat("\n")

  cat("Treatment Rankings (Top 5):\n")
  if (!is.null(x$rankings)) {
    print(head(x$rankings, 5), row.names = FALSE)
  }
  cat("\n")

  cat("Inconsistency:\n")
  if (!is.null(x$inconsistency$global)) {
    cat(sprintf("  Global I²: %.1f%%\n", x$inconsistency$global$I2 * 100))
    cat(sprintf("  Tau: %.4f\n", sqrt(x$inconsistency$global$tau2)))
  }
  cat("\n")

  if (!is.null(x$leave_one_out)) {
    influential <- sum(x$leave_one_out$influential, na.rm = TRUE)
    cat(sprintf("Sensitivity: %d influential study(ies) detected\n\n", influential))
  }

  cat("Use summary() for more details\n")

  invisible(x)
}

#' Summary of Comprehensive CNMA Results
#' @param object cnma_comprehensive object
#' @param ... Additional arguments
#' @export
summary.cnma_comprehensive <- function(object, ...) {
  cat("========================================\n")
  cat("COMPREHENSIVE NMA SUMMARY\n")
  cat("========================================\n\n")

  # Network characteristics
  cat("NETWORK CHARACTERISTICS\n")
  cat("------------------------\n")
  print(object$network_summary)
  cat("\n")

  # Main results
  cat("MAIN RESULTS\n")
  cat("------------\n")
  summary(object$main_results)
  cat("\n")

  # Rankings
  cat("TREATMENT RANKINGS\n")
  cat("------------------\n")
  if (!is.null(object$rankings)) {
    print(object$rankings, row.names = FALSE)
  }
  cat("\n")

  # Inconsistency
  cat("INCONSISTENCY ASSESSMENT\n")
  cat("------------------------\n")
  if (!is.null(object$inconsistency)) {
    print(object$inconsistency)
  }
  cat("\n")

  # Sensitivity
  if (!is.null(object$leave_one_out)) {
    cat("SENSITIVITY ANALYSIS\n")
    cat("--------------------\n")
    influential <- object$leave_one_out[object$leave_one_out$influential == TRUE, ]
    if (nrow(influential) > 0) {
      cat("Influential studies:\n")
      print(influential[, c("excluded_study", "tau_pct_change")], row.names = FALSE)
    } else {
      cat("No influential studies detected\n")
    }
    cat("\n")
  }

  invisible(object)
}
