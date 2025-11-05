# =========================================================
# Integrated AI-Powered Workflow with Rules Validation
# Complete analysis with automated quality assurance
# =========================================================

#' Run AI-Powered Comprehensive NMA with Rules Validation
#'
#' Performs complete network meta-analysis with automated rules-based validation,
#' AI-powered quality checks, recommendations, and interpretation. Combines all
#' advanced features into a single intelligent workflow.
#'
#' @param data Data frame with NMA data
#' @param ref_treatment Reference treatment (NULL for automatic)
#' @param config Configuration object from setup_cnma()
#' @param ai_model LLama 3 model name (default "llama3")
#' @param validate_rules Run rules validation (default TRUE)
#' @param generate_scenarios Test against scenario database (default FALSE)
#' @param ai_interpretation Generate AI interpretation (default TRUE)
#' @param output_dir Directory for outputs (NULL for no saving)
#' @return Comprehensive results with validation and AI insights
#' @export
#' @examples
#' \dontrun{
#' # Install and start Ollama first:
#' # 1. Download from https://ollama.com/download
#' # 2. ollama pull llama3
#' # 3. ollama serve
#'
#' data <- simulate_cnma_data(40)
#'
#' # Complete AI-powered analysis
#' results <- run_ai_powered_nma(
#'   data,
#'   ai_model = "llama3",
#'   validate_rules = TRUE,
#'   ai_interpretation = TRUE
#' )
#'
#' # Access components
#' print(results$validation_report)
#' print(results$ai_quality_assessment)
#' print(results$ai_interpretation)
#' print(results$nma_results)
#' }
run_ai_powered_nma <- function(data,
                               ref_treatment = NULL,
                               config = setup_cnma(use_bayesian = FALSE),
                               ai_model = "llama3",
                               validate_rules = TRUE,
                               generate_scenarios = FALSE,
                               ai_interpretation = TRUE,
                               output_dir = NULL) {

  cat("\n")
  cat("========================================\n")
  cat("AI-POWERED NMA WITH RULES VALIDATION\n")
  cat("Intelligent Quality Assurance Workflow\n")
  cat("========================================\n\n")

  results <- list()

  # 1. Configure AI assistant if requested
  if (ai_interpretation) {
    msg("Step 1: Configuring local AI assistant (LLama 3)...")

    ai_config <- .safe_try({
      configure_llama3(
        model = ai_model,
        temperature = 0.7,
        seed = 42,
        check_ollama = TRUE
      )
    }, context = "AI configuration", silent = TRUE)

    if (inherits(ai_config, "try-error")) {
      msg("  WARNING: AI assistant unavailable. Continuing without AI features.")
      msg("  Install Ollama and pull model to enable AI: ollama pull %s", ai_model)
      ai_interpretation <- FALSE
    } else {
      results$ai_config <- ai_config
    }
    cat("\n")
  }

  # 2. AI-powered data quality check
  if (ai_interpretation) {
    msg("Step 2: AI-powered data quality assessment...")
    results$ai_quality_assessment <- ai_quality_check(data, verbose = TRUE)
    cat("\n")
  }

  # 3. Rules-based validation
  if (validate_rules) {
    msg("Step 3: Comprehensive rules validation (500+ rules)...")

    engine <- initialize_rules_engine()
    results$rules_engine <- engine

    validation_report <- run_rules_validation(
      data = data,
      nma_results = NULL,
      engine = engine,
      severity_threshold = "warning",
      ai_assist = ai_interpretation
    )

    results$validation_report <- validation_report

    # Stop if critical errors
    if (validation_report$n_errors > 0) {
      msg("\nCRITICAL ERRORS DETECTED!")
      msg("Please address %d error(s) before proceeding.", validation_report$n_errors)
      msg("See results$validation_report for details.\n")

      # Optionally continue anyway
      response <- readline("Continue analysis anyway? (yes/no): ")
      if (!tolower(response) %in% c("yes", "y")) {
        return(results)
      }
    }
    cat("\n")
  }

  # 4. Main NMA analysis
  msg("Step 4: Running comprehensive network meta-analysis...")

  nma_results <- run_comprehensive_nma(
    data = data,
    ref_treatment = ref_treatment,
    config = config,
    output_dir = output_dir,
    run_sensitivity = TRUE,
    generate_plots = TRUE,
    generate_report = TRUE
  )

  results$nma_results <- nma_results
  cat("\n")

  # 5. Post-analysis validation
  if (validate_rules) {
    msg("Step 5: Post-analysis rules validation...")

    post_validation <- run_rules_validation(
      data = data,
      nma_results = nma_results,
      engine = results$rules_engine,
      severity_threshold = "info",
      ai_assist = FALSE  # Already got AI recommendations
    )

    results$post_validation <- post_validation
    cat("\n")
  }

  # 6. AI-powered interpretation
  if (ai_interpretation) {
    msg("Step 6: AI-powered result interpretation...")

    ai_interp <- ai_interpret_results(
      nma_results,
      context = "general clinical",
      target_audience = "clinical"
    )

    results$ai_interpretation <- ai_interp
    cat("\n")
  }

  # 7. AI analysis recommendations
  if (ai_interpretation) {
    msg("Step 7: AI recommendations for further analyses...")

    ai_rec <- ai_recommend_analysis(
      data,
      research_question = "Comparative effectiveness of treatments",
      constraints = "Standard timeline"
    )

    results$ai_recommendations <- ai_rec
    cat("\n")
  }

  # 8. Scenario testing (optional)
  if (generate_scenarios) {
    msg("Step 8: Scenario database testing...")

    scenarios <- generate_scenario_database(n_scenarios = 1000, seed = 42)
    scenario_results <- run_scenario_testing(scenarios, verbose = TRUE)

    results$scenario_database <- scenarios
    results$scenario_test_results <- scenario_results
    cat("\n")
  }

  # 9. Generate comprehensive report
  msg("Step 9: Generating comprehensive quality report...")

  quality_report <- .generate_quality_report(results)
  results$quality_report <- quality_report

  if (!is.null(output_dir)) {
    report_file <- file.path(output_dir, "quality_assurance_report.txt")
    writeLines(quality_report, report_file)
    msg("  Quality report saved: %s", report_file)
  }

  cat("\n")

  # Final summary
  cat("========================================\n")
  cat("ANALYSIS COMPLETE!\n")
  cat("========================================\n\n")

  cat("Quality Assurance Summary:\n")
  if (!is.null(results$validation_report)) {
    cat(sprintf("  Rules checked: %d\n", results$validation_report$n_rules_checked))
    cat(sprintf("  Errors: %d\n", results$validation_report$n_errors))
    cat(sprintf("  Warnings: %d\n", results$validation_report$n_warnings))
  }

  if (!is.null(results$nma_results)) {
    nma <- results$nma_results$main_results$results$main_nma
    cat(sprintf("\nMain Results:\n"))
    cat(sprintf("  Studies: %d\n", length(unique(data$studlab))))
    cat(sprintf("  Treatments: %d\n", length(unique(c(data$treat1, data$treat2)))))
    cat(sprintf("  Heterogeneity (I²): %.1f%%\n", nma$I2.random * 100))
  }

  if (ai_interpretation) {
    cat("\nAI-Powered Insights: Available\n")
    cat("  - Data quality assessment\n")
    cat("  - Result interpretation\n")
    cat("  - Analysis recommendations\n")
  }

  if (!is.null(output_dir)) {
    cat(sprintf("\nAll outputs saved to: %s\n", output_dir))
  }

  cat("\n✓ High-quality, publication-ready analysis complete!\n\n")

  class(results) <- "cnma_ai_powered"

  return(results)
}

#' AI-Assisted Troubleshooting
#'
#' When an analysis fails, use AI to diagnose the issue and provide solutions.
#'
#' @param error_object Error object from failed analysis
#' @param data Data used in analysis
#' @param context Additional context about what you were trying to do
#' @return AI diagnosis and recommendations
#' @export
#' @examples
#' \dontrun{
#' data <- data.frame(studlab = 1:5, treat1 = "A", treat2 = "B",
#'                    TE = c(1, 2, NA, 4, 5), seTE = c(0.5, 0.6, 0.7, 0.8, 0.9))
#'
#' result <- tryCatch(
#'   run_cnma_analysis(data),
#'   error = function(e) e
#' )
#'
#' if (inherits(result, "error")) {
#'   diagnosis <- ai_troubleshoot(result, data, "Running basic NMA")
#' }
#' }
ai_troubleshoot <- function(error_object, data, context = "NMA analysis") {

  if (is.null(.cnma_env$ai_config)) {
    msg("AI assistant not configured.")
    msg("Configure with: configure_llama3()")
    return(NULL)
  }

  error_msg <- if (inherits(error_object, "error")) {
    error_object$message
  } else {
    as.character(error_object)
  }

  cat("\n")
  cat("========================================\n")
  cat("AI-ASSISTED TROUBLESHOOTING\n")
  cat("========================================\n\n")

  # Gather diagnostic information
  data_summary <- sprintf(
    "Data characteristics:\n- %d rows\n- %d studies\n- %d treatments\n- Missing values: %d",
    nrow(data),
    length(unique(data$studlab)),
    length(unique(c(data$treat1, data$treat2))),
    sum(is.na(data))
  )

  # Get AI diagnosis
  diagnosis <- ai_diagnose_error(
    error_message = error_msg,
    context = paste(context, "\n\n", data_summary)
  )

  # Run rules validation to identify issues
  msg("\nRunning rules validation to identify issues...")
  validation <- run_rules_validation(
    data,
    severity_threshold = "error",
    ai_assist = FALSE
  )

  # Combine insights
  result <- list(
    error = error_msg,
    ai_diagnosis = diagnosis,
    validation_issues = validation$violations,
    data_summary = data_summary
  )

  class(result) <- "cnma_ai_troubleshoot"

  return(result)
}

# ========== Internal Helper Functions ==========

.generate_quality_report <- function(results) {
  report <- character()

  report <- c(report, "========================================")
  report <- c(report, "COMPREHENSIVE QUALITY ASSURANCE REPORT")
  report <- c(report, "========================================")
  report <- c(report, "")
  report <- c(report, paste("Generated:", Sys.time()))
  report <- c(report, "")

  # Validation results
  if (!is.null(results$validation_report)) {
    report <- c(report, "RULES VALIDATION")
    report <- c(report, "----------------")
    report <- c(report, sprintf("Rules checked: %d", results$validation_report$n_rules_checked))
    report <- c(report, sprintf("Errors: %d", results$validation_report$n_errors))
    report <- c(report, sprintf("Warnings: %d", results$validation_report$n_warnings))
    report <- c(report, sprintf("Info: %d", results$validation_report$n_info))
    report <- c(report, "")

    if (nrow(results$validation_report$violations) > 0) {
      report <- c(report, "Top Issues:")
      top_issues <- head(results$validation_report$violations, 10)
      for (i in 1:nrow(top_issues)) {
        report <- c(report, sprintf("  - [%s] %s",
                                   top_issues$severity[i],
                                   top_issues$message[i]))
      }
      report <- c(report, "")
    }
  }

  # AI assessment
  if (!is.null(results$ai_quality_assessment)) {
    report <- c(report, "AI QUALITY ASSESSMENT")
    report <- c(report, "---------------------")
    report <- c(report, results$ai_quality_assessment$ai_assessment)
    report <- c(report, "")
  }

  # Analysis results
  if (!is.null(results$nma_results)) {
    report <- c(report, "ANALYSIS RESULTS")
    report <- c(report, "----------------")

    nma <- results$nma_results$main_results$results$main_nma

    report <- c(report, sprintf("Studies: %d",
                               length(unique(results$nma_results$data_clean$studlab))))
    report <- c(report, sprintf("Treatments: %d",
                               length(unique(c(results$nma_results$data_clean$treat1,
                                             results$nma_results$data_clean$treat2)))))
    report <- c(report, sprintf("Heterogeneity (I²): %.1f%%", nma$I2.random * 100))
    report <- c(report, sprintf("Between-study SD (Tau): %.4f", nma$tau))
    report <- c(report, "")
  }

  # AI interpretation
  if (!is.null(results$ai_interpretation)) {
    report <- c(report, "AI INTERPRETATION")
    report <- c(report, "-----------------")
    report <- c(report, results$ai_interpretation$interpretation)
    report <- c(report, "")
  }

  # AI recommendations
  if (!is.null(results$ai_recommendations)) {
    report <- c(report, "AI RECOMMENDATIONS")
    report <- c(report, "------------------")
    report <- c(report, results$ai_recommendations$recommendations)
    report <- c(report, "")
  }

  report <- c(report, "========================================")
  report <- c(report, "END OF REPORT")
  report <- c(report, "========================================")

  paste(report, collapse = "\n")
}

#' Print AI-Powered NMA Results
#' @param x AI-powered NMA results object
#' @param ... Additional arguments
#' @export
print.cnma_ai_powered <- function(x, ...) {
  cat("<AI-Powered CNMA Results>\n\n")

  cat("Components:\n")
  cat(sprintf("  - NMA Results: %s\n",
             ifelse(!is.null(x$nma_results), "✓", "✗")))
  cat(sprintf("  - Rules Validation: %s\n",
             ifelse(!is.null(x$validation_report), "✓", "✗")))
  cat(sprintf("  - AI Quality Assessment: %s\n",
             ifelse(!is.null(x$ai_quality_assessment), "✓", "✗")))
  cat(sprintf("  - AI Interpretation: %s\n",
             ifelse(!is.null(x$ai_interpretation), "✓", "✗")))
  cat(sprintf("  - AI Recommendations: %s\n",
             ifelse(!is.null(x$ai_recommendations), "✓", "✗")))
  cat("\n")

  if (!is.null(x$validation_report)) {
    cat("Validation Summary:\n")
    cat(sprintf("  Errors: %d, Warnings: %d, Info: %d\n",
               x$validation_report$n_errors,
               x$validation_report$n_warnings,
               x$validation_report$n_info))
    cat("\n")
  }

  cat("Access components with $nma_results, $validation_report, etc.\n")

  invisible(x)
}

#' Print AI Troubleshooting Results
#' @param x AI troubleshooting object
#' @param ... Additional arguments
#' @export
print.cnma_ai_troubleshoot <- function(x, ...) {
  cat("<AI Troubleshooting Results>\n\n")

  cat("Error:\n")
  cat(x$error)
  cat("\n\n")

  cat("Data Summary:\n")
  cat(x$data_summary)
  cat("\n\n")

  if (nrow(x$validation_issues) > 0) {
    cat(sprintf("Validation Issues Found: %d\n", nrow(x$validation_issues)))
    cat("See $validation_issues for details\n\n")
  }

  cat("AI Diagnosis:\n")
  cat(x$ai_diagnosis)
  cat("\n")

  invisible(x)
}
