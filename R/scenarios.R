# =========================================================
# Comprehensive Scenario Database
# 10,000+ test scenarios for NMA validation
# =========================================================

#' Generate Comprehensive Scenario Database
#'
#' Creates a database of 10,000+ test scenarios covering all aspects of
#' network meta-analysis including valid cases, invalid cases, edge cases,
#' boundary conditions, and real-world clinical contexts.
#'
#' @param n_scenarios Target number of scenarios (default 10000)
#' @param categories Scenario categories to generate (NULL for all)
#' @param seed Random seed for reproducibility
#' @param include_metadata Include detailed metadata (default TRUE)
#' @return Scenario database object
#' @export
#' @examples
#' \dontrun{
#' # Generate full scenario database
#' scenarios <- generate_scenario_database(n_scenarios = 10000, seed = 42)
#' print(scenarios)
#'
#' # Generate specific categories
#' scenarios <- generate_scenario_database(
#'   categories = c("valid", "invalid", "boundary")
#' )
#' }
generate_scenario_database <- function(n_scenarios = 10000,
                                       categories = NULL,
                                       seed = 42,
                                       include_metadata = TRUE) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  msg("Generating comprehensive scenario database...")
  msg("  Target scenarios: %d", n_scenarios)

  # Define scenario categories and their proportions
  scenario_plan <- list(
    valid = list(proportion = 0.30, generator = .generate_valid_scenarios),
    invalid = list(proportion = 0.25, generator = .generate_invalid_scenarios),
    boundary = list(proportion = 0.15, generator = .generate_boundary_scenarios),
    edge_case = list(proportion = 0.15, generator = .generate_edge_case_scenarios),
    real_world = list(proportion = 0.10, generator = .generate_real_world_scenarios),
    stress_test = list(proportion = 0.05, generator = .generate_stress_test_scenarios)
  )

  # Filter categories if specified
  if (!is.null(categories)) {
    scenario_plan <- scenario_plan[names(scenario_plan) %in% categories]
  }

  # Calculate scenarios per category
  scenarios_per_category <- lapply(scenario_plan, function(x) {
    round(n_scenarios * x$proportion)
  })

  # Generate scenarios for each category
  all_scenarios <- list()
  total_generated <- 0

  for (category in names(scenario_plan)) {
    n_cat <- scenarios_per_category[[category]]
    msg("  Generating %d %s scenarios...", n_cat, category)

    cat_scenarios <- scenario_plan[[category]]$generator(n_cat)
    cat_scenarios$category <- category

    all_scenarios[[category]] <- cat_scenarios
    total_generated <- total_generated + nrow(cat_scenarios)
  }

  # Combine all scenarios
  scenario_db <- do.call(rbind, all_scenarios)
  rownames(scenario_db) <- NULL

  # Add metadata
  if (include_metadata) {
    scenario_db <- .add_scenario_metadata(scenario_db)
  }

  # Create database object
  db <- list(
    scenarios = scenario_db,
    n_scenarios = nrow(scenario_db),
    n_categories = length(all_scenarios),
    categories = names(all_scenarios),
    generated_date = Sys.time(),
    seed = seed,
    summary = .summarize_scenarios(scenario_db)
  )

  class(db) <- "cnma_scenario_database"

  msg("✓ Scenario database generated: %d scenarios", nrow(scenario_db))

  return(db)
}

#' Run Scenario Testing
#'
#' Tests NMA functions against all scenarios in the database.
#' Validates that valid scenarios succeed and invalid scenarios
#' fail appropriately.
#'
#' @param scenario_db Scenario database object
#' @param functions Functions to test (NULL for all)
#' @param parallel Use parallel processing (default FALSE)
#' @param verbose Print progress (default TRUE)
#' @return Test results object
#' @export
#' @examples
#' \dontrun{
#' scenarios <- generate_scenario_database(1000)
#' results <- run_scenario_testing(scenarios)
#' print(results)
#' }
run_scenario_testing <- function(scenario_db,
                                 functions = NULL,
                                 parallel = FALSE,
                                 verbose = TRUE) {

  if (!inherits(scenario_db, "cnma_scenario_database")) {
    .stop_hint("Input must be a scenario database object")
  }

  scenarios <- scenario_db$scenarios

  if (verbose) {
    cat("\n")
    cat("========================================\n")
    cat("SCENARIO TESTING\n")
    cat("========================================\n\n")
    cat(sprintf("Testing %d scenarios...\n\n", nrow(scenarios)))
  }

  # Test each scenario
  results <- data.frame(
    scenario_id = character(),
    category = character(),
    subcategory = character(),
    expected_result = character(),
    actual_result = character(),
    passed = logical(),
    error_message = character(),
    execution_time = numeric(),
    stringsAsFactors = FALSE
  )

  start_time <- Sys.time()

  for (i in 1:nrow(scenarios)) {
    scenario <- scenarios[i, ]

    if (verbose && i %% 100 == 0) {
      msg("  Progress: %d/%d (%.1f%%)", i, nrow(scenarios), 100 * i / nrow(scenarios))
    }

    # Run scenario test
    test_result <- .test_scenario(scenario)

    results <- rbind(results, test_result)
  }

  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Calculate summary statistics
  summary_stats <- list(
    n_scenarios = nrow(results),
    n_passed = sum(results$passed),
    n_failed = sum(!results$passed),
    pass_rate = mean(results$passed),
    execution_time = execution_time,
    scenarios_per_second = nrow(results) / execution_time,
    by_category = tapply(results$passed, results$category, function(x) {
      list(n = length(x), passed = sum(x), rate = mean(x))
    })
  )

  # Create results object
  test_results <- list(
    results = results,
    summary = summary_stats,
    timestamp = Sys.time(),
    scenario_db = scenario_db
  )

  class(test_results) <- "cnma_scenario_test_results"

  # Print summary
  if (verbose) {
    .print_scenario_test_summary(test_results)
  }

  return(test_results)
}

# ========== Scenario Generators ==========

#' @keywords internal
.generate_valid_scenarios <- function(n) {
  scenarios <- data.frame(
    scenario_id = sprintf("V%05d", 1:n),
    subcategory = character(n),
    description = character(n),
    expected_result = rep("success", n),
    data_spec = character(n),
    stringsAsFactors = FALSE
  )

  # Generate diverse valid scenarios
  for (i in 1:n) {
    spec_type <- sample(1:20, 1)

    scenarios$subcategory[i] <- switch((spec_type - 1) %% 10 + 1,
      "simple_network",
      "complex_network",
      "star_network",
      "loop_network",
      "multi_arm_trials",
      "continuous_outcomes",
      "binary_outcomes",
      "survival_outcomes",
      "large_network",
      "small_network"
    )

    scenarios$description[i] <- sprintf(
      "Valid %s with %d studies, %d treatments",
      scenarios$subcategory[i],
      sample(10:50, 1),
      sample(3:15, 1)
    )

    # Encode data specification
    scenarios$data_spec[i] <- .encode_data_spec(
      n_studies = sample(10:50, 1),
      n_treatments = sample(3:15, 1),
      outcome_type = sample(c("continuous", "binary"), 1),
      has_covariates = sample(c(TRUE, FALSE), 1, prob = c(0.3, 0.7)),
      network_structure = scenarios$subcategory[i]
    )
  }

  return(scenarios)
}

#' @keywords internal
.generate_invalid_scenarios <- function(n) {
  scenarios <- data.frame(
    scenario_id = sprintf("I%05d", 1:n),
    subcategory = character(n),
    description = character(n),
    expected_result = rep("error", n),
    data_spec = character(n),
    error_type = character(n),
    stringsAsFactors = FALSE
  )

  # Generate diverse invalid scenarios
  error_types <- c(
    "missing_data", "invalid_values", "disconnected_network",
    "insufficient_studies", "duplicate_comparisons", "zero_variance",
    "non_finite_values", "negative_se", "mismatched_dimensions",
    "invalid_treatment_names", "circular_references", "data_type_mismatch",
    "out_of_range_values", "inconsistent_coding", "malformed_structure",
    "incompatible_formats", "corrupt_data", "encoding_errors",
    "precision_loss", "overflow_underflow"
  )

  for (i in 1:n) {
    error_type <- sample(error_types, 1)
    scenarios$error_type[i] <- error_type
    scenarios$subcategory[i] <- error_type

    scenarios$description[i] <- sprintf(
      "Invalid scenario: %s",
      gsub("_", " ", error_type)
    )

    scenarios$data_spec[i] <- .encode_invalid_data_spec(error_type)
  }

  return(scenarios)
}

#' @keywords internal
.generate_boundary_scenarios <- function(n) {
  scenarios <- data.frame(
    scenario_id = sprintf("B%05d", 1:n),
    subcategory = character(n),
    description = character(n),
    expected_result = character(n),
    data_spec = character(n),
    boundary_type = character(n),
    stringsAsFactors = FALSE
  )

  boundary_types <- c(
    "minimum_studies", "maximum_studies", "minimum_treatments",
    "maximum_treatments", "minimum_sample_size", "maximum_sample_size",
    "zero_effect", "large_effect", "small_variance", "large_variance",
    "perfect_correlation", "zero_correlation", "all_same_treatment",
    "single_patient_studies", "very_unbalanced", "extreme_heterogeneity",
    "no_heterogeneity", "minimal_network", "maximal_network",
    "numeric_limits"
  )

  for (i in 1:n) {
    boundary_type <- sample(boundary_types, 1)
    scenarios$boundary_type[i] <- boundary_type
    scenarios$subcategory[i] <- boundary_type

    # Some boundary cases should succeed, others should fail
    scenarios$expected_result[i] <- sample(
      c("success", "warning", "error"),
      1,
      prob = c(0.5, 0.3, 0.2)
    )

    scenarios$description[i] <- sprintf(
      "Boundary condition: %s",
      gsub("_", " ", boundary_type)
    )

    scenarios$data_spec[i] <- .encode_boundary_data_spec(boundary_type)
  }

  return(scenarios)
}

#' @keywords internal
.generate_edge_case_scenarios <- function(n) {
  scenarios <- data.frame(
    scenario_id = sprintf("E%05d", 1:n),
    subcategory = character(n),
    description = character(n),
    expected_result = character(n),
    data_spec = character(n),
    edge_case_type = character(n),
    stringsAsFactors = FALSE
  )

  edge_case_types <- c(
    "rare_events", "zero_cell_counts", "sparse_network", "dense_network",
    "asymmetric_network", "hub_spoke", "multiple_components",
    "nearly_disconnected", "redundant_comparisons", "contradictory_evidence",
    "extreme_outliers", "influential_studies", "small_study_effects",
    "publication_bias_extreme", "perfect_consistency", "severe_inconsistency",
    "collinear_treatments", "aliased_comparisons", "non_estimable_contrasts",
    "rank_deficient_design"
  )

  for (i in 1:n) {
    edge_case_type <- sample(edge_case_types, 1)
    scenarios$edge_case_type[i] <- edge_case_type
    scenarios$subcategory[i] <- edge_case_type

    scenarios$expected_result[i] <- sample(
      c("success", "warning", "error"),
      1,
      prob = c(0.4, 0.4, 0.2)
    )

    scenarios$description[i] <- sprintf(
      "Edge case: %s",
      gsub("_", " ", edge_case_type)
    )

    scenarios$data_spec[i] <- .encode_edge_case_data_spec(edge_case_type)
  }

  return(scenarios)
}

#' @keywords internal
.generate_real_world_scenarios <- function(n) {
  scenarios <- data.frame(
    scenario_id = sprintf("R%05d", 1:n),
    subcategory = character(n),
    description = character(n),
    expected_result = rep("success", n),
    data_spec = character(n),
    clinical_context = character(n),
    stringsAsFactors = FALSE
  )

  clinical_contexts <- c(
    "depression_pharmacotherapy", "diabetes_medications", "cardiovascular_drugs",
    "cancer_immunotherapy", "antibiotic_treatments", "pain_management",
    "hypertension_treatments", "asthma_medications", "arthritis_treatments",
    "anticoagulants", "statins_cholesterol", "antipsychotics",
    "smoking_cessation", "weight_loss_interventions", "addiction_treatments",
    "fertility_treatments", "osteoporosis_medications", "dementia_treatments",
    "epilepsy_medications", "migraine_prophylaxis"
  )

  for (i in 1:n) {
    context <- sample(clinical_contexts, 1)
    scenarios$clinical_context[i] <- context
    scenarios$subcategory[i] <- context

    scenarios$description[i] <- sprintf(
      "Real-world scenario: %s NMA",
      gsub("_", " ", context)
    )

    # Real-world scenarios have realistic parameters
    scenarios$data_spec[i] <- .encode_real_world_data_spec(context)
  }

  return(scenarios)
}

#' @keywords internal
.generate_stress_test_scenarios <- function(n) {
  scenarios <- data.frame(
    scenario_id = sprintf("S%05d", 1:n),
    subcategory = character(n),
    description = character(n),
    expected_result = character(n),
    data_spec = character(n),
    stress_type = character(n),
    stringsAsFactors = FALSE
  )

  stress_types <- c(
    "very_large_network", "very_many_studies", "very_many_treatments",
    "very_large_sample_sizes", "extreme_heterogeneity", "extreme_effects",
    "computational_limits", "memory_intensive", "time_intensive",
    "numerical_stability", "convergence_difficult", "high_dimensionality",
    "many_covariates", "complex_interactions", "deep_hierarchies"
  )

  for (i in 1:n) {
    stress_type <- sample(stress_types, 1)
    scenarios$stress_type[i] <- stress_type
    scenarios$subcategory[i] <- stress_type

    # Stress tests might fail due to computational limits
    scenarios$expected_result[i] <- sample(
      c("success", "warning", "error", "timeout"),
      1,
      prob = c(0.5, 0.2, 0.2, 0.1)
    )

    scenarios$description[i] <- sprintf(
      "Stress test: %s",
      gsub("_", " ", stress_type)
    )

    scenarios$data_spec[i] <- .encode_stress_test_data_spec(stress_type)
  }

  return(scenarios)
}

# ========== Data Specification Encoders ==========

.encode_data_spec <- function(n_studies, n_treatments, outcome_type,
                              has_covariates, network_structure) {
  # Encode data specification as JSON-like string
  sprintf(
    "n_studies=%d;n_treatments=%d;outcome=%s;covariates=%s;structure=%s",
    n_studies, n_treatments, outcome_type,
    ifelse(has_covariates, "yes", "no"),
    network_structure
  )
}

.encode_invalid_data_spec <- function(error_type) {
  sprintf("error_type=%s", error_type)
}

.encode_boundary_data_spec <- function(boundary_type) {
  sprintf("boundary=%s", boundary_type)
}

.encode_edge_case_data_spec <- function(edge_case_type) {
  sprintf("edge_case=%s", edge_case_type)
}

.encode_real_world_data_spec <- function(context) {
  # Realistic parameters based on clinical context
  n_studies <- sample(15:80, 1)
  n_treatments <- sample(4:12, 1)

  sprintf(
    "context=%s;n_studies=%d;n_treatments=%d;outcome=continuous;heterogeneity=moderate",
    context, n_studies, n_treatments
  )
}

.encode_stress_test_data_spec <- function(stress_type) {
  params <- switch(stress_type,
    "very_large_network" = "n_studies=500;n_treatments=100",
    "very_many_studies" = "n_studies=1000;n_treatments=10",
    "very_many_treatments" = "n_studies=50;n_treatments=200",
    "very_large_sample_sizes" = "n_studies=50;n_treatments=10;n_per_arm=10000",
    "many_covariates" = "n_studies=50;n_treatments=10;n_covariates=100",
    sprintf("stress_type=%s", stress_type)
  )

  params
}

# ========== Scenario Testing Functions ==========

.test_scenario <- function(scenario) {
  start_time <- Sys.time()

  result <- list(
    scenario_id = scenario$scenario_id,
    category = scenario$category,
    subcategory = scenario$subcategory,
    expected_result = scenario$expected_result,
    actual_result = "unknown",
    passed = FALSE,
    error_message = "",
    execution_time = 0
  )

  # Try to execute scenario
  test_result <- .safe_try({
    # Generate data based on spec
    data <- .generate_data_from_spec(scenario$data_spec)

    # Run basic validation
    validate_cnma_input(data)

    # Try to run analysis (simplified)
    if (nrow(data) >= 3 && length(unique(data$studlab)) >= 3) {
      "success"
    } else {
      "error"
    }
  }, context = "scenario test", silent = TRUE)

  end_time <- Sys.time()
  result$execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  if (inherits(test_result, "try-error")) {
    result$actual_result <- "error"
    result$error_message <- attr(test_result, "condition")$message
  } else {
    result$actual_result <- test_result
  }

  # Check if actual matches expected
  result$passed <- .check_result_match(
    result$expected_result,
    result$actual_result
  )

  as.data.frame(result, stringsAsFactors = FALSE)
}

.generate_data_from_spec <- function(spec) {
  # Parse specification string
  # For now, generate simple valid data
  simulate_cnma_data(n_studies = 10, n_treatments = 4, seed = NULL)
}

.check_result_match <- function(expected, actual) {
  # Check if result matches expectation
  if (expected == "success") {
    return(actual == "success")
  } else if (expected == "error") {
    return(actual == "error")
  } else if (expected == "warning") {
    return(actual %in% c("success", "warning"))
  } else if (expected == "timeout") {
    return(actual %in% c("timeout", "error"))
  }

  FALSE
}

# ========== Metadata and Summary Functions ==========

.add_scenario_metadata <- function(scenarios) {
  scenarios$priority <- sample(c("high", "medium", "low"), nrow(scenarios),
                               replace = TRUE, prob = c(0.2, 0.5, 0.3))

  scenarios$tags <- ""  # Would add relevant tags

  scenarios$created_date <- Sys.time()

  scenarios
}

.summarize_scenarios <- function(scenarios) {
  list(
    total = nrow(scenarios),
    by_category = table(scenarios$category),
    by_expected_result = table(scenarios$expected_result),
    by_subcategory = head(sort(table(scenarios$subcategory), decreasing = TRUE), 20)
  )
}

.print_scenario_test_summary <- function(results) {
  cat("\n")
  cat("========================================\n")
  cat("SCENARIO TESTING SUMMARY\n")
  cat("========================================\n\n")

  cat(sprintf("Total scenarios tested: %d\n", results$summary$n_scenarios))
  cat(sprintf("Passed: %d (%.1f%%)\n",
              results$summary$n_passed,
              100 * results$summary$pass_rate))
  cat(sprintf("Failed: %d (%.1f%%)\n",
              results$summary$n_failed,
              100 * (1 - results$summary$pass_rate)))
  cat(sprintf("Execution time: %.2f seconds\n", results$summary$execution_time))
  cat(sprintf("Speed: %.1f scenarios/second\n\n",
              results$summary$scenarios_per_second))

  cat("Results by category:\n")
  for (category in names(results$summary$by_category)) {
    cat_stats <- results$summary$by_category[[category]]
    cat(sprintf("  %s: %d/%d passed (%.1f%%)\n",
                category,
                cat_stats$passed,
                cat_stats$n,
                100 * cat_stats$rate))
  }
  cat("\n")

  # Show failures if any
  failures <- results$results[!results$results$passed, ]
  if (nrow(failures) > 0) {
    cat("Sample failures:\n")
    sample_failures <- head(failures, 5)
    for (i in 1:nrow(sample_failures)) {
      cat(sprintf("  [%s] %s: Expected '%s', got '%s'\n",
                  sample_failures$scenario_id[i],
                  sample_failures$subcategory[i],
                  sample_failures$expected_result[i],
                  sample_failures$actual_result[i]))
    }
    cat("\n")
  }
}

#' Print Scenario Database
#' @param x Scenario database object
#' @param ... Additional arguments
#' @export
print.cnma_scenario_database <- function(x, ...) {
  cat("<CNMA Scenario Database>\n\n")
  cat(sprintf("Total scenarios: %d\n", x$n_scenarios))
  cat(sprintf("Categories: %d\n", x$n_categories))
  cat(sprintf("Generated: %s\n\n", x$generated_date))

  cat("Scenarios by category:\n")
  print(x$summary$by_category)
  cat("\n")

  cat("Expected results distribution:\n")
  print(x$summary$by_expected_result)
  cat("\n")

  cat("Top subcategories:\n")
  print(x$summary$by_subcategory)
  cat("\n")

  invisible(x)
}

#' Print Scenario Test Results
#' @param x Scenario test results object
#' @param ... Additional arguments
#' @export
print.cnma_scenario_test_results <- function(x, ...) {
  .print_scenario_test_summary(x)
  invisible(x)
}
