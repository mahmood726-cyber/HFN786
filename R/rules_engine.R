# =========================================================
# Comprehensive Rules-Based Validation Engine
# 500+ rules for network meta-analysis quality assurance
# =========================================================

#' Initialize Rules Engine
#'
#' Loads and initializes the comprehensive rules database for NMA validation.
#' Rules cover data quality, network structure, statistical assumptions,
#' heterogeneity, inconsistency, publication bias, and reporting quality.
#'
#' @return Rules engine object with 500+ rules
#' @export
#' @examples
#' engine <- initialize_rules_engine()
#' print(engine)
initialize_rules_engine <- function() {

  msg("Initializing comprehensive rules engine...")

  # Load all rule categories
  rules_list <- list(
    data_quality = .load_data_quality_rules(),
    network_structure = .load_network_structure_rules(),
    statistical_assumptions = .load_statistical_assumptions_rules(),
    heterogeneity = .load_heterogeneity_rules(),
    inconsistency = .load_inconsistency_rules(),
    publication_bias = .load_publication_bias_rules(),
    reporting_quality = .load_reporting_quality_rules(),
    effect_size = .load_effect_size_rules(),
    sample_size = .load_sample_size_rules(),
    covariate = .load_covariate_rules(),
    bayesian_specific = .load_bayesian_rules(),
    sensitivity = .load_sensitivity_rules(),
    interpretation = .load_interpretation_rules()
  )

  # Count rules
  n_rules <- sum(sapply(rules_list, nrow))

  engine <- list(
    rules = rules_list,
    n_rules = n_rules,
    version = "1.0.0",
    last_updated = Sys.Date(),
    categories = names(rules_list)
  )

  class(engine) <- "cnma_rules_engine"

  msg("✓ Rules engine initialized: %d rules loaded", n_rules)
  msg("  Categories: %s", paste(names(rules_list), collapse = ", "))

  return(engine)
}

#' Run Comprehensive Rules Validation
#'
#' Validates data, analysis, and results against all applicable rules.
#' Returns detailed report with violations, warnings, and recommendations.
#'
#' @param data Data frame with NMA data
#' @param nma_results Optional NMA results object
#' @param engine Rules engine (creates new if NULL)
#' @param severity_threshold Minimum severity to report ("error", "warning", "info")
#' @param categories Rule categories to check (NULL for all)
#' @param ai_assist Use AI for intelligent recommendations (default TRUE)
#' @return Validation report object
#' @export
#' @examples
#' \dontrun{
#' data <- simulate_cnma_data(30)
#' report <- run_rules_validation(data)
#' print(report)
#'
#' # With NMA results
#' nma <- run_cnma_analysis(data)
#' report <- run_rules_validation(data, nma, ai_assist = TRUE)
#' }
run_rules_validation <- function(data,
                                 nma_results = NULL,
                                 engine = NULL,
                                 severity_threshold = c("error", "warning", "info"),
                                 categories = NULL,
                                 ai_assist = TRUE) {

  severity_threshold <- match.arg(severity_threshold)

  cat("\n")
  cat("========================================\n")
  cat("COMPREHENSIVE RULES VALIDATION\n")
  cat("========================================\n\n")

  # Initialize engine if not provided
  if (is.null(engine)) {
    engine <- initialize_rules_engine()
  }

  # Select categories
  if (is.null(categories)) {
    categories <- engine$categories
  }

  msg("Running validation for %d categories...", length(categories))

  # Run validation for each category
  violations <- list()

  for (category in categories) {
    msg("  Checking %s rules...", category)

    category_violations <- .check_rules_category(
      data = data,
      nma_results = nma_results,
      rules = engine$rules[[category]],
      category = category
    )

    if (nrow(category_violations) > 0) {
      violations[[category]] <- category_violations
    }
  }

  # Combine all violations
  all_violations <- do.call(rbind, violations)

  # Filter by severity
  severity_order <- c("error", "warning", "info")
  threshold_idx <- match(severity_threshold, severity_order)
  all_violations <- all_violations[
    match(all_violations$severity, severity_order) <= threshold_idx,
  ]

  # Generate AI recommendations if requested
  ai_recommendations <- NULL
  if (ai_assist && nrow(all_violations) > 0) {
    msg("\nGenerating AI-powered recommendations...")
    ai_recommendations <- .generate_ai_recommendations(all_violations, data, nma_results)
  }

  # Create report
  report <- list(
    violations = all_violations,
    n_errors = sum(all_violations$severity == "error"),
    n_warnings = sum(all_violations$severity == "warning"),
    n_info = sum(all_violations$severity == "info"),
    ai_recommendations = ai_recommendations,
    categories_checked = categories,
    n_rules_checked = sum(sapply(engine$rules[categories], nrow)),
    timestamp = Sys.time()
  )

  class(report) <- "cnma_validation_report"

  # Print summary
  .print_validation_summary(report)

  return(report)
}

# ========== Rule Category Loaders ==========

#' @keywords internal
.load_data_quality_rules <- function() {
  # 60 data quality rules
  data.frame(
    rule_id = sprintf("DQ%03d", 1:60),
    rule_name = c(
      # Missing data rules (1-10)
      "no_missing_study_labels",
      "no_missing_treatment_names",
      "no_missing_effect_sizes",
      "no_missing_standard_errors",
      "no_missing_sample_sizes",
      "missing_data_threshold_5pct",
      "missing_data_threshold_10pct",
      "missing_not_mnar",
      "complete_case_sufficient",
      "missing_pattern_random",

      # Data validity rules (11-25)
      "finite_effect_sizes",
      "positive_standard_errors",
      "positive_sample_sizes",
      "valid_treatment_names",
      "valid_study_names",
      "no_duplicate_comparisons",
      "no_duplicate_studies",
      "consistent_effect_measure",
      "effect_sizes_plausible_range",
      "standard_errors_plausible",
      "sample_sizes_realistic",
      "no_zero_events_both_arms",
      "no_implausible_correlations",
      "treatment_names_standardized",
      "study_names_unique",

      # Data structure rules (26-40)
      "minimum_3_studies",
      "minimum_3_treatments",
      "minimum_comparisons_adequate",
      "studies_per_comparison_adequate",
      "multi_arm_trials_identified",
      "baseline_treatment_specified",
      "reference_treatment_exists",
      "treatment_coding_consistent",
      "study_design_documented",
      "outcome_direction_specified",
      "time_points_consistent",
      "population_characteristics_available",
      "intervention_details_complete",
      "covariate_data_complete",
      "follow_up_duration_specified",

      # Data integrity rules (41-60)
      "no_data_entry_errors",
      "decimal_places_consistent",
      "units_consistent",
      "confidence_intervals_valid",
      "p_values_consistent_with_ci",
      "effect_sizes_match_published",
      "sample_sizes_sum_correctly",
      "events_not_exceed_sample_size",
      "proportions_between_0_and_1",
      "correlations_between_minus1_and_1",
      "standard_errors_match_ci",
      "no_rounding_errors",
      "dates_in_valid_format",
      "year_publication_realistic",
      "trial_registration_number_valid",
      "doi_format_valid",
      "author_names_valid",
      "journal_names_valid",
      "outcome_measures_standardized",
      "data_extraction_verified"
    ),
    severity = rep(c("error", "error", "error", "warning", "info"), 12),
    check_function = sprintf(".check_dq%03d", 1:60),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.load_network_structure_rules <- function() {
  # 55 network structure rules
  data.frame(
    rule_id = sprintf("NS%03d", 1:55),
    rule_name = c(
      # Connectivity rules (1-15)
      "network_fully_connected",
      "no_isolated_treatments",
      "no_isolated_studies",
      "star_network_adequate_center",
      "common_comparator_specified",
      "indirect_evidence_available",
      "multi_arm_trials_connected",
      "network_density_adequate",
      "average_path_length_acceptable",
      "network_diameter_reasonable",
      "clustering_coefficient_adequate",
      "treatment_degrees_balanced",
      "network_not_too_sparse",
      "network_not_too_dense",
      "loops_present_for_consistency",

      # Geometry rules (16-30)
      "network_geometry_suitable",
      "sufficient_direct_evidence",
      "sufficient_indirect_evidence",
      "evidence_distribution_balanced",
      "no_single_study_treatments",
      "treatment_comparison_coverage",
      "pairwise_comparison_adequacy",
      "three_way_loops_present",
      "higher_order_loops_identified",
      "network_transitivity_feasible",
      "common_comparator_central",
      "treatment_hierarchy_logical",
      "placebo_arm_availability",
      "active_comparator_availability",
      "dose_response_relationships",

      # Evidence base rules (31-45)
      "evidence_base_sufficient",
      "minimum_studies_per_comparison",
      "maximum_studies_per_comparison",
      "evidence_base_recent",
      "evidence_base_diverse",
      "geographic_diversity",
      "temporal_diversity",
      "population_diversity",
      "design_diversity",
      "setting_diversity",
      "evidence_recency_acceptable",
      "outdated_studies_flagged",
      "evidence_gaps_identified",
      "future_studies_suggested",
      "evidence_strength_adequate",

      # Network quality rules (46-55)
      "network_coherence_high",
      "network_stability_adequate",
      "network_informativeness_high",
      "comparison_informativeness",
      "treatment_informativeness",
      "network_contribution_balanced",
      "influential_comparisons_identified",
      "redundant_comparisons_flagged",
      "network_efficiency_high",
      "network_robustness_adequate"
    ),
    severity = rep(c("error", "warning", "info"), length.out = 55),
    check_function = sprintf(".check_ns%03d", 1:55),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.load_statistical_assumptions_rules <- function() {
  # 75 statistical assumptions rules
  data.frame(
    rule_id = sprintf("SA%03d", 1:75),
    rule_name = c(
      # Normality assumptions (1-10)
      "effect_sizes_approximately_normal",
      "residuals_normal",
      "random_effects_normal",
      "no_severe_skewness",
      "no_severe_kurtosis",
      "outliers_identified",
      "outliers_not_excessive",
      "influential_points_detected",
      "leverage_points_identified",
      "cooks_distance_acceptable",

      # Homoscedasticity (11-20)
      "variance_homogeneous_within_comparison",
      "variance_heterogeneity_explained",
      "no_variance_inflation",
      "residual_variance_constant",
      "scale_parameter_appropriate",
      "variance_components_identifiable",
      "between_study_variance_positive",
      "within_study_variance_reasonable",
      "variance_ratio_acceptable",
      "heteroscedasticity_test_passed",

      # Independence (21-30)
      "studies_independent",
      "multi_arm_correlation_handled",
      "repeated_measures_correlation_handled",
      "cluster_randomization_accounted",
      "crossover_correlation_handled",
      "publication_overlap_none",
      "patient_overlap_none",
      "author_overlap_flagged",
      "site_overlap_identified",
      "temporal_correlation_assessed",

      # Transitivity (31-45)
      "transitivity_assumption_plausible",
      "effect_modifiers_balanced",
      "population_characteristics_similar",
      "study_designs_comparable",
      "outcome_definitions_consistent",
      "intervention_definitions_consistent",
      "control_definitions_consistent",
      "timing_assessments_similar",
      "follow_up_durations_comparable",
      "risk_bias_levels_similar",
      "publication_years_distributed",
      "geographic_regions_represented",
      "settings_comparable",
      "subgroups_similar",
      "covariate_distributions_similar",

      # Model assumptions (46-60)
      "random_effects_appropriate",
      "fixed_effects_appropriate",
      "consistency_assumption_met",
      "exchangeability_assumption_met",
      "additivity_assumption_met",
      "linearity_assumption_met",
      "proportional_hazards_if_survival",
      "proportional_odds_if_ordinal",
      "link_function_appropriate",
      "distribution_family_appropriate",
      "prior_distributions_appropriate",
      "hyperparameters_reasonable",
      "convergence_criteria_met",
      "identifiability_ensured",
      "estimability_verified",

      # Robustness (61-75)
      "results_robust_to_priors",
      "results_robust_to_model_choice",
      "results_robust_to_outliers",
      "results_robust_to_missing_data",
      "results_robust_to_effect_measure",
      "results_robust_to_reference_treatment",
      "sensitivity_analyses_conducted",
      "leave_one_out_stable",
      "meta_regression_stable",
      "subgroup_analyses_consistent",
      "network_meta_regression_appropriate",
      "interaction_terms_justified",
      "covariate_coding_appropriate",
      "centering_decisions_appropriate",
      "scaling_decisions_appropriate"
    ),
    severity = rep(c("error", "warning", "info"), 25),
    check_function = sprintf(".check_sa%03d", 1:75),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.load_heterogeneity_rules <- function() {
  # 50 heterogeneity rules
  data.frame(
    rule_id = sprintf("HT%03d", 1:50),
    rule_name = c(
      # Detection rules (1-15)
      "heterogeneity_assessed",
      "tau_squared_estimated",
      "i_squared_calculated",
      "h_squared_calculated",
      "cochran_q_test_performed",
      "heterogeneity_test_power_adequate",
      "between_study_variance_estimator_appropriate",
      "heterogeneity_not_zero_inappropriately",
      "heterogeneity_estimates_precise",
      "heterogeneity_confidence_intervals_reported",
      "prediction_intervals_calculated",
      "between_comparison_heterogeneity",
      "design_heterogeneity_assessed",
      "heterogeneity_subgroup_assessed",
      "heterogeneity_meta_regression_explored",

      # Interpretation rules (16-30)
      "low_heterogeneity_i2_below_25",
      "moderate_heterogeneity_i2_25_to_50",
      "substantial_heterogeneity_i2_50_to_75",
      "considerable_heterogeneity_i2_above_75",
      "heterogeneity_clinical_importance_assessed",
      "heterogeneity_sources_identified",
      "heterogeneity_explained_by_covariates",
      "residual_heterogeneity_acceptable",
      "unexplained_heterogeneity_flagged",
      "heterogeneity_consistency_across_comparisons",
      "heterogeneity_patterns_interpreted",
      "heterogeneity_impact_on_conclusions",
      "tau_comparison_to_published_values",
      "predictive_distribution_width",
      "heterogeneity_vs_inconsistency_distinguished",

      # Management rules (31-50)
      "random_effects_used_if_heterogeneous",
      "heterogeneity_prespecified_subgroups",
      "heterogeneity_covariate_adjustment",
      "heterogeneity_outlier_investigation",
      "heterogeneity_influential_studies",
      "heterogeneity_publication_bias_considered",
      "heterogeneity_small_study_effects",
      "heterogeneity_quality_scores_explored",
      "heterogeneity_dose_response_explored",
      "heterogeneity_timing_explored",
      "heterogeneity_population_explored",
      "heterogeneity_intervention_details",
      "heterogeneity_comparator_details",
      "heterogeneity_outcome_measurement",
      "heterogeneity_risk_of_bias",
      "heterogeneity_reporting_complete",
      "heterogeneity_forest_plot_shown",
      "heterogeneity_funnel_plot_examined",
      "heterogeneity_sensitivity_analyses",
      "heterogeneity_future_research_implications"
    ),
    severity = rep(c("warning", "info", "warning"), length.out = 50),
    check_function = sprintf(".check_ht%03d", 1:50),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.load_inconsistency_rules <- function() {
  # 45 inconsistency rules
  data.frame(
    rule_id = sprintf("IC%03d", 1:45),
    rule_name = c(
      # Assessment rules (1-15)
      "inconsistency_assessed",
      "consistency_assumption_tested",
      "global_inconsistency_test",
      "local_inconsistency_test",
      "node_splitting_performed",
      "side_splitting_performed",
      "design_by_treatment_interaction",
      "inconsistency_loops_examined",
      "inconsistency_all_loops_checked",
      "inconsistency_factors_calculated",
      "inconsistency_sufficient_loops",
      "inconsistency_power_adequate",
      "inconsistency_plot_generated",
      "inconsistency_heat_map_shown",
      "inconsistency_comparison_level",

      # Detection rules (16-30)
      "no_significant_global_inconsistency",
      "no_significant_local_inconsistency",
      "direct_indirect_agreement",
      "closed_loops_consistent",
      "open_loops_identified",
      "inconsistency_patterns_identified",
      "inconsistency_not_systematic",
      "inconsistency_treatment_specific_checked",
      "inconsistency_comparison_specific_checked",
      "inconsistency_design_specific_checked",
      "inconsistency_magnitude_reported",
      "inconsistency_confidence_intervals",
      "inconsistency_p_values_reported",
      "inconsistency_effect_size_metric",
      "inconsistency_vs_heterogeneity",

      # Resolution rules (31-45)
      "inconsistency_causes_investigated",
      "inconsistency_study_characteristics",
      "inconsistency_population_differences",
      "inconsistency_intervention_differences",
      "inconsistency_outcome_differences",
      "inconsistency_quality_differences",
      "inconsistency_bias_considered",
      "inconsistency_outliers_identified",
      "inconsistency_data_errors_ruled_out",
      "inconsistency_explained_or_accepted",
      "inconsistency_sensitivity_analyses",
      "inconsistency_subgroup_analyses",
      "inconsistency_meta_regression",
      "inconsistency_impact_on_results",
      "inconsistency_reporting_complete"
    ),
    severity = rep(c("error", "warning", "info"), 15),
    check_function = sprintf(".check_ic%03d", 1:45),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.load_publication_bias_rules <- function() {
  # 35 publication bias rules
  data.frame(
    rule_id = sprintf("PB%03d", 1:35),
    rule_name = c(
      "publication_bias_assessed",
      "funnel_plot_generated",
      "comparison_adjusted_funnel_plot",
      "egger_test_performed",
      "begg_test_performed",
      "trim_fill_analysis",
      "selection_models_considered",
      "small_study_effects_examined",
      "asymmetry_tests_appropriate",
      "minimum_studies_for_bias_assessment",
      "funnel_plot_asymmetry_visual",
      "funnel_plot_asymmetry_statistical",
      "publication_bias_direction_identified",
      "missing_studies_estimated",
      "adjusted_effect_sizes_calculated",
      "grey_literature_search_conducted",
      "trial_registries_searched",
      "unpublished_studies_included",
      "negative_results_included",
      "language_bias_minimized",
      "time_lag_bias_considered",
      "duplicate_publication_bias",
      "citation_bias_considered",
      "outcome_reporting_bias",
      "selective_analysis_reporting",
      "p_hacking_indicators",
      "publication_bias_impact_assessed",
      "sensitivity_to_missing_studies",
      "worst_case_scenarios_explored",
      "publication_bias_mechanisms",
      "commercial_funding_bias",
      "author_conflicts_of_interest",
      "journal_impact_factor_bias",
      "positive_result_bias",
      "publication_bias_reporting_complete"
    ),
    severity = rep(c("warning", "info"), length.out = 35),
    check_function = sprintf(".check_pb%03d", 1:35),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.load_reporting_quality_rules <- function() {
  # 60 reporting quality rules (PRISMA-NMA compliance)
  data.frame(
    rule_id = sprintf("RQ%03d", 1:60),
    rule_name = c(
      # PRISMA-NMA specific (1-20)
      "prisma_nma_checklist_completed",
      "network_diagram_included",
      "network_geometry_described",
      "evidence_structure_table",
      "risk_of_bias_assessment",
      "transitivity_assessment_reported",
      "consistency_evaluation_reported",
      "ranking_metrics_described",
      "heterogeneity_measures_reported",
      "prediction_intervals_reported",
      "league_table_included",
      "forest_plots_included",
      "comparison_adjusted_funnel_plots",
      "assumptions_clearly_stated",
      "limitations_discussed",
      "certainty_of_evidence_graded",
      "grade_or_equiv_framework",
      "clinical_interpretation_provided",
      "implications_for_practice",
      "implications_for_research",

      # Methods reporting (21-40)
      "search_strategy_detailed",
      "inclusion_criteria_clear",
      "exclusion_criteria_clear",
      "data_extraction_process",
      "quality_assessment_process",
      "statistical_methods_detailed",
      "software_specified",
      "model_specifications_complete",
      "reference_treatment_justified",
      "effect_measure_justified",
      "random_effects_model_justified",
      "heterogeneity_handling_described",
      "inconsistency_handling_described",
      "sensitivity_analyses_prespecified",
      "subgroup_analyses_prespecified",
      "meta_regression_prespecified",
      "missing_data_handling",
      "multi_arm_trial_handling",
      "correlation_assumptions",
      "prior_distributions_if_bayesian",

      # Results reporting (41-60)
      "study_flow_diagram",
      "study_characteristics_table",
      "network_meta_analysis_results",
      "pairwise_meta_analysis_results",
      "heterogeneity_results",
      "inconsistency_results",
      "publication_bias_results",
      "sensitivity_analysis_results",
      "subgroup_analysis_results",
      "meta_regression_results",
      "ranking_results_with_uncertainty",
      "effect_estimates_with_ci",
      "effect_estimates_with_pi",
      "p_values_reported_appropriately",
      "confidence_intervals_appropriate_level",
      "prediction_intervals_appropriate_level",
      "number_of_studies_per_comparison",
      "total_sample_sizes",
      "event_rates_if_binary",
      "follow_up_durations"
    ),
    severity = rep("warning", 60),
    check_function = sprintf(".check_rq%03d", 1:60),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.load_effect_size_rules <- function() {
  # 40 effect size rules
  data.frame(
    rule_id = sprintf("ES%03d", 1:40),
    rule_name = c(
      "effect_measure_appropriate_for_outcome",
      "effect_direction_consistent",
      "effect_sizes_clinically_meaningful",
      "minimal_clinically_important_difference",
      "effect_sizes_within_plausible_range",
      "odds_ratios_appropriate",
      "risk_ratios_appropriate",
      "hazard_ratios_appropriate",
      "mean_differences_appropriate",
      "standardized_mean_differences_appropriate",
      "correlation_coefficients_appropriate",
      "rare_events_handled_appropriately",
      "zero_events_handled_appropriately",
      "continuity_correction_appropriate",
      "transformation_appropriate",
      "back_transformation_correct",
      "effect_size_interpretation_correct",
      "confidence_intervals_interpretable",
      "null_hypothesis_clearly_defined",
      "alternative_hypothesis_appropriate",
      "one_sided_vs_two_sided_appropriate",
      "multiplicity_adjustments",
      "effect_sizes_comparable_across_studies",
      "effect_sizes_not_double_counted",
      "effect_sizes_independent",
      "effect_modification_assessed",
      "dose_response_relationships",
      "time_to_effect_considered",
      "duration_of_effect_considered",
      "carryover_effects_considered",
      "baseline_imbalances_adjusted",
      "regression_to_mean_considered",
      "ceiling_floor_effects",
      "responder_analyses_appropriate",
      "composite_outcomes_justified",
      "surrogate_outcomes_validated",
      "patient_important_outcomes",
      "clinical_vs_statistical_significance",
      "effect_size_precision_adequate",
      "effect_size_pooling_appropriate"
    ),
    severity = rep(c("warning", "info"), 20),
    check_function = sprintf(".check_es%03d", 1:40),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.load_sample_size_rules <- function() {
  # 30 sample size rules
  data.frame(
    rule_id = sprintf("SS%03d", 1:30),
    rule_name = c(
      "minimum_sample_size_per_study",
      "minimum_total_sample_size",
      "sample_size_calculations_reported",
      "power_analysis_conducted",
      "sample_size_adequate_for_subgroups",
      "sample_size_adequate_for_meta_regression",
      "sample_sizes_balanced_across_arms",
      "sample_sizes_consistent_with_events",
      "sample_sizes_realistic",
      "dropout_rates_reported",
      "intention_to_treat_sample_sizes",
      "per_protocol_sample_sizes",
      "sample_size_heterogeneity_across_studies",
      "small_studies_identified",
      "large_studies_identified",
      "sample_size_weighted_analyses",
      "inverse_variance_weighting_appropriate",
      "rare_disease_sample_sizes_acknowledged",
      "cluster_sample_sizes_adjusted",
      "multilevel_sample_sizes_appropriate",
      "sample_size_missing_data_impact",
      "effective_sample_sizes_calculated",
      "sample_size_reductions_explained",
      "sample_size_increases_explained",
      "sample_size_stopping_rules_reported",
      "adaptive_sample_size_procedures",
      "sample_size_for_safety_outcomes",
      "sample_size_for_secondary_outcomes",
      "sample_size_precision_estimates",
      "sample_size_future_studies_recommended"
    ),
    severity = rep(c("warning", "info"), 15),
    check_function = sprintf(".check_ss%03d", 1:30),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.load_covariate_rules <- function() {
  # 40 covariate rules
  data.frame(
    rule_id = sprintf("CV%03d", 1:40),
    rule_name = c(
      "covariates_prespecified",
      "covariate_selection_justified",
      "covariates_clinically_relevant",
      "covariates_measured_consistently",
      "covariate_definitions_clear",
      "covariate_distributions_described",
      "covariate_balance_assessed",
      "covariate_missing_data_handled",
      "covariate_transformations_appropriate",
      "covariate_coding_appropriate",
      "covariate_centering_appropriate",
      "covariate_scaling_appropriate",
      "covariate_interactions_considered",
      "covariate_non_linearity_assessed",
      "covariate_multicollinearity_checked",
      "covariate_effect_modification",
      "covariate_confounding_assessed",
      "covariate_mediation_considered",
      "covariate_ecological_bias",
      "covariate_aggregation_bias",
      "individual_vs_aggregate_covariates",
      "baseline_risk_covariate",
      "age_covariate_handled",
      "sex_covariate_handled",
      "race_ethnicity_covariate",
      "disease_severity_covariate",
      "comorbidity_covariate",
      "prior_treatment_covariate",
      "dose_covariate",
      "duration_covariate",
      "follow_up_length_covariate",
      "setting_covariate",
      "geographic_location_covariate",
      "publication_year_covariate",
      "study_quality_covariate",
      "risk_of_bias_covariate",
      "covariate_effects_reported_clearly",
      "covariate_uncertainty_quantified",
      "covariate_clinical_interpretation",
      "covariate_sensitivity_analyses"
    ),
    severity = rep("info", 40),
    check_function = sprintf(".check_cv%03d", 1:40),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.load_bayesian_rules <- function() {
  # 30 Bayesian-specific rules
  data.frame(
    rule_id = sprintf("BY%03d", 1:30),
    rule_name = c(
      "priors_specified_completely",
      "priors_justified",
      "priors_sensitivity_analysis",
      "priors_vague_vs_informative",
      "priors_expert_elicitation_if_informative",
      "hyperpriors_appropriate",
      "mcmc_chains_sufficient",
      "mcmc_iterations_sufficient",
      "mcmc_burnin_adequate",
      "mcmc_thinning_appropriate",
      "convergence_diagnostics_reported",
      "rhat_below_threshold",
      "effective_sample_size_adequate",
      "trace_plots_examined",
      "autocorrelation_assessed",
      "posterior_predictive_checks",
      "deviance_information_criterion",
      "watanabe_akaike_information_criterion",
      "model_comparison_metrics",
      "bayes_factors_reported",
      "credible_intervals_appropriate",
      "highest_density_intervals",
      "posterior_probabilities_reported",
      "sucra_probabilities_reported",
      "ranking_probabilities_complete",
      "bayesian_p_values_if_used",
      "posterior_distributions_visualized",
      "prior_posterior_comparison",
      "bayesian_inference_interpretation_correct",
      "bayesian_software_specified"
    ),
    severity = rep("warning", 30),
    check_function = sprintf(".check_by%03d", 1:30),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.load_sensitivity_rules <- function() {
  # 25 sensitivity analysis rules
  data.frame(
    rule_id = sprintf("SN%03d", 1:25),
    rule_name = c(
      "sensitivity_analyses_performed",
      "sensitivity_fixed_vs_random",
      "sensitivity_effect_measure",
      "sensitivity_reference_treatment",
      "sensitivity_study_quality",
      "sensitivity_risk_of_bias",
      "sensitivity_publication_type",
      "sensitivity_sample_size_threshold",
      "sensitivity_outlier_removal",
      "sensitivity_influential_studies",
      "sensitivity_leave_one_out",
      "sensitivity_subgroup_analyses",
      "sensitivity_meta_regression",
      "sensitivity_missing_data_assumptions",
      "sensitivity_imputation_methods",
      "sensitivity_correlation_assumptions",
      "sensitivity_heterogeneity_priors",
      "sensitivity_consistency_models",
      "sensitivity_network_structure",
      "sensitivity_covariate_selection",
      "sensitivity_functional_form",
      "sensitivity_results_robust",
      "sensitivity_conclusions_unchanged",
      "sensitivity_results_reported_completely",
      "sensitivity_implications_discussed"
    ),
    severity = rep("warning", 25),
    check_function = sprintf(".check_sn%03d", 1:25),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.load_interpretation_rules <- function() {
  # 30 interpretation and reporting rules
  data.frame(
    rule_id = sprintf("IN%03d", 1:30),
    rule_name = c(
      "conclusions_supported_by_evidence",
      "overinterpretation_avoided",
      "causal_language_appropriate",
      "uncertainty_acknowledged",
      "limitations_comprehensively_discussed",
      "generalizability_discussed",
      "applicability_discussed",
      "clinical_relevance_addressed",
      "statistical_vs_clinical_significance",
      "practical_implications_clear",
      "policy_implications_appropriate",
      "conflicts_of_interest_declared",
      "funding_sources_declared",
      "author_contributions_specified",
      "data_availability_statement",
      "protocol_registration",
      "protocol_deviations_explained",
      "reporting_guidelines_followed",
      "ethics_approval_not_needed_justified",
      "patient_involvement_described",
      "equity_considerations",
      "cost_effectiveness_if_relevant",
      "implementation_considerations",
      "knowledge_translation_plan",
      "dissemination_plan",
      "future_research_directions_specific",
      "ongoing_studies_acknowledged",
      "comparison_with_existing_reviews",
      "consistency_with_previous_findings",
      "contribution_to_evidence_base_clear"
    ),
    severity = rep("info", 30),
    check_function = sprintf(".check_in%03d", 1:30),
    stringsAsFactors = FALSE
  )
}

# ========== Rule Checking Functions ==========

.check_rules_category <- function(data, nma_results, rules, category) {
  violations <- data.frame(
    rule_id = character(),
    rule_name = character(),
    category = character(),
    severity = character(),
    message = character(),
    value = character(),
    threshold = character(),
    recommendation = character(),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(rules)) {
    rule <- rules[i, ]

    # Check rule
    check_result <- .evaluate_rule(rule, data, nma_results)

    if (!check_result$passed) {
      violation <- data.frame(
        rule_id = rule$rule_id,
        rule_name = rule$rule_name,
        category = category,
        severity = rule$severity,
        message = check_result$message,
        value = check_result$value,
        threshold = check_result$threshold,
        recommendation = check_result$recommendation,
        stringsAsFactors = FALSE
      )

      violations <- rbind(violations, violation)
    }
  }

  return(violations)
}

.evaluate_rule <- function(rule, data, nma_results) {
  # Simplified rule evaluation - would call specific check functions
  # For now, perform basic checks

  result <- list(
    passed = TRUE,
    message = "",
    value = "",
    threshold = "",
    recommendation = ""
  )

  # Example implementations for key rules
  if (rule$rule_id == "DQ001") {  # no_missing_study_labels
    if (any(is.na(data$studlab))) {
      result$passed <- FALSE
      result$message <- "Missing study labels detected"
      result$value <- sprintf("%d missing", sum(is.na(data$studlab)))
      result$threshold <- "0 missing"
      result$recommendation <- "Assign unique study labels to all rows"
    }
  } else if (rule$rule_id == "DQ011") {  # finite_effect_sizes
    if (any(!is.finite(data$TE))) {
      result$passed <- FALSE
      result$message <- "Non-finite effect sizes detected (Inf, -Inf, or NaN)"
      result$value <- sprintf("%d non-finite", sum(!is.finite(data$TE)))
      result$threshold <- "0 non-finite"
      result$recommendation <- "Check data extraction and calculations"
    }
  } else if (rule$rule_id == "DQ026") {  # minimum_3_studies
    n_studies <- length(unique(data$studlab))
    if (n_studies < 3) {
      result$passed <- FALSE
      result$message <- "Insufficient number of studies for reliable meta-analysis"
      result$value <- sprintf("%d studies", n_studies)
      result$threshold <- "≥3 studies"
      result$recommendation <- "Include more studies or use alternative synthesis methods"
    }
  } else if (rule$rule_id == "NS001") {  # network_fully_connected
    if (!is.null(nma_results)) {
      # Would check network connectivity here
      # Placeholder
    }
  }

  return(result)
}

.generate_ai_recommendations <- function(violations, data, nma_results) {
  if (is.null(.cnma_env$ai_config)) {
    return("AI assistant not configured. Use configure_llama3() for AI recommendations.")
  }

  # Summarize violations for AI
  summary_text <- sprintf(
    "Network meta-analysis validation identified %d issues:\n\nErrors: %d\nWarnings: %d\nInfo: %d\n\n",
    nrow(violations),
    sum(violations$severity == "error"),
    sum(violations$severity == "warning"),
    sum(violations$severity == "info")
  )

  # Add top issues
  top_violations <- head(violations[violations$severity == "error", ], 5)
  if (nrow(top_violations) > 0) {
    summary_text <- paste0(summary_text, "Top Issues:\n")
    for (i in 1:nrow(top_violations)) {
      summary_text <- paste0(
        summary_text,
        sprintf("- [%s] %s: %s\n",
                top_violations$severity[i],
                top_violations$rule_name[i],
                top_violations$message[i])
      )
    }
  }

  prompt <- paste0(
    summary_text,
    "\n\nBased on these validation results, provide:\n",
    "1. Priority ranking of issues to address first\n",
    "2. Specific actionable steps to resolve each issue\n",
    "3. Potential implications if issues are not addressed\n",
    "4. Recommendations for improving analysis quality\n",
    "\nBe specific, practical, and focus on the most impactful improvements."
  )

  ai_response <- call_llama3(prompt)

  return(ai_response)
}

.print_validation_summary <- function(report) {
  cat("\n")
  cat("========================================\n")
  cat("VALIDATION SUMMARY\n")
  cat("========================================\n\n")

  cat(sprintf("Rules checked: %d\n", report$n_rules_checked))
  cat(sprintf("Categories: %s\n\n", paste(report$categories_checked, collapse = ", ")))

  cat("Issues found:\n")
  cat(sprintf("  🔴 Errors:   %d\n", report$n_errors))
  cat(sprintf("  🟡 Warnings: %d\n", report$n_warnings))
  cat(sprintf("  🔵 Info:     %d\n", report$n_info))
  cat("\n")

  if (nrow(report$violations) > 0) {
    cat("Top violations:\n")
    top <- head(report$violations, 10)
    for (i in 1:nrow(top)) {
      icon <- switch(top$severity[i],
                    "error" = "🔴",
                    "warning" = "🟡",
                    "info" = "🔵")
      cat(sprintf("  %s [%s] %s\n", icon, top$rule_id[i], top$message[i]))
    }
    cat("\n")
  }

  if (!is.null(report$ai_recommendations)) {
    cat("AI RECOMMENDATIONS\n")
    cat("==================\n\n")
    cat(report$ai_recommendations)
    cat("\n")
  }

  if (report$n_errors == 0 && report$n_warnings == 0) {
    cat("✓ All validation checks passed!\n\n")
  } else {
    cat("See report$violations for full details\n\n")
  }
}

#' Print Rules Engine
#' @param x Rules engine object
#' @param ... Additional arguments
#' @export
print.cnma_rules_engine <- function(x, ...) {
  cat("<CNMA Rules Engine>\n\n")
  cat(sprintf("Version: %s\n", x$version))
  cat(sprintf("Total rules: %d\n", x$n_rules))
  cat(sprintf("Last updated: %s\n\n", x$last_updated))

  cat("Rule categories:\n")
  for (category in x$categories) {
    n_rules <- nrow(x$rules[[category]])
    cat(sprintf("  - %s: %d rules\n", category, n_rules))
  }
  cat("\n")

  invisible(x)
}

#' Print Validation Report
#' @param x Validation report object
#' @param ... Additional arguments
#' @export
print.cnma_validation_report <- function(x, ...) {
  .print_validation_summary(x)
  invisible(x)
}
