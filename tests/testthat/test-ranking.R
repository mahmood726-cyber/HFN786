# Tests for treatment ranking functions

test_that("calculate_rankings works with valid netmeta object", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(20, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
  results <- run_cnma_analysis(data, config = config)

  rankings <- calculate_rankings(results$results$main_nma)

  expect_s3_class(rankings, "cnma_ranking")
  expect_s3_class(rankings, "data.frame")
  expect_true("treatment" %in% names(rankings))
  expect_true("p_score" %in% names(rankings))
  expect_true("rank" %in% names(rankings))
  expect_true(all(rankings$p_score >= 0 & rankings$p_score <= 1))
})

test_that("calculate_rankings handles small_values argument", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(15, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
  results <- run_cnma_analysis(data, config = config)

  rankings_undesirable <- calculate_rankings(results$results$main_nma,
                                             small_values = "undesirable")
  rankings_desirable <- calculate_rankings(results$results$main_nma,
                                           small_values = "desirable")

  expect_s3_class(rankings_undesirable, "cnma_ranking")
  expect_s3_class(rankings_desirable, "cnma_ranking")
})

test_that("create_league_table generates matrix", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(20, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
  results <- run_cnma_analysis(data, config = config)

  league <- create_league_table(results$results$main_nma)

  expect_s3_class(league, "cnma_league_table")
  expect_true(is.matrix(league))
  expect_equal(nrow(league), ncol(league))
})

test_that("calculate_prediction_intervals returns data frame", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(20, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
  results <- run_cnma_analysis(data, config = config)

  pred_int <- calculate_prediction_intervals(results$results$main_nma)

  expect_s3_class(pred_int, "cnma_prediction_intervals")
  expect_true("comparison" %in% names(pred_int))
  expect_true("ci_lower" %in% names(pred_int))
  expect_true("ci_upper" %in% names(pred_int))
  expect_true("pi_lower" %in% names(pred_int))
  expect_true("pi_upper" %in% names(pred_int))

  # Prediction intervals should be wider than confidence intervals
  expect_true(all(pred_int$pi_lower <= pred_int$ci_lower, na.rm = TRUE))
  expect_true(all(pred_int$pi_upper >= pred_int$ci_upper, na.rm = TRUE))
})

test_that("print methods work for ranking objects", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(15, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
  results <- run_cnma_analysis(data, config = config)

  rankings <- calculate_rankings(results$results$main_nma)
  league <- create_league_table(results$results$main_nma)
  pred_int <- calculate_prediction_intervals(results$results$main_nma)

  expect_output(print(rankings), "Treatment Rankings")
  expect_output(print(league), "League Table")
  expect_output(print(pred_int), "Prediction Intervals")
})
