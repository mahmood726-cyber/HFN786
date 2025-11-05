# Tests for sensitivity analysis functions

test_that("leave_one_out_analysis identifies influential studies", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(15, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)

  loo <- leave_one_out_analysis(data, config, progress = FALSE)

  expect_s3_class(loo, "cnma_loo")
  expect_true("excluded_study" %in% names(loo))
  expect_true("tau" %in% names(loo))
  expect_true("influential" %in% names(loo))
  expect_true(nrow(loo) > 1)  # Should have full model + LOO results
})

test_that("compare_models returns comparison results", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(20, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE)

  comparison <- compare_models(data, config)

  expect_s3_class(comparison, "cnma_model_comparison")
  expect_true("heterogeneity" %in% names(comparison))
})

test_that("assess_publication_bias generates assessment", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(25, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
  results <- run_cnma_analysis(data, config = config)

  pub_bias <- assess_publication_bias(
    results$results$main_nma,
    methods = c("comparison", "egger")
  )

  expect_s3_class(pub_bias, "cnma_pub_bias")
  expect_true("summary" %in% names(pub_bias))
})

test_that("print methods work for sensitivity objects", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(15, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)

  loo <- leave_one_out_analysis(data, config, progress = FALSE)
  comparison <- compare_models(data, config)

  expect_output(print(loo), "Leave-One-Out")
  expect_output(print(comparison), "Model Comparison")
})
