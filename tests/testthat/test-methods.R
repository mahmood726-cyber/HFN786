# Tests for S3 methods

test_that("print.cnma works", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(n_studies = 15, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
  results <- run_cnma_analysis(data, config = config)

  expect_output(print(results), "CNMA Results")
  expect_output(print(results), "Summary measure")
  expect_output(print(results), "Reference")
})

test_that("summary.cnma works", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(n_studies = 15, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
  results <- run_cnma_analysis(data, config = config)

  expect_output(summary(results), "CNMA Analysis Summary")
  expect_output(summary(results), "Network characteristics")
  expect_output(summary(results), "Studies")
  expect_output(summary(results), "Treatments")
})

test_that("summary.cnma shows heterogeneity", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(n_studies = 20, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
  results <- run_cnma_analysis(data, config = config)

  expect_output(summary(results), "Heterogeneity")
  expect_output(summary(results), "Tau")
  expect_output(summary(results), "I²")
})
