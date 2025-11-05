# Tests for analysis functions

test_that("run_cnma_analysis completes successfully", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(n_studies = 20, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)

  expect_silent({
    results <- run_cnma_analysis(data, config = config)
  })

  expect_s3_class(results, "cnma")
  expect_true("results" %in% names(results))
  expect_true("config" %in% names(results))
  expect_true("ref_treatment" %in% names(results))
})

test_that("run_cnma_analysis auto-selects reference", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(n_studies = 15, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)

  results <- run_cnma_analysis(data, config = config)

  expect_type(results$ref_treatment, "character")
  expect_true(nchar(results$ref_treatment) > 0)
})

test_that("run_cnma_analysis accepts specified reference", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(n_studies = 15, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)

  results <- run_cnma_analysis(
    data,
    ref_treatment = "Placebo",
    config = config
  )

  expect_equal(results$ref_treatment, "Placebo")
})

test_that("run_cnma_analysis validates config", {
  data <- simulate_cnma_data(n_studies = 15, seed = 42)

  expect_error(
    run_cnma_analysis(data, config = list()),
    "cnma_config object"
  )
})

test_that("run_cnma_analysis handles invalid reference", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(n_studies = 15, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)

  expect_error(
    run_cnma_analysis(data, ref_treatment = "NonExistentTreatment", config = config),
    "not found in data"
  )
})

test_that("cnma_quickstart works", {
  skip_if_not_installed("netmeta")

  expect_silent({
    results <- cnma_quickstart(n_studies = 15)
  })

  expect_s3_class(results, "cnma")
})

test_that("cnma_quickstart accepts custom data", {
  skip_if_not_installed("netmeta")

  data <- data.frame(
    studlab = c("S1", "S1", "S2"),
    treat1 = c("A", "A", "A"),
    treat2 = c("B", "C", "B"),
    TE = c(0.5, 0.3, 0.4),
    seTE = c(0.1, 0.2, 0.15)
  )

  results <- cnma_quickstart(data = data)

  expect_s3_class(results, "cnma")
  expect_equal(nrow(results$data), 3)
})

test_that("analysis results contain netmeta object", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(n_studies = 20, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)

  results <- run_cnma_analysis(data, config = config)

  expect_false(inherits(results$results$main_nma, "try-error"))
  expect_s3_class(results$results$main_nma, "netmeta")
})
