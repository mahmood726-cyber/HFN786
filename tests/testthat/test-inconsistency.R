# Tests for inconsistency assessment functions

test_that("assess_inconsistency works with all methods", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(30, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
  results <- run_cnma_analysis(data, config = config)

  inconsistency <- assess_inconsistency(
    results$results$main_nma,
    methods = c("global", "local", "design")
  )

  expect_s3_class(inconsistency, "cnma_inconsistency")
  expect_true("global" %in% names(inconsistency))
})

test_that("assess_inconsistency handles global method", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(25, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
  results <- run_cnma_analysis(data, config = config)

  inconsistency <- assess_inconsistency(
    results$results$main_nma,
    methods = "global"
  )

  expect_true(!is.null(inconsistency$global))
  expect_true("Q" %in% names(inconsistency$global))
  expect_true("I2" %in% names(inconsistency$global))
  expect_true("tau2" %in% names(inconsistency$global))
})

test_that("assess_transitivity works with covariates", {
  data <- simulate_cnma_data(30, seed = 42)

  transitivity <- assess_transitivity(
    data,
    variables = c("age_mean", "female_pct"),
    by_comparison = TRUE
  )

  expect_s3_class(transitivity, "cnma_transitivity")
  expect_true("age_mean" %in% names(transitivity))
  expect_true("female_pct" %in% names(transitivity))
})

test_that("assess_transitivity handles missing variables gracefully", {
  data <- simulate_cnma_data(20, seed = 42)

  expect_error(
    assess_transitivity(data, variables = c("nonexistent_var")),
    "not found in data"
  )
})

test_that("assess_transitivity auto-detects numeric variables", {
  data <- simulate_cnma_data(20, seed = 42)

  transitivity <- assess_transitivity(data, variables = NULL)

  # Should auto-detect and work or return NULL
  expect_true(is.null(transitivity) || inherits(transitivity, "cnma_transitivity"))
})

test_that("print method works for inconsistency", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(25, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
  results <- run_cnma_analysis(data, config = config)

  inconsistency <- assess_inconsistency(
    results$results$main_nma,
    methods = "global"
  )

  expect_output(print(inconsistency), "Inconsistency Assessment")
})
