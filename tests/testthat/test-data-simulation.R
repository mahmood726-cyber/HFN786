# Tests for data simulation

test_that("simulate_cnma_data generates valid data", {
  data <- simulate_cnma_data(n_studies = 20, seed = 42)

  expect_s3_class(data, "data.frame")
  expect_true(nrow(data) > 0)

  required_cols <- c("studlab", "treat1", "treat2", "TE", "seTE")
  expect_true(all(required_cols %in% names(data)))
})

test_that("simulate_cnma_data respects seed", {
  data1 <- simulate_cnma_data(n_studies = 20, seed = 123)
  data2 <- simulate_cnma_data(n_studies = 20, seed = 123)
  data3 <- simulate_cnma_data(n_studies = 20, seed = 456)

  expect_equal(data1$TE, data2$TE)
  expect_false(all(data1$TE == data3$TE))
})

test_that("simulate_cnma_data includes covariates", {
  data <- simulate_cnma_data(n_studies = 20, seed = 42)

  expected_covars <- c(
    "year", "is_rct", "study_design", "grade", "rob2",
    "age_mean", "female_pct", "bmi_mean", "charlson", "baseline_risk"
  )

  for (var in expected_covars) {
    expect_true(var %in% names(data),
                info = sprintf("Missing covariate: %s", var))
  }
})

test_that("simulate_cnma_data generates plausible values", {
  data <- simulate_cnma_data(n_studies = 50, seed = 42)

  # Check treatment effects are log-scale
  expect_true(all(abs(data$TE) < 5))

  # Check standard errors are positive and reasonable
  expect_true(all(data$seTE > 0))
  expect_true(all(data$seTE < 1))

  # Check covariates are in plausible ranges
  expect_true(all(data$age_mean > 40 & data$age_mean < 90))
  expect_true(all(data$female_pct >= 0 & data$female_pct <= 1))
  expect_true(all(data$bmi_mean > 15 & data$bmi_mean < 45))
})

test_that("simulate_cnma_data validates input", {
  expect_error(simulate_cnma_data(n_studies = 0), "positive integer")
  expect_error(simulate_cnma_data(n_studies = -5), "positive integer")
})

test_that("simulated data passes validation", {
  data <- simulate_cnma_data(n_studies = 30, seed = 42)
  expect_silent(validate_cnma_input(data))
})

test_that("simulated data has expected network structure", {
  data <- simulate_cnma_data(n_studies = 40, seed = 42)

  # Should have multiple treatments
  treatments <- unique(c(data$treat1, data$treat2))
  expect_true(length(treatments) >= 3)

  # Should include Placebo
  expect_true("Placebo" %in% treatments)

  # Should have multiple studies
  studies <- unique(data$studlab)
  expect_equal(length(studies), 40)
})
