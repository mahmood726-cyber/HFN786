test_that("simulate_cnma_data works", {
  # Test basic simulation
  data <- simulate_cnma_data(n_studies = 10, seed = 123)

  expect_s3_class(data, "data.frame")
  expect_true(nrow(data) > 0)
  expect_true(all(c("studlab", "treat1", "treat2", "TE", "seTE") %in% names(data)))
  expect_true(all(is.finite(data$TE)))
  expect_true(all(is.finite(data$seTE)))
  expect_true(all(data$seTE > 0))
})

test_that("cnma_clean_data works", {
  # Test data cleaning
  data <- data.frame(
    studlab = c("S1", "S2", "S3"),
    treat1 = c("A", "A", "A"),
    treat2 = c("B", "C", "B"),
    TE = c(0.5, 0.3, NA),
    seTE = c(0.1, 0.2, 0.15)
  )

  cleaned <- cnma_clean_data(data)

  expect_equal(nrow(cleaned), 2)
  expect_true(all(is.finite(cleaned$TE)))
})

test_that("validate_cnma_input detects missing columns", {
  data <- data.frame(
    studlab = "S1",
    treat1 = "A"
  )

  expect_error(validate_cnma_input(data), "missing columns")
})

test_that("validate_cnma_input detects invalid seTE", {
  data <- data.frame(
    studlab = "S1",
    treat1 = "A",
    treat2 = "B",
    TE = 0.5,
    seTE = -0.1
  )

  expect_error(validate_cnma_input(data), "strictly positive")
})

test_that("setup_cnma creates valid configuration", {
  config <- setup_cnma(sm = "OR", use_bayesian = FALSE)

  expect_s3_class(config, "cnma_config")
  expect_equal(config$sm, "OR")
  expect_false(config$use_bayesian)
})

test_that("cnma_parallel_on/off work", {
  # Test parallel controls
  result_off <- cnma_parallel_off()
  expect_null(result_off)

  result_on <- cnma_parallel_on("sequential")
  expect_type(result_on, "logical")
})

test_that("run_cnma_analysis works with minimal data", {
  skip_on_cran()

  data <- simulate_cnma_data(n_studies = 10, seed = 42)
  config <- setup_cnma(
    use_bayesian = FALSE,
    report_html = FALSE,
    export_plots = FALSE,
    run_metareg = FALSE
  )

  result <- run_cnma_analysis(data, config = config)

  expect_s3_class(result, "cnma")
  expect_true("results" %in% names(result))
  expect_true("config" %in% names(result))
  expect_true("ref_treatment" %in% names(result))
})

test_that("print.cnma works", {
  skip_on_cran()

  data <- simulate_cnma_data(n_studies = 5, seed = 123)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
  result <- run_cnma_analysis(data, config = config)

  expect_output(print(result), "CNMA Results")
})

test_that("summary.cnma works", {
  skip_on_cran()

  data <- simulate_cnma_data(n_studies = 5, seed = 123)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
  result <- run_cnma_analysis(data, config = config)

  expect_output(summary(result), "CNMA Analysis Summary")
})
