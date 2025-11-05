# Tests for configuration functions

test_that("setup_cnma creates config object", {
  config <- setup_cnma()

  expect_s3_class(config, "cnma_config")
  expect_type(config, "list")
  expect_true("sm" %in% names(config))
  expect_equal(config$sm, "HR")
})

test_that("setup_cnma accepts custom parameters", {
  config <- setup_cnma(
    sm = "OR",
    use_bayesian = FALSE,
    seed = 123
  )

  expect_equal(config$sm, "OR")
  expect_false(config$use_bayesian)
  expect_equal(config$seed, 123)
})

test_that("setup_cnma validates input", {
  expect_error(setup_cnma(bayes_chains = 0), "bayes_chains must be positive")
  expect_error(setup_cnma(bayes_iter = 0), "bayes_iter must be positive")
  expect_error(setup_cnma(bayes_chains = -1), "bayes_chains must be positive")
  expect_error(setup_cnma(n_cores = 0), "n_cores must be positive")
})

test_that("setup_cnma validates bayes_warmup", {
  expect_error(
    setup_cnma(bayes_iter = 1000, bayes_warmup = 1500),
    "bayes_warmup must be less than bayes_iter"
  )
})

test_that("print.cnma_config works", {
  config <- setup_cnma()
  expect_output(print(config), "CNMA Configuration")
  expect_output(print(config), "Summary measure")
  expect_output(print(config), "Bayesian analysis")
})

test_that("config has all expected fields", {
  config <- setup_cnma()

  expected_fields <- c(
    "sm", "use_transport", "use_bayesian",
    "seed", "parallel_strategy", "n_cores"
  )

  for (field in expected_fields) {
    expect_true(field %in% names(config),
                info = sprintf("Missing field: %s", field))
  }
})
