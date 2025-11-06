test_that("setup_cnma creates config with default values", {
  config <- setup_cnma()

  expect_s3_class(config, "cnma_config")
  expect_equal(config$sm, "HR")
  expect_true(config$use_bayesian)
  expect_true(config$use_transport)
  expect_equal(config$seed, 42)
})

test_that("setup_cnma accepts custom summary measure", {
  measures <- c("OR", "RR", "HR", "MD", "SMD")

  for (sm in measures) {
    config <- setup_cnma(sm = sm)
    expect_equal(config$sm, sm)
  }
})

test_that("setup_cnma accepts custom transport metrics", {
  config_maha <- setup_cnma(transport_metric = "mahalanobis")
  expect_equal(config_maha$transport_metric, "mahalanobis")

  config_eucl <- setup_cnma(transport_metric = "euclidean")
  expect_equal(config_eucl$transport_metric, "euclidean")
})

test_that("setup_cnma accepts custom transport kernels", {
  config_gauss <- setup_cnma(transport_kernel = "gaussian")
  expect_equal(config_gauss$transport_kernel, "gaussian")

  config_tri <- setup_cnma(transport_kernel = "tricube")
  expect_equal(config_tri$transport_kernel, "tricube")
})

test_that("setup_cnma accepts custom Bayesian settings", {
  config <- setup_cnma(
    use_bayesian = TRUE,
    bayes_chains = 4,
    bayes_iter = 5000,
    bayes_warmup = 2500
  )

  expect_true(config$use_bayesian)
  expect_equal(config$bayes_chains, 4)
  expect_equal(config$bayes_iter, 5000)
  expect_equal(config$bayes_warmup, 2500)
})

test_that("setup_cnma accepts custom parallel settings", {
  config <- setup_cnma(
    parallel_strategy = "multisession",
    n_cores = 2
  )

  expect_equal(config$parallel_strategy, "multisession")
  expect_equal(config$n_cores, 2)
})

test_that("setup_cnma accepts custom output settings", {
  config <- setup_cnma(
    export_plots = TRUE,
    plot_dir = "my_plots",
    output_dir = "my_output",
    report_html = TRUE,
    report_file = "my_report.html"
  )

  expect_true(config$export_plots)
  expect_equal(config$plot_dir, "my_plots")
  expect_equal(config$output_dir, "my_output")
  expect_true(config$report_html)
  expect_equal(config$report_file, "my_report.html")
})

test_that("print.cnma_config displays config information", {
  config <- setup_cnma(sm = "OR", use_bayesian = FALSE)

  expect_output(print(config), "CNMA Configuration")
  expect_output(print(config), "Summary measure: OR")
  expect_output(print(config), "Bayesian analysis: FALSE")
})

test_that("setup_cnma accepts weighting options", {
  config <- setup_cnma(
    use_grade_weighting = TRUE,
    use_design_weighting = TRUE,
    use_rob2_weighting = TRUE
  )

  expect_true(config$use_grade_weighting)
  expect_true(config$use_design_weighting)
  expect_true(config$use_rob2_weighting)
})

test_that("setup_cnma accepts meta-regression options", {
  config <- setup_cnma(
    run_metareg = TRUE,
    metareg_covariates = c("age", "sex"),
    metareg_spline_covars = c("age"),
    metareg_spline_df = 4,
    metareg_cr2 = TRUE
  )

  expect_true(config$run_metareg)
  expect_equal(config$metareg_covariates, c("age", "sex"))
  expect_equal(config$metareg_spline_covars, c("age"))
  expect_equal(config$metareg_spline_df, 4)
  expect_true(config$metareg_cr2)
})

test_that("setup_cnma accepts diagnostic options", {
  config <- setup_cnma(
    run_loo = TRUE,
    run_loto = TRUE,
    enable_ume = TRUE,
    enable_egger = TRUE,
    run_netheat = TRUE
  )

  expect_true(config$run_loo)
  expect_true(config$run_loto)
  expect_true(config$enable_ume)
  expect_true(config$enable_egger)
  expect_true(config$run_netheat)
})

test_that("setup_cnma accepts publication bias options", {
  config <- setup_cnma(
    run_selection_models = TRUE,
    run_copas = TRUE,
    run_trimfill = TRUE,
    pet_peese_min_k = 15
  )

  expect_true(config$run_selection_models)
  expect_true(config$run_copas)
  expect_true(config$run_trimfill)
  expect_equal(config$pet_peese_min_k, 15)
})

test_that("setup_cnma sets glmm_outcome correctly", {
  config_bin <- setup_cnma(glmm_outcome = "binary")
  expect_equal(config_bin$glmm_outcome, "binary")

  config_cont <- setup_cnma(glmm_outcome = "continuous")
  expect_equal(config_cont$glmm_outcome, "continuous")
})

test_that("setup_cnma preserves all settings", {
  config <- setup_cnma(
    sm = "RR",
    use_transport = FALSE,
    use_bayesian = FALSE,
    seed = 999
  )

  expect_equal(config$sm, "RR")
  expect_false(config$use_transport)
  expect_false(config$use_bayesian)
  expect_equal(config$seed, 999)
})
