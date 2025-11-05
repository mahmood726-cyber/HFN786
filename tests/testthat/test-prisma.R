# Tests for PRISMA-NMA compliance functions

test_that("network_characteristics_summary generates summary", {
  data <- simulate_cnma_data(30, seed = 42)

  summary <- network_characteristics_summary(data)

  expect_s3_class(summary, "cnma_network_summary")
  expect_true("characteristic" %in% names(summary))
  expect_true("value" %in% names(summary))
  expect_true(nrow(summary) > 0)
})

test_that("network_characteristics_summary includes heterogeneity with netmeta object", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(25, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
  results <- run_cnma_analysis(data, config = config)

  summary <- network_characteristics_summary(data, results$results$main_nma)

  expect_true(any(grepl("Tau", summary$characteristic)))
  expect_true(any(grepl("I²", summary$characteristic)))
})

test_that("generate_prisma_report creates report", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(25, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
  results <- run_cnma_analysis(data, config = config)

  report <- generate_prisma_report(results, output_format = "text")

  expect_type(report, "character")
  expect_true(grepl("PRISMA-NMA", report))
  expect_true(grepl("Network Characteristics", report))
})

test_that("generate_prisma_report supports multiple formats", {
  skip_if_not_installed("netmeta")

  data <- simulate_cnma_data(20, seed = 42)
  config <- setup_cnma(use_bayesian = FALSE, export_plots = FALSE)
  results <- run_cnma_analysis(data, config = config)

  text_report <- generate_prisma_report(results, output_format = "text")
  md_report <- generate_prisma_report(results, output_format = "markdown")
  html_report <- generate_prisma_report(results, output_format = "html")

  expect_type(text_report, "character")
  expect_type(md_report, "character")
  expect_type(html_report, "character")

  expect_true(grepl("##", md_report))  # Markdown headers
  expect_true(grepl("<h", html_report))  # HTML tags
})

test_that("print method works for network summary", {
  data <- simulate_cnma_data(25, seed = 42)
  summary <- network_characteristics_summary(data)

  expect_output(print(summary), "Network Characteristics")
})
