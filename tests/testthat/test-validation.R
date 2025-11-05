# Tests for validation functions

test_that("validate_cnma_input accepts valid data", {
  data <- data.frame(
    studlab = c("S1", "S2"),
    treat1 = c("A", "A"),
    treat2 = c("B", "C"),
    TE = c(0.5, 0.3),
    seTE = c(0.1, 0.2)
  )

  expect_silent(validate_cnma_input(data))
  expect_true(validate_cnma_input(data))
})

test_that("validate_cnma_input rejects missing columns", {
  data <- data.frame(
    studlab = c("S1", "S2"),
    treat1 = c("A", "A"),
    treat2 = c("B", "C")
  )

  expect_error(validate_cnma_input(data), "missing columns")
})

test_that("validate_cnma_input rejects non-finite values", {
  data <- data.frame(
    studlab = c("S1", "S2"),
    treat1 = c("A", "A"),
    treat2 = c("B", "C"),
    TE = c(NA, 0.3),
    seTE = c(0.1, 0.2)
  )

  expect_error(validate_cnma_input(data), "non-finite")
})

test_that("validate_cnma_input rejects non-positive seTE", {
  data <- data.frame(
    studlab = c("S1", "S2"),
    treat1 = c("A", "A"),
    treat2 = c("B", "C"),
    TE = c(0.5, 0.3),
    seTE = c(0, 0.2)
  )

  expect_error(validate_cnma_input(data), "strictly positive")
})

test_that("validate_cnma_input rejects empty data", {
  data <- data.frame(
    studlab = character(0),
    treat1 = character(0),
    treat2 = character(0),
    TE = numeric(0),
    seTE = numeric(0)
  )

  expect_error(validate_cnma_input(data), "zero rows")
})

test_that("validate_cnma_input rejects non-data.frame input", {
  expect_error(validate_cnma_input(list(a = 1)), "data frame")
  expect_error(validate_cnma_input(matrix(1:10, 5, 2)), "data frame")
})

test_that("cnma_clean_data removes invalid rows", {
  data <- data.frame(
    studlab = c("S1", "S2", "S3", "S4"),
    treat1 = c("A", "A", "A", "A"),
    treat2 = c("B", "C", "D", "E"),
    TE = c(0.5, NA, Inf, 0.3),
    seTE = c(0.1, 0.2, 0.15, 0)
  )

  clean <- cnma_clean_data(data)

  expect_equal(nrow(clean), 1)
  expect_equal(clean$studlab, "S1")
})

test_that("cnma_clean_data converts types", {
  data <- data.frame(
    studlab = factor(c("S1", "S2")),
    treat1 = factor(c("A", "A")),
    treat2 = factor(c("B", "C")),
    TE = c(0.5, 0.3),
    seTE = c(0.1, 0.2)
  )

  clean <- cnma_clean_data(data)

  expect_type(clean$studlab, "character")
  expect_type(clean$treat1, "character")
  expect_type(clean$treat2, "character")
})
