test_that("validate_cnma_input detects non-finite TE values", {
  data <- data.frame(
    studlab = c("S1", "S2"),
    treat1 = c("A", "A"),
    treat2 = c("B", "C"),
    TE = c(0.5, Inf),
    seTE = c(0.1, 0.2)
  )

  expect_error(validate_cnma_input(data), "non-finite")
})

test_that("validate_cnma_input detects non-finite seTE values", {
  data <- data.frame(
    studlab = c("S1", "S2"),
    treat1 = c("A", "A"),
    treat2 = c("B", "C"),
    TE = c(0.5, 0.3),
    seTE = c(0.1, NaN)
  )

  expect_error(validate_cnma_input(data), "non-finite")
})

test_that("validate_cnma_input detects zero seTE values", {
  data <- data.frame(
    studlab = c("S1", "S2"),
    treat1 = c("A", "A"),
    treat2 = c("B", "C"),
    TE = c(0.5, 0.3),
    seTE = c(0.1, 0)
  )

  expect_error(validate_cnma_input(data), "strictly positive")
})

test_that("validate_cnma_input detects negative seTE values", {
  data <- data.frame(
    studlab = "S1",
    treat1 = "A",
    treat2 = "B",
    TE = 0.5,
    seTE = -0.1
  )

  expect_error(validate_cnma_input(data), "strictly positive")
})

test_that("validate_cnma_input passes valid data", {
  data <- data.frame(
    studlab = c("S1", "S2"),
    treat1 = c("A", "A"),
    treat2 = c("B", "C"),
    TE = c(0.5, 0.3),
    seTE = c(0.1, 0.2)
  )

  expect_invisible(validate_cnma_input(data))
})

test_that("cnma_clean_data removes rows with NA in TE", {
  data <- data.frame(
    studlab = c("S1", "S2", "S3"),
    treat1 = c("A", "A", "A"),
    treat2 = c("B", "C", "D"),
    TE = c(0.5, NA, 0.3),
    seTE = c(0.1, 0.2, 0.15)
  )

  cleaned <- cnma_clean_data(data)
  expect_equal(nrow(cleaned), 2)
  expect_false(any(is.na(cleaned$TE)))
})

test_that("cnma_clean_data removes rows with NA in seTE", {
  data <- data.frame(
    studlab = c("S1", "S2", "S3"),
    treat1 = c("A", "A", "A"),
    treat2 = c("B", "C", "D"),
    TE = c(0.5, 0.3, 0.4),
    seTE = c(0.1, NA, 0.15)
  )

  cleaned <- cnma_clean_data(data)
  expect_equal(nrow(cleaned), 2)
  expect_false(any(is.na(cleaned$seTE)))
})

test_that("cnma_clean_data removes rows with zero or negative seTE", {
  data <- data.frame(
    studlab = c("S1", "S2", "S3", "S4"),
    treat1 = c("A", "A", "A", "A"),
    treat2 = c("B", "C", "D", "E"),
    TE = c(0.5, 0.3, 0.4, 0.2),
    seTE = c(0.1, 0, -0.1, 0.2)
  )

  cleaned <- cnma_clean_data(data)
  expect_equal(nrow(cleaned), 2)
  expect_true(all(cleaned$seTE > 0))
})

test_that("cnma_clean_data converts character columns correctly", {
  data <- data.frame(
    studlab = factor(c("S1", "S2")),
    treat1 = factor(c("A", "A")),
    treat2 = factor(c("B", "C")),
    TE = c(0.5, 0.3),
    seTE = c(0.1, 0.2)
  )

  cleaned <- cnma_clean_data(data)
  expect_type(cleaned$studlab, "character")
  expect_type(cleaned$treat1, "character")
  expect_type(cleaned$treat2, "character")
})

test_that("cnma_clean_data preserves valid rows", {
  data <- data.frame(
    studlab = c("S1", "S2"),
    treat1 = c("A", "A"),
    treat2 = c("B", "C"),
    TE = c(0.5, 0.3),
    seTE = c(0.1, 0.2)
  )

  cleaned <- cnma_clean_data(data)
  expect_equal(nrow(cleaned), 2)
  expect_equal(cleaned$TE, c(0.5, 0.3))
})

test_that("validate_cnma_input detects multiple missing columns", {
  data <- data.frame(
    studlab = "S1",
    treat1 = "A"
  )

  expect_error(validate_cnma_input(data), "missing columns")
  expect_error(validate_cnma_input(data), "treat2")
})

test_that("cnma_clean_data handles Inf values", {
  data <- data.frame(
    studlab = c("S1", "S2", "S3"),
    treat1 = c("A", "A", "A"),
    treat2 = c("B", "C", "D"),
    TE = c(0.5, Inf, 0.3),
    seTE = c(0.1, 0.2, 0.15)
  )

  cleaned <- cnma_clean_data(data)
  expect_equal(nrow(cleaned), 2)
  expect_true(all(is.finite(cleaned$TE)))
})
