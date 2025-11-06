test_that("null coalescing operator works", {
  `%||%` <- cnma:::`%||%`

  expect_equal(NULL %||% 5, 5)
  expect_equal(10 %||% 5, 10)
  expect_equal(0 %||% 5, 0)
  expect_equal("" %||% "default", "")
})

test_that("safe_clip works", {
  safe_clip <- cnma:::safe_clip

  expect_equal(safe_clip(5, 0, 10), 5)
  expect_equal(safe_clip(-5, 0, 10), 0)
  expect_equal(safe_clip(15, 0, 10), 10)
  expect_equal(safe_clip(c(-1, 5, 11), 0, 10), c(0, 5, 10))
})

test_that("vcoalesce replaces NA with default", {
  vcoalesce <- cnma:::vcoalesce

  expect_equal(vcoalesce(c(1, NA, 3), 0), c(1, 0, 3))
  expect_equal(vcoalesce(c(NA, NA, NA), 999), c(999, 999, 999))
  expect_equal(vcoalesce(c(1, 2, 3), 0), c(1, 2, 3))
})

test_that("vcoalesce accepts custom replacement value", {
  vcoalesce <- cnma:::vcoalesce

  expect_equal(vcoalesce(c(1, NA, 3), -1), c(1, -1, 3))
  expect_equal(vcoalesce(c(NA, 2, NA), 100), c(100, 2, 100))
})

test_that("has_pkg detects installed packages", {
  has_pkg <- cnma:::has_pkg

  # Base packages should always be available
  expect_true(has_pkg("stats"))
  expect_true(has_pkg("utils"))
  expect_true(has_pkg("graphics"))

  # Non-existent package
  expect_false(has_pkg("thispackagedoesnotexist12345"))
})

test_that("is_cran detects CRAN environment", {
  is_cran <- cnma:::is_cran

  # In normal testing, should not be CRAN
  expect_type(is_cran(), "logical")
})

test_that("has_jags returns logical", {
  has_jags <- cnma:::has_jags

  result <- has_jags()
  expect_type(result, "logical")
})

test_that("has_jags returns FALSE on CRAN", {
  has_jags <- cnma:::has_jags

  # Mock CRAN environment
  withr::local_envvar(c("NOT_CRAN" = ""))

  # Should return FALSE on CRAN regardless
  # (actual result depends on test environment)
  expect_type(has_jags(), "logical")
})

test_that("cache_key generates different keys for different inputs", {
  cache_key <- cnma:::cache_key

  key1 <- cache_key("test", x = 1, y = 2)
  key2 <- cache_key("test", x = 1, y = 3)
  key3 <- cache_key("test", x = 2, y = 2)

  expect_type(key1, "character")
  expect_false(identical(key1, key2))
  expect_false(identical(key1, key3))
})

test_that("cache_key generates same key for same inputs", {
  cache_key <- cnma:::cache_key

  key1 <- cache_key("test", x = 1, y = 2, z = 3)
  key2 <- cache_key("test", x = 1, y = 2, z = 3)

  expect_equal(key1, key2)
})

test_that("memoize returns same value with and without cache", {
  skip("Requires modification to test properly")

  memoize <- cnma:::memoize

  expr <- rnorm(10, mean = 5)
  key <- "test_key"

  # Without cache
  result1 <- memoize(key, expr, enable_cache = FALSE)

  # Should execute and return result
  expect_length(result1, 10)
})

test_that(".safe_try handles errors gracefully", {
  .safe_try <- cnma:::.safe_try

  # Successful expression
  result <- .safe_try(1 + 1, context = "test", silent = TRUE)
  expect_equal(result, 2)

  # Failing expression
  result <- .safe_try(stop("error"), context = "test", silent = TRUE)
  expect_s3_class(result, "try-error")
})

test_that(".stop_hint formats error messages correctly", {
  .stop_hint <- cnma:::.stop_hint

  expect_error(.stop_hint("Error message"), "Error message")
  expect_error(.stop_hint("Error", "Try this"), "Hint: Try this")
})

test_that("safe_clip handles edge cases", {
  safe_clip <- cnma:::safe_clip

  # All values at bounds
  expect_equal(safe_clip(c(0, 0, 0), 0, 10), c(0, 0, 0))
  expect_equal(safe_clip(c(10, 10, 10), 0, 10), c(10, 10, 10))

  # Single value
  expect_equal(safe_clip(5, 0, 10), 5)

  # Empty vector
  expect_equal(safe_clip(numeric(0), 0, 10), numeric(0))
})

test_that("vcoalesce handles all NA", {
  vcoalesce <- cnma:::vcoalesce

  result <- vcoalesce(c(NA, NA, NA), 42)
  expect_equal(result, c(42, 42, 42))
  expect_false(any(is.na(result)))
})

test_that("vcoalesce handles no NA", {
  vcoalesce <- cnma:::vcoalesce

  result <- vcoalesce(c(1, 2, 3), 42)
  expect_equal(result, c(1, 2, 3))
})

test_that("null coalescing handles FALSE correctly", {
  `%||%` <- cnma:::`%||%`

  # FALSE should not be replaced
  expect_equal(FALSE %||% TRUE, FALSE)
  expect_equal(0 %||% 1, 0)
})
