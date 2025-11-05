# Tests for utility functions

test_that("null coalescing operator works", {
  expect_equal(NULL %||% 5, 5)
  expect_equal(10 %||% 5, 10)
  expect_equal(0 %||% 5, 0)
  expect_equal(FALSE %||% TRUE, FALSE)
})

test_that("safe_clip works correctly", {
  expect_equal(safe_clip(5, 0, 10), 5)
  expect_equal(safe_clip(-5, 0, 10), 0)
  expect_equal(safe_clip(15, 0, 10), 10)
  expect_equal(safe_clip(c(-1, 5, 11), 0, 10), c(0, 5, 10))
})

test_that("vcoalesce replaces NA values", {
  x <- c(1, NA, 3, NA, 5)
  expect_equal(vcoalesce(x, 0), c(1, 0, 3, 0, 5))
  expect_equal(vcoalesce(x, -999), c(1, -999, 3, -999, 5))
  expect_equal(vcoalesce(c(1, 2, 3), 0), c(1, 2, 3))
})

test_that("has_pkg detects packages", {
  expect_true(has_pkg("stats"))
  expect_true(has_pkg("utils"))
  expect_false(has_pkg("nonexistent_package_xyz"))
})

test_that("is_cran works", {
  expect_type(is_cran(), "logical")
  expect_length(is_cran(), 1)
})

test_that("cache_key generates keys", {
  key1 <- cache_key("test", 1, 2, 3)
  key2 <- cache_key("test", 1, 2, 3)
  key3 <- cache_key("test", 1, 2, 4)

  expect_type(key1, "character")

  # Same inputs should generate same key if digest is available
  if (requireNamespace("digest", quietly = TRUE)) {
    expect_equal(key1, key2)
    expect_false(key1 == key3)
  }
})
