# Tests for parallel processing

test_that("cnma_parallel_on accepts valid strategies", {
  expect_silent(cnma_parallel_on("sequential"))
})

test_that("cnma_parallel_off works", {
  cnma_parallel_on("sequential")
  expect_silent(cnma_parallel_off())
})

test_that("parallel functions don't error without future", {
  expect_silent(cnma_parallel_on("sequential"))
  expect_silent(cnma_parallel_off())
})

test_that("cnma_parallel_on validates strategy", {
  expect_error(cnma_parallel_on("invalid_strategy"))
})

test_that("cnma_parallel_on returns logical", {
  result <- cnma_parallel_on("sequential")
  expect_type(result, "logical")
})
