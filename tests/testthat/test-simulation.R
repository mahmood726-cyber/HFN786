test_that("simulate_cnma_data generates correct structure", {
  data <- simulate_cnma_data(n_studies = 10, seed = 123)

  # Check class
  expect_s3_class(data, "data.frame")

  # Check required columns
  required_cols <- c("studlab", "treat1", "treat2", "TE", "seTE")
  expect_true(all(required_cols %in% names(data)))

  # Check optional columns
  optional_cols <- c("year", "is_rct", "study_design", "grade", "rob2",
                    "age_mean", "female_pct", "bmi_mean", "charlson", "baseline_risk")
  expect_true(all(optional_cols %in% names(data)))
})

test_that("simulate_cnma_data generates valid values", {
  data <- simulate_cnma_data(n_studies = 15, seed = 42)

  # All TE should be finite
  expect_true(all(is.finite(data$TE)))

  # All seTE should be positive and finite
  expect_true(all(is.finite(data$seTE)))
  expect_true(all(data$seTE > 0))

  # Study labels should be unique per comparison
  expect_true(length(unique(data$studlab)) <= 15)

  # Treatments should be characters
  expect_type(data$treat1, "character")
  expect_type(data$treat2, "character")
})

test_that("simulate_cnma_data respects seed", {
  data1 <- simulate_cnma_data(n_studies = 10, seed = 999)
  data2 <- simulate_cnma_data(n_studies = 10, seed = 999)

  expect_equal(data1$TE, data2$TE)
  expect_equal(data1$seTE, data2$seTE)
})

test_that("simulate_cnma_data generates different data with different seeds", {
  data1 <- simulate_cnma_data(n_studies = 10, seed = 111)
  data2 <- simulate_cnma_data(n_studies = 10, seed = 222)

  expect_false(identical(data1$TE, data2$TE))
})

test_that("simulate_cnma_data generates expected treatments", {
  data <- simulate_cnma_data(n_studies = 20, seed = 42)

  treatments <- unique(c(data$treat1, data$treat2))
  expected_treatments <- c("Placebo", "DrugA", "DrugB", "DrugC", "DrugD")

  expect_true(all(treatments %in% expected_treatments))
})

test_that("simulate_cnma_data generates valid GRADE levels", {
  data <- simulate_cnma_data(n_studies = 20, seed = 42)

  expected_levels <- c("High", "Moderate", "Low", "Very low")
  expect_true(all(levels(data$grade) %in% expected_levels))
  expect_true(all(data$grade %in% expected_levels))
})

test_that("simulate_cnma_data generates valid RoB2 levels", {
  data <- simulate_cnma_data(n_studies = 20, seed = 42)

  expected_rob <- c("low", "some concerns", "high")
  expect_true(all(data$rob2 %in% expected_rob))
})

test_that("simulate_cnma_data generates valid is_rct values", {
  data <- simulate_cnma_data(n_studies = 20, seed = 42)

  expect_true(all(data$is_rct %in% c(0, 1)))
})

test_that("simulate_cnma_data generates reasonable covariate values", {
  data <- simulate_cnma_data(n_studies = 30, seed = 42)

  # Age should be reasonable (adult population)
  expect_true(all(data$age_mean >= 18))
  expect_true(all(data$age_mean <= 100))

  # Female percentage should be 0-1
  expect_true(all(data$female_pct >= 0))
  expect_true(all(data$female_pct <= 1))

  # BMI should be reasonable
  expect_true(all(data$bmi_mean >= 15))
  expect_true(all(data$bmi_mean <= 50))

  # Charlson should be non-negative
  expect_true(all(data$charlson >= 0))

  # Baseline risk should be 0-1
  expect_true(all(data$baseline_risk >= 0))
  expect_true(all(data$baseline_risk <= 1))
})

test_that("simulate_cnma_data generates at least some rows", {
  data <- simulate_cnma_data(n_studies = 5, seed = 42)

  expect_true(nrow(data) > 0)
  expect_true(nrow(data) >= 5)  # At least 5 comparisons
})

test_that("simulate_cnma_data handles small n_studies", {
  data <- simulate_cnma_data(n_studies = 3, seed = 42)

  expect_s3_class(data, "data.frame")
  expect_true(nrow(data) > 0)
})

test_that("simulate_cnma_data handles large n_studies", {
  data <- simulate_cnma_data(n_studies = 100, seed = 42)

  expect_s3_class(data, "data.frame")
  expect_true(nrow(data) > 0)
  expect_true(length(unique(data$studlab)) == 100)
})

test_that("simulate_cnma_data generates multi-arm studies", {
  data <- simulate_cnma_data(n_studies = 30, seed = 42)

  # Count comparisons per study
  comparisons_per_study <- table(data$studlab)

  # Some studies should have multiple comparisons (multi-arm)
  expect_true(any(comparisons_per_study > 1))
})

test_that("simulate_cnma_data creates valid study designs", {
  data <- simulate_cnma_data(n_studies = 20, seed = 42)

  expect_true(is.factor(data$study_design))
  expect_true(all(data$study_design %in% c("RCT", "Non-RCT")))
})

test_that("simulate_cnma_data year values are reasonable", {
  data <- simulate_cnma_data(n_studies = 20, seed = 42)

  expect_true(all(data$year >= 2010))
  expect_true(all(data$year <= 2025))
})
