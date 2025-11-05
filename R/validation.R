# =========================================================
# Data Validation Functions
# =========================================================

#' Validate CNMA Input Data
#'
#' Validates that input data contains required columns and valid values.
#'
#' @param data Data frame with network meta-analysis data
#' @return Invisible TRUE if valid
#' @export
#' @examples
#' \dontrun{
#' data <- simulate_cnma_data(20)
#' validate_cnma_input(data)
#' }
validate_cnma_input <- function(data) {
  # Check if data is a data frame
  if (!is.data.frame(data)) {
    .stop_hint("Input must be a data frame.",
               "Convert your data to a data.frame using as.data.frame().")
  }

  # Check for empty data
  if (nrow(data) == 0) {
    .stop_hint("Input data has zero rows.",
               "Provide a dataset with at least one study comparison.")
  }

  # Check required columns
  need <- c("studlab", "treat1", "treat2", "TE", "seTE")
  miss <- setdiff(need, names(data))

  if (length(miss)) {
    .stop_hint(sprintf("Input data missing columns: %s",
                      paste(miss, collapse = ", ")),
              "If you have arm-level data, use make_pairwise_from_arms() first.")
  }

  # Check for non-finite values
  if (any(!is.finite(data$TE) | !is.finite(data$seTE))) {
    .stop_hint("TE/seTE contain non-finite values (NA, NaN, or Inf).",
               "Remove or impute missing values before analysis.")
  }

  # Check for non-positive standard errors
  if (any(data$seTE <= 0)) {
    n_invalid <- sum(data$seTE <= 0)
    .stop_hint(sprintf("seTE must be strictly positive. Found %d invalid value(s).", n_invalid),
               "Check your data for zero or negative standard errors.")
  }

  # Check for duplicate comparisons within studies
  study_treat_combos <- paste(data$studlab, data$treat1, data$treat2)
  if (any(duplicated(study_treat_combos))) {
    n_dups <- sum(duplicated(study_treat_combos))
    .stop_hint(sprintf("Found %d duplicate treatment comparisons within studies.", n_dups),
               "Each study-treatment pair should appear only once.")
  }

  invisible(TRUE)
}

#' Clean CNMA Data
#'
#' Cleans input data by ensuring proper types and removing invalid rows.
#'
#' @param data Input data frame
#' @return Cleaned data frame
#' @export
#' @examples
#' data <- data.frame(
#'   studlab = c("S1", "S2"),
#'   treat1 = c("A", "A"),
#'   treat2 = c("B", "C"),
#'   TE = c(0.5, 0.3),
#'   seTE = c(0.1, 0.2)
#' )
#' clean_data <- cnma_clean_data(data)
cnma_clean_data <- function(data) {
  # Store original row count
  n_orig <- nrow(data)

  # Clean data
  data_clean <- data %>%
    mutate(
      studlab = as.character(studlab),
      treat1 = as.character(treat1),
      treat2 = as.character(treat2)
    ) %>%
    filter(is.finite(TE), is.finite(seTE), seTE > 0)

  # Report if rows were removed
  n_removed <- n_orig - nrow(data_clean)
  if (n_removed > 0) {
    msg("Removed %d row(s) with invalid values", n_removed)
  }

  data_clean
}
