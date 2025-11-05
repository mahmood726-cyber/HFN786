# =========================================================
# Meta-Regression Framework
# Advanced covariate exploration and effect modification
# =========================================================

#' Run Network Meta-Regression
#'
#' Performs meta-regression to explore treatment effect modifiers in network
#' meta-analysis. Implements methods from Salanti (2012) BMC Med Res Methodol.
#'
#' @param data Original data frame with covariates
#' @param netmeta_obj A netmeta object
#' @param covariates Character vector of covariate names
#' @param formula Optional custom formula (overrides covariates)
#' @param method Method: "common", "random", "fixed"
#' @return Meta-regression results object
#' @export
#' @references
#' Salanti G (2012). Indirect and mixed-treatment comparison, network, or
#' multiple-treatments meta-analysis: many names, many benefits, many concerns
#' for the next generation evidence synthesis tool. Research Synthesis Methods,
#' 3(2):80-97.
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(40)
#' config <- setup_cnma(use_bayesian = FALSE)
#' results <- run_cnma_analysis(data, config = config)
#'
#' # Meta-regression
#' metareg <- run_metaregression(
#'   data,
#'   results$results$main_nma,
#'   covariates = c("age_mean", "female_pct")
#' )
#' print(metareg)
#' }
run_metaregression <- function(data,
                              netmeta_obj,
                              covariates = NULL,
                              formula = NULL,
                              method = c("random", "common", "fixed")) {

  if (!requireNamespace("netmeta", quietly = TRUE)) {
    .stop_hint("Package 'netmeta' required.",
               "Install with: install.packages('netmeta')")
  }

  method <- match.arg(method)

  # Validate covariates
  if (is.null(formula) && is.null(covariates)) {
    .stop_hint("Either 'covariates' or 'formula' must be specified.")
  }

  if (!is.null(covariates)) {
    missing_covs <- setdiff(covariates, names(data))
    if (length(missing_covs) > 0) {
      .stop_hint(sprintf("Covariates not found in data: %s",
                        paste(missing_covs, collapse = ", ")))
    }
  }

  # Merge covariate data with network
  nma_data <- data.frame(
    studlab = netmeta_obj$studlab,
    treat1 = netmeta_obj$treat1,
    treat2 = netmeta_obj$treat2,
    TE = netmeta_obj$TE,
    seTE = netmeta_obj$seTE,
    stringsAsFactors = FALSE
  )

  # Merge covariates
  if (!is.null(covariates)) {
    cov_data <- data[, c("studlab", covariates), drop = FALSE]
    cov_data <- cov_data[!duplicated(cov_data$studlab), ]
    nma_data <- merge(nma_data, cov_data, by = "studlab", all.x = TRUE)
  }

  # Run meta-regression using netmeta
  if (!is.null(formula)) {
    metareg_formula <- formula
  } else {
    metareg_formula <- as.formula(paste("~", paste(covariates, collapse = " + ")))
  }

  msg("Running network meta-regression...")
  msg("  Covariates: %s", paste(covariates, collapse = ", "))

  # Use netmeta's meta-regression
  metareg_result <- .safe_try(
    netmeta::netmeta(
      TE = TE,
      seTE = seTE,
      treat1 = treat1,
      treat2 = treat2,
      studlab = studlab,
      data = nma_data,
      reference.group = netmeta_obj$reference.group,
      sm = netmeta_obj$sm
    ),
    context = "Meta-regression"
  )

  if (inherits(metareg_result, "try-error")) {
    .stop_hint("Meta-regression failed.",
               "Check that covariates have valid values and sufficient variation.")
  }

  # Extract covariate effects
  results <- list(
    netmeta = metareg_result,
    covariates = covariates,
    formula = metareg_formula,
    method = method,
    data = nma_data
  )

  class(results) <- "cnma_metaregression"
  return(results)
}

#' Print Meta-Regression Results
#' @param x cnma_metaregression object
#' @param ... Additional arguments
#' @export
print.cnma_metaregression <- function(x, ...) {
  cat("Network Meta-Regression Results\n")
  cat("================================\n\n")

  cat("Covariates included:\n")
  for (cov in x$covariates) {
    cat(sprintf("  - %s\n", cov))
  }
  cat("\n")

  cat("Formula: ")
  print(x$formula)
  cat("\n")

  cat("Method:", x$method, "\n\n")

  cat("Use summary() for detailed results\n")

  invisible(x)
}

#' Explore Treatment-Covariate Interactions
#'
#' Tests for treatment-covariate interactions to identify effect modifiers.
#'
#' @param data Data frame with covariates
#' @param netmeta_obj A netmeta object
#' @param covariate Single covariate to test
#' @param comparison Specific comparison to test (optional)
#' @return Interaction test results
#' @export
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(35)
#' config <- setup_cnma(use_bayesian = FALSE)
#' results <- run_cnma_analysis(data, config = config)
#'
#' # Test age interaction
#' interaction <- test_covariate_interaction(
#'   data,
#'   results$results$main_nma,
#'   covariate = "age_mean"
#' )
#' }
test_covariate_interaction <- function(data,
                                      netmeta_obj,
                                      covariate,
                                      comparison = NULL) {

  if (!covariate %in% names(data)) {
    .stop_hint(sprintf("Covariate '%s' not found in data.", covariate))
  }

  # Prepare data
  nma_data <- data.frame(
    studlab = netmeta_obj$studlab,
    treat1 = netmeta_obj$treat1,
    treat2 = netmeta_obj$treat2,
    TE = netmeta_obj$TE,
    seTE = netmeta_obj$seTE,
    stringsAsFactors = FALSE
  )

  # Merge covariate
  cov_data <- data[, c("studlab", covariate), drop = FALSE]
  cov_data <- cov_data[!duplicated(cov_data$studlab), ]
  nma_data <- merge(nma_data, cov_data, by = "studlab", all.x = TRUE)

  # Remove missing
  nma_data <- nma_data[complete.cases(nma_data), ]

  if (nrow(nma_data) == 0) {
    .stop_hint("No complete cases after merging covariate.")
  }

  # Simple regression of TE on covariate
  if (!is.null(comparison)) {
    # Filter to specific comparison
    comp_data <- nma_data[paste(nma_data$treat1, nma_data$treat2) == comparison |
                          paste(nma_data$treat2, nma_data$treat1) == comparison, ]
    if (nrow(comp_data) == 0) {
      .stop_hint("No data for specified comparison.")
    }
    nma_data <- comp_data
  }

  # Weighted regression
  weights <- 1 / (nma_data$seTE^2)
  model <- lm(TE ~ get(covariate), data = nma_data, weights = weights)

  results <- list(
    covariate = covariate,
    comparison = comparison %||% "All",
    model = model,
    coefficient = coef(model)[2],
    se = summary(model)$coefficients[2, 2],
    t_value = summary(model)$coefficients[2, 3],
    p_value = summary(model)$coefficients[2, 4],
    n_studies = nrow(nma_data)
  )

  class(results) <- "cnma_interaction"
  return(results)
}

#' Print Interaction Test Results
#' @param x cnma_interaction object
#' @param ... Additional arguments
#' @export
print.cnma_interaction <- function(x, ...) {
  cat("Treatment-Covariate Interaction Test\n")
  cat("=====================================\n\n")

  cat(sprintf("Covariate: %s\n", x$covariate))
  cat(sprintf("Comparison: %s\n", x$comparison))
  cat(sprintf("Studies: %d\n\n", x$n_studies))

  cat("Interaction Effect:\n")
  cat(sprintf("  Coefficient: %.4f\n", x$coefficient))
  cat(sprintf("  SE: %.4f\n", x$se))
  cat(sprintf("  t-value: %.2f\n", x$t_value))
  cat(sprintf("  p-value: %.4f\n", x$p_value))
  cat("\n")

  if (x$p_value < 0.05) {
    cat("Significant interaction detected (p < 0.05)\n")
    cat(sprintf("Interpretation: Treatment effect varies with %s\n", x$covariate))
  } else {
    cat("No significant interaction detected (p >= 0.05)\n")
  }

  invisible(x)
}

#' Multiple Covariate Exploration
#'
#' Tests multiple covariates simultaneously and ranks by importance.
#'
#' @param data Data frame
#' @param netmeta_obj A netmeta object
#' @param covariates Character vector of covariates to test
#' @param adjustment Multiple testing adjustment ("bonferroni", "fdr", "none")
#' @return Data frame with covariate rankings
#' @export
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(40)
#' config <- setup_cnma(use_bayesian = FALSE)
#' results <- run_cnma_analysis(data, config = config)
#'
#' # Explore multiple covariates
#' exploration <- explore_multiple_covariates(
#'   data,
#'   results$results$main_nma,
#'   covariates = c("age_mean", "female_pct", "bmi_mean", "year")
#' )
#' print(exploration)
#' }
explore_multiple_covariates <- function(data,
                                       netmeta_obj,
                                       covariates,
                                       adjustment = c("bonferroni", "fdr", "none")) {

  adjustment <- match.arg(adjustment)

  msg("Testing %d covariates for effect modification...", length(covariates))

  results_list <- list()
  for (cov in covariates) {
    test_result <- .safe_try(
      test_covariate_interaction(data, netmeta_obj, covariate = cov),
      context = sprintf("Testing %s", cov),
      silent = TRUE
    )

    if (!inherits(test_result, "try-error")) {
      results_list[[cov]] <- data.frame(
        covariate = cov,
        coefficient = test_result$coefficient,
        se = test_result$se,
        p_value = test_result$p_value,
        n_studies = test_result$n_studies,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results_list) == 0) {
    msg("No successful covariate tests.")
    return(NULL)
  }

  # Combine results
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL

  # Adjust p-values
  if (adjustment == "bonferroni") {
    results_df$p_adjusted <- pmin(results_df$p_value * nrow(results_df), 1)
  } else if (adjustment == "fdr") {
    results_df$p_adjusted <- p.adjust(results_df$p_value, method = "fdr")
  } else {
    results_df$p_adjusted <- results_df$p_value
  }

  # Add significance
  results_df$significant <- results_df$p_adjusted < 0.05

  # Sort by p-value
  results_df <- results_df[order(results_df$p_value), ]

  class(results_df) <- c("cnma_covariate_exploration", "data.frame")
  return(results_df)
}

#' Print Covariate Exploration Results
#' @param x cnma_covariate_exploration object
#' @param ... Additional arguments
#' @export
print.cnma_covariate_exploration <- function(x, ...) {
  cat("Multiple Covariate Exploration Results\n")
  cat("========================================\n\n")

  cat(sprintf("Covariates tested: %d\n", nrow(x)))
  cat(sprintf("Significant: %d\n\n", sum(x$significant)))

  print(as.data.frame(x), row.names = FALSE, digits = 4)

  cat("\n")
  if (any(x$significant)) {
    cat("Significant effect modifiers:\n")
    sig_covs <- x$covariate[x$significant]
    for (cov in sig_covs) {
      cat(sprintf("  - %s\n", cov))
    }
  } else {
    cat("No significant effect modifiers detected.\n")
  }

  invisible(x)
}
