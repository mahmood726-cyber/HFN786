# =========================================================
# Inconsistency Assessment Functions
# Aligned with PRISMA-NMA and Statistics in Medicine guidelines
# =========================================================

#' Assess Network Inconsistency
#'
#' Performs comprehensive inconsistency assessment using multiple methods
#' recommended by PRISMA-NMA guidelines and recent methodological papers.
#'
#' @param netmeta_obj A netmeta object from netmeta::netmeta()
#' @param methods Character vector of methods to use. Options:
#'   "global" (Cochran's Q), "local" (netsplit), "design" (design-by-treatment)
#' @return List with inconsistency assessment results
#' @export
#' @references
#' Dias S, Welton NJ, Caldwell DM, Ades AE (2010). Checking consistency in
#' mixed treatment comparison meta-analysis. Statistics in Medicine, 29(7-8):932-44.
#'
#' Krahn U, Binder H, König J (2013). A graphical tool for locating
#' inconsistency in network meta-analyses. BMC Medical Research Methodology, 13:35.
#'
#' Higgins JP, Jackson D, Barrett JK, Lu G, Ades AE, White IR (2012).
#' Consistency and inconsistency in network meta-analysis: concepts and models
#' for multi-arm studies. Research Synthesis Methods, 3(2):98-110.
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(40)
#' config <- setup_cnma(use_bayesian = FALSE)
#' results <- run_cnma_analysis(data, config = config)
#' inconsistency <- assess_inconsistency(results$results$main_nma)
#' print(inconsistency)
#' }
assess_inconsistency <- function(netmeta_obj,
                                 methods = c("global", "local", "design")) {
  if (!requireNamespace("netmeta", quietly = TRUE)) {
    .stop_hint("Package 'netmeta' is required.",
               "Install it with: install.packages('netmeta')")
  }

  if (inherits(netmeta_obj, "try-error") || is.null(netmeta_obj)) {
    .stop_hint("Invalid netmeta object provided.")
  }

  results <- list()

  # Global inconsistency (Cochran's Q)
  if ("global" %in% methods) {
    results$global <- list(
      Q = netmeta_obj$Q,
      df = netmeta_obj$df.Q,
      p_value = netmeta_obj$pval.Q,
      I2 = netmeta_obj$I2.random,
      tau2 = netmeta_obj$tau^2,
      interpretation = ifelse(
        netmeta_obj$pval.Q < 0.05,
        "Significant heterogeneity/inconsistency detected",
        "No significant heterogeneity/inconsistency"
      )
    )
  }

  # Local inconsistency (node-splitting / netsplit)
  if ("local" %in% methods) {
    split_result <- .safe_try(
      netmeta::netsplit(netmeta_obj),
      context = "Node-splitting analysis"
    )

    if (!inherits(split_result, "try-error")) {
      # Extract comparisons with significant inconsistency
      p_values <- split_result$compare.random$p
      inconsistent <- p_values < 0.05

      results$local <- list(
        split_object = split_result,
        n_comparisons = length(p_values),
        n_inconsistent = sum(inconsistent, na.rm = TRUE),
        inconsistent_comparisons = if (any(inconsistent, na.rm = TRUE)) {
          names(p_values)[inconsistent]
        } else {
          "None"
        },
        interpretation = ifelse(
          sum(inconsistent, na.rm = TRUE) == 0,
          "No significant local inconsistency detected",
          sprintf("%d comparison(s) show significant inconsistency",
                  sum(inconsistent, na.rm = TRUE))
        )
      )
    } else {
      results$local <- list(
        error = "Node-splitting analysis failed",
        reason = "Insufficient data or network structure"
      )
    }
  }

  # Design-by-treatment inconsistency
  if ("design" %in% methods) {
    design_result <- .safe_try(
      netmeta::decomp.design(netmeta_obj),
      context = "Design-by-treatment decomposition"
    )

    if (!inherits(design_result, "try-error")) {
      results$design <- list(
        Q_inconsistency = design_result$Q.inconsistency,
        df_inconsistency = design_result$df.Q.inconsistency,
        p_inconsistency = design_result$pval.Q.inconsistency,
        Q_heterogeneity = design_result$Q.heterogeneity,
        interpretation = ifelse(
          design_result$pval.Q.inconsistency < 0.05,
          "Significant design-by-treatment inconsistency",
          "No significant design-by-treatment inconsistency"
        )
      )
    } else {
      results$design <- list(
        error = "Design decomposition failed",
        reason = "Insufficient data or network structure"
      )
    }
  }

  class(results) <- "cnma_inconsistency"
  return(results)
}

#' Print Inconsistency Assessment
#' @param x cnma_inconsistency object
#' @param ... Additional arguments
#' @export
print.cnma_inconsistency <- function(x, ...) {
  cat("Network Inconsistency Assessment\n")
  cat("==================================\n\n")

  # Global
  if (!is.null(x$global)) {
    cat("1. GLOBAL HETEROGENEITY/INCONSISTENCY\n")
    cat("   (Cochran's Q statistic)\n")
    cat("   ------------------------------------\n")
    cat(sprintf("   Q = %.2f (df = %d, p = %.4f)\n",
                x$global$Q, x$global$df, x$global$p_value))
    cat(sprintf("   I² = %.1f%%\n", x$global$I2 * 100))
    cat(sprintf("   Tau² = %.4f\n", x$global$tau2))
    cat(sprintf("   Interpretation: %s\n\n", x$global$interpretation))
  }

  # Local
  if (!is.null(x$local)) {
    cat("2. LOCAL INCONSISTENCY\n")
    cat("   (Node-splitting approach)\n")
    cat("   ------------------------------------\n")
    if (!is.null(x$local$error)) {
      cat(sprintf("   Error: %s\n", x$local$error))
      cat(sprintf("   Reason: %s\n\n", x$local$reason))
    } else {
      cat(sprintf("   Comparisons evaluated: %d\n", x$local$n_comparisons))
      cat(sprintf("   Inconsistent comparisons: %d\n", x$local$n_inconsistent))
      if (x$local$n_inconsistent > 0) {
        cat("   Comparisons with inconsistency:\n")
        for (comp in x$local$inconsistent_comparisons) {
          cat(sprintf("     - %s\n", comp))
        }
      }
      cat(sprintf("   Interpretation: %s\n\n", x$local$interpretation))
    }
  }

  # Design
  if (!is.null(x$design)) {
    cat("3. DESIGN-BY-TREATMENT INCONSISTENCY\n")
    cat("   ------------------------------------\n")
    if (!is.null(x$design$error)) {
      cat(sprintf("   Error: %s\n", x$design$error))
      cat(sprintf("   Reason: %s\n\n", x$design$reason))
    } else {
      cat(sprintf("   Q (inconsistency) = %.2f (df = %d, p = %.4f)\n",
                  x$design$Q_inconsistency,
                  x$design$df_inconsistency,
                  x$design$p_inconsistency))
      cat(sprintf("   Q (heterogeneity) = %.2f\n", x$design$Q_heterogeneity))
      cat(sprintf("   Interpretation: %s\n\n", x$design$interpretation))
    }
  }

  cat("Recommendations:\n")
  cat("  - If inconsistency is detected, investigate potential causes\n")
  cat("  - Consider effect modifiers, study quality, or network structure\n")
  cat("  - Sensitivity analyses may help identify sources of inconsistency\n")

  invisible(x)
}

#' Calculate Contribution Matrix
#'
#' Computes the contribution of each direct comparison to the network
#' meta-analysis estimates. Based on methods from Papakonstantinou et al. (2018).
#'
#' @param netmeta_obj A netmeta object
#' @return Matrix showing contribution percentages
#' @export
#' @references
#' Papakonstantinou T, Nikolakopoulou A, Rücker G, et al. (2018).
#' Estimating the contribution of studies in network meta-analysis:
#' paths, flows and streams. F1000Research, 7:610.
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(30)
#' config <- setup_cnma(use_bayesian = FALSE)
#' results <- run_cnma_analysis(data, config = config)
#' contrib <- calculate_contribution_matrix(results$results$main_nma)
#' print(contrib)
#' }
calculate_contribution_matrix <- function(netmeta_obj) {
  if (!requireNamespace("netmeta", quietly = TRUE)) {
    .stop_hint("Package 'netmeta' is required.",
               "Install it with: install.packages('netmeta')")
  }

  if (inherits(netmeta_obj, "try-error") || is.null(netmeta_obj)) {
    .stop_hint("Invalid netmeta object provided.")
  }

  # netmeta 2.0+ has netcontrib function
  contrib_result <- .safe_try(
    netmeta::netcontrib(netmeta_obj),
    context = "Contribution matrix calculation"
  )

  if (inherits(contrib_result, "try-error")) {
    msg("Contribution matrix calculation not available.")
    msg("This requires netmeta version >= 2.0.0 and sufficient network structure.")
    return(NULL)
  }

  class(contrib_result) <- c("cnma_contribution", class(contrib_result))
  return(contrib_result)
}

#' Assess Transitivity Assumption
#'
#' Evaluates the transitivity assumption (similarity of studies across comparisons)
#' which is fundamental to network meta-analysis validity.
#'
#' @param data Original data frame with study characteristics
#' @param variables Character vector of variables to assess
#' @param by_comparison Logical; assess by treatment comparison
#' @return Data frame with transitivity assessment
#' @export
#' @references
#' Jansen JP, Naci H (2013). Is network meta-analysis as valid as standard
#' pairwise meta-analysis? It all depends on the distribution of effect modifiers.
#' BMC Medicine, 11:159.
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(40)
#' transitivity <- assess_transitivity(
#'   data,
#'   variables = c("age_mean", "female_pct", "year")
#' )
#' print(transitivity)
#' }
assess_transitivity <- function(data, variables = NULL, by_comparison = TRUE) {
  if (!is.data.frame(data)) {
    .stop_hint("data must be a data.frame")
  }

  # Auto-detect numeric variables if not specified
  if (is.null(variables)) {
    numeric_vars <- names(data)[sapply(data, is.numeric)]
    # Exclude TE and seTE
    variables <- setdiff(numeric_vars, c("TE", "seTE"))

    if (length(variables) == 0) {
      msg("No numeric covariates found for transitivity assessment.")
      return(NULL)
    }

    msg("Auto-detected variables: %s", paste(variables, collapse = ", "))
  }

  # Check variables exist
  missing_vars <- setdiff(variables, names(data))
  if (length(missing_vars) > 0) {
    .stop_hint(sprintf("Variables not found in data: %s",
                       paste(missing_vars, collapse = ", ")))
  }

  results <- list()

  if (by_comparison) {
    # Create comparison identifier
    data$comparison <- paste(data$treat1, "vs", data$treat2)

    # Summary statistics by comparison
    for (var in variables) {
      var_data <- data[, c("comparison", var)]
      var_data <- var_data[!is.na(var_data[[var]]), ]

      summary_by_comp <- var_data %>%
        group_by(comparison) %>%
        summarise(
          mean = mean(.data[[var]], na.rm = TRUE),
          sd = sd(.data[[var]], na.rm = TRUE),
          min = min(.data[[var]], na.rm = TRUE),
          max = max(.data[[var]], na.rm = TRUE),
          n = n(),
          .groups = "drop"
        )

      # Overall statistics
      overall <- data.frame(
        comparison = "Overall",
        mean = mean(var_data[[var]], na.rm = TRUE),
        sd = sd(var_data[[var]], na.rm = TRUE),
        min = min(var_data[[var]], na.rm = TRUE),
        max = max(var_data[[var]], na.rm = TRUE),
        n = nrow(var_data)
      )

      results[[var]] <- rbind(summary_by_comp, overall)
    }
  } else {
    # Overall summary only
    for (var in variables) {
      results[[var]] <- data.frame(
        variable = var,
        mean = mean(data[[var]], na.rm = TRUE),
        sd = sd(data[[var]], na.rm = TRUE),
        min = min(data[[var]], na.rm = TRUE),
        max = max(data[[var]], na.rm = TRUE),
        n = sum(!is.na(data[[var]]))
      )
    }
  }

  class(results) <- "cnma_transitivity"
  attr(results, "by_comparison") <- by_comparison
  return(results)
}

#' Print Transitivity Assessment
#' @param x cnma_transitivity object
#' @param ... Additional arguments
#' @export
print.cnma_transitivity <- function(x, ...) {
  cat("Transitivity Assessment\n")
  cat("========================\n\n")
  cat("Evaluating similarity of study characteristics across comparisons.\n")
  cat("Large differences may indicate violation of transitivity assumption.\n\n")

  by_comp <- attr(x, "by_comparison")

  if (by_comp) {
    for (var_name in names(x)) {
      cat(sprintf("\nVariable: %s\n", var_name))
      cat("-----------------------------------\n")
      print(x[[var_name]], row.names = FALSE, digits = 2)
    }
  } else {
    for (var_name in names(x)) {
      print(x[[var_name]], row.names = FALSE, digits = 2)
    }
  }

  cat("\nInterpretation:\n")
  cat("  - Compare mean values across comparisons\n")
  cat("  - Large differences suggest potential effect modification\n")
  cat("  - Consider meta-regression if substantial imbalance exists\n")

  invisible(x)
}
