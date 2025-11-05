# =========================================================
# Treatment Ranking Functions
# Based on best practices from Statistics in Medicine
# and Research Synthesis Methods
# =========================================================

#' Calculate Treatment Rankings with P-scores and SUCRA
#'
#' Computes treatment rankings using P-scores (frequentist analogue of SUCRA)
#' as described in Rücker & Schwarzer (2015) Statistics in Medicine.
#'
#' @param netmeta_obj A netmeta object from netmeta::netmeta()
#' @param small_values Character indicating if small values are beneficial
#'   ("desirable") or harmful ("undesirable"). Default: "undesirable"
#' @return Data frame with treatment rankings, P-scores, and probabilities
#' @export
#' @references
#' Rücker G, Schwarzer G (2015). Ranking treatments in frequentist network
#' meta-analysis works without resampling methods. BMC Medical Research
#' Methodology, 15:58.
#'
#' Salanti G, Ades AE, Ioannidis JP (2011). Graphical methods and numerical
#' summaries for presenting results from multiple-treatment meta-analysis:
#' an overview and tutorial. Journal of Clinical Epidemiology, 64(2):163-71.
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(30)
#' config <- setup_cnma(use_bayesian = FALSE)
#' results <- run_cnma_analysis(data, config = config)
#' rankings <- calculate_rankings(results$results$main_nma)
#' print(rankings)
#' }
calculate_rankings <- function(netmeta_obj, small_values = c("undesirable", "desirable")) {
  if (!requireNamespace("netmeta", quietly = TRUE)) {
    .stop_hint("Package 'netmeta' is required for ranking calculations.",
               "Install it with: install.packages('netmeta')")
  }

  if (inherits(netmeta_obj, "try-error") || is.null(netmeta_obj)) {
    .stop_hint("Invalid netmeta object provided.",
               "Ensure your analysis completed successfully.")
  }

  small_values <- match.arg(small_values)

  # Calculate P-scores using netmeta's netrank function
  rank_results <- netmeta::netrank(netmeta_obj, small.values = small_values)

  # Extract P-scores (analogous to SUCRA)
  pscores <- rank_results$Pscore.random
  treatments <- names(pscores)

  # Calculate ranks
  ranks <- rank(-pscores)

  # Create comprehensive ranking table
  ranking_table <- data.frame(
    treatment = treatments,
    p_score = round(pscores, 4),
    rank = ranks,
    interpretation = ifelse(
      pscores >= 0.75, "Very likely best",
      ifelse(pscores >= 0.50, "Likely effective",
             ifelse(pscores >= 0.25, "Possibly effective", "Likely inferior"))
    ),
    stringsAsFactors = FALSE
  )

  # Sort by P-score descending
  ranking_table <- ranking_table[order(-ranking_table$p_score), ]
  rownames(ranking_table) <- NULL

  # Add class for S3 methods
  class(ranking_table) <- c("cnma_ranking", "data.frame")

  return(ranking_table)
}

#' Print Treatment Rankings
#' @param x cnma_ranking object
#' @param ... Additional arguments
#' @export
print.cnma_ranking <- function(x, ...) {
  cat("Treatment Rankings (P-scores)\n")
  cat("==============================\n\n")
  cat("P-score: Probability that treatment is best (0-1 scale)\n")
  cat("Higher P-scores indicate better treatments\n\n")
  print(as.data.frame(x), row.names = FALSE)
  cat("\n")
  cat("Note: P-scores are frequentist analogues of SUCRA (Surface Under\n")
  cat("the Cumulative RAnking curve) from Bayesian network meta-analysis.\n")
  invisible(x)
}

#' Create League Table for Pairwise Comparisons
#'
#' Generates a league table showing all pairwise treatment comparisons with
#' effect estimates and confidence intervals. Standard format used in
#' BMJ, Lancet, and other top-tier journals.
#'
#' @param netmeta_obj A netmeta object from netmeta::netmeta()
#' @param digits Number of decimal places (default: 2)
#' @return Matrix with pairwise comparisons
#' @export
#' @references
#' BMJ Best Practice guidelines for network meta-analysis reporting.
#'
#' Salanti G (2012). Indirect and mixed-treatment comparison, network, or
#' multiple-treatments meta-analysis: many names, many benefits, many concerns
#' for the next generation evidence synthesis tool. Research Synthesis Methods,
#' 3(2):80-97.
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(30)
#' config <- setup_cnma(use_bayesian = FALSE)
#' results <- run_cnma_analysis(data, config = config)
#' league <- create_league_table(results$results$main_nma)
#' print(league)
#' }
create_league_table <- function(netmeta_obj, digits = 2) {
  if (!requireNamespace("netmeta", quietly = TRUE)) {
    .stop_hint("Package 'netmeta' is required.",
               "Install it with: install.packages('netmeta')")
  }

  if (inherits(netmeta_obj, "try-error") || is.null(netmeta_obj)) {
    .stop_hint("Invalid netmeta object provided.")
  }

  # Extract treatment effects and confidence intervals
  te <- netmeta_obj$TE.random
  lower <- netmeta_obj$lower.random
  upper <- netmeta_obj$upper.random
  treatments <- rownames(te)
  n_treat <- length(treatments)

  # Create league table matrix
  league <- matrix("", nrow = n_treat, ncol = n_treat)
  rownames(league) <- treatments
  colnames(league) <- treatments

  # Fill upper triangle with effect estimates (row vs column)
  for (i in 1:(n_treat - 1)) {
    for (j in (i + 1):n_treat) {
      # Effect of treatment i vs j
      effect <- te[i, j]
      ci_low <- lower[i, j]
      ci_high <- upper[i, j]

      # Format: "effect (CI_low, CI_high)"
      league[i, j] <- sprintf("%.*f (%.*f, %.*f)",
                             digits, effect,
                             digits, ci_low,
                             digits, ci_high)
    }
  }

  # Fill lower triangle with inverse comparisons
  for (i in 2:n_treat) {
    for (j in 1:(i - 1)) {
      # Effect of treatment i vs j (reversed)
      effect <- -te[j, i]
      ci_low <- -upper[j, i]
      ci_high <- -lower[j, i]

      league[i, j] <- sprintf("%.*f (%.*f, %.*f)",
                             digits, effect,
                             digits, ci_low,
                             digits, ci_high)
    }
  }

  # Diagonal: reference
  diag(league) <- treatments

  # Add class
  class(league) <- c("cnma_league_table", "matrix")

  return(league)
}

#' Print League Table
#' @param x cnma_league_table object
#' @param ... Additional arguments
#' @export
print.cnma_league_table <- function(x, ...) {
  cat("League Table: Pairwise Treatment Comparisons\n")
  cat("=============================================\n\n")
  cat("Upper triangle: Row treatment vs Column treatment\n")
  cat("Lower triangle: Column treatment vs Row treatment\n")
  cat("Format: Effect estimate (95% CI lower, upper)\n\n")

  # Print as data frame for better formatting
  print(noquote(x))

  cat("\n")
  cat("Note: Positive values favor row treatment (upper triangle)\n")
  cat("      or column treatment (lower triangle).\n")

  invisible(x)
}

#' Calculate Prediction Intervals
#'
#' Computes prediction intervals for network meta-analysis, accounting for
#' between-study heterogeneity. Recommended by Cochrane and PRISMA-NMA.
#'
#' @param netmeta_obj A netmeta object
#' @param level Confidence level (default: 0.95)
#' @return Data frame with prediction intervals for all comparisons
#' @export
#' @references
#' Riley RD, Higgins JP, Deeks JJ (2011). Interpretation of random effects
#' meta-analyses. BMJ, 342:d549.
#'
#' IntHout J, Ioannidis JP, Rovers MM, Goeman JJ (2016). Plea for routinely
#' presenting prediction intervals in meta-analysis. BMJ Open, 6(7):e010247.
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(30)
#' config <- setup_cnma(use_bayesian = FALSE)
#' results <- run_cnma_analysis(data, config = config)
#' pred_int <- calculate_prediction_intervals(results$results$main_nma)
#' print(pred_int)
#' }
calculate_prediction_intervals <- function(netmeta_obj, level = 0.95) {
  if (inherits(netmeta_obj, "try-error") || is.null(netmeta_obj)) {
    .stop_hint("Invalid netmeta object provided.")
  }

  # Extract components
  te <- netmeta_obj$TE.random
  tau <- netmeta_obj$tau
  treatments <- rownames(te)
  n_treat <- length(treatments)

  # Standard errors for each comparison
  seTE_matrix <- netmeta_obj$seTE.random

  # Calculate prediction intervals
  # PI = TE ± t * sqrt(seTE^2 + tau^2)
  alpha <- 1 - level
  # Use t-distribution (conservative)
  # df approximated by number of studies minus number of treatments
  k <- length(unique(netmeta_obj$studlab))
  df <- max(k - n_treat, 3)  # Minimum df of 3
  t_crit <- qt(1 - alpha / 2, df)

  pred_results <- list()

  for (i in 1:(n_treat - 1)) {
    for (j in (i + 1):n_treat) {
      comparison <- paste(treatments[i], "vs", treatments[j])
      effect <- te[i, j]
      se <- seTE_matrix[i, j]

      # Prediction interval
      pred_se <- sqrt(se^2 + tau^2)
      pi_lower <- effect - t_crit * pred_se
      pi_upper <- effect + t_crit * pred_se

      # Confidence interval (for comparison)
      ci_lower <- netmeta_obj$lower.random[i, j]
      ci_upper <- netmeta_obj$upper.random[i, j]

      pred_results[[comparison]] <- data.frame(
        comparison = comparison,
        treatment_1 = treatments[i],
        treatment_2 = treatments[j],
        effect = effect,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        pi_lower = pi_lower,
        pi_upper = pi_upper,
        stringsAsFactors = FALSE
      )
    }
  }

  pred_table <- do.call(rbind, pred_results)
  rownames(pred_table) <- NULL

  class(pred_table) <- c("cnma_prediction_intervals", "data.frame")

  return(pred_table)
}

#' Print Prediction Intervals
#' @param x cnma_prediction_intervals object
#' @param digits Number of decimal places (default: 3)
#' @param ... Additional arguments
#' @export
print.cnma_prediction_intervals <- function(x, digits = 3, ...) {
  cat("Prediction Intervals for Treatment Comparisons\n")
  cat("===============================================\n\n")
  cat("CI: Confidence Interval (uncertainty in mean effect)\n")
  cat("PI: Prediction Interval (range for future study)\n\n")

  # Format output
  output <- x
  output$effect <- round(output$effect, digits)
  output$ci_lower <- round(output$ci_lower, digits)
  output$ci_upper <- round(output$ci_upper, digits)
  output$pi_lower <- round(output$pi_lower, digits)
  output$pi_upper <- round(output$pi_upper, digits)

  print(output, row.names = FALSE)

  cat("\n")
  cat("Note: Prediction intervals are wider than confidence intervals\n")
  cat("      because they account for between-study heterogeneity.\n")
  cat("      PI represents the expected range of true effects in future studies.\n")

  invisible(x)
}
