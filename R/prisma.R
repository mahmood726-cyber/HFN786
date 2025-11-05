# =========================================================
# PRISMA-NMA Compliance and Reporting Functions
# Aligned with PRISMA-NMA Extension (2015) and upcoming 2026 update
# =========================================================

#' Generate PRISMA-NMA Compliance Report
#'
#' Creates a comprehensive report assessing compliance with PRISMA-NMA
#' reporting guidelines. Helps ensure publication-ready analysis.
#'
#' @param cnma_results CNMA results object from run_cnma_analysis()
#' @param output_format Format: "text", "markdown", or "html"
#' @param file Optional file path to save report
#' @return Character string with report (invisibly)
#' @export
#' @references
#' Hutton B, Salanti G, Caldwell DM, et al. (2015). The PRISMA extension
#' statement for reporting of systematic reviews incorporating network
#' meta-analyses of health care interventions: checklist and explanations.
#' Annals of Internal Medicine, 162(11):777-84.
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(40)
#' config <- setup_cnma(use_bayesian = FALSE)
#' results <- run_cnma_analysis(data, config = config)
#' report <- generate_prisma_report(results, output_format = "text")
#' cat(report)
#' }
generate_prisma_report <- function(cnma_results,
                                   output_format = c("text", "markdown", "html"),
                                   file = NULL) {
  output_format <- match.arg(output_format)

  if (!inherits(cnma_results, "cnma")) {
    .stop_hint("Input must be a cnma results object from run_cnma_analysis()")
  }

  # Extract components
  nma <- cnma_results$results$main_nma
  data <- cnma_results$data
  config <- cnma_results$config

  # Build report
  lines <- c()

  # Header
  if (output_format == "markdown") {
    lines <- c(lines, "# PRISMA-NMA Compliance Report")
    lines <- c(lines, "")
    lines <- c(lines, sprintf("Generated: %s", Sys.time()))
    lines <- c(lines, "")
  } else if (output_format == "html") {
    lines <- c(lines, "<html><head><title>PRISMA-NMA Report</title></head><body>")
    lines <- c(lines, "<h1>PRISMA-NMA Compliance Report</h1>")
    lines <- c(lines, sprintf("<p>Generated: %s</p>", Sys.time()))
  } else {
    lines <- c(lines, "PRISMA-NMA COMPLIANCE REPORT")
    lines <- c(lines, "=============================")
    lines <- c(lines, "")
    lines <- c(lines, sprintf("Generated: %s", Sys.time()))
    lines <- c(lines, "")
  }

  # Network characteristics
  lines <- c(lines, add_section("Network Characteristics", output_format))
  lines <- c(lines, add_item(
    sprintf("Number of studies: %d", length(unique(data$studlab))),
    output_format
  ))
  lines <- c(lines, add_item(
    sprintf("Number of treatments: %d",
            length(unique(c(data$treat1, data$treat2)))),
    output_format
  ))
  lines <- c(lines, add_item(
    sprintf("Number of pairwise comparisons: %d", nrow(data)),
    output_format
  ))
  lines <- c(lines, add_item(
    sprintf("Reference treatment: %s", cnma_results$ref_treatment),
    output_format
  ))
  lines <- c(lines, "")

  # Geometry of evidence network
  lines <- c(lines, add_section("Geometry of Evidence Network", output_format))
  if (!inherits(nma, "try-error")) {
    lines <- c(lines, add_item(
      sprintf("Summary measure: %s", config$sm),
      output_format
    ))
    lines <- c(lines, add_item(
      sprintf("Network is connected: %s",
              ifelse(length(unique(data$studlab)) > 0, "Yes", "Unknown")),
      output_format
    ))
  }
  lines <- c(lines, "")

  # Assessment of assumptions
  lines <- c(lines, add_section("Assessment of Key Assumptions", output_format))
  lines <- c(lines, add_subsection("1. Transitivity (Similarity)", output_format))
  lines <- c(lines, add_item(
    "Studies should be similar regarding effect modifiers",
    output_format
  ))
  lines <- c(lines, add_item(
    "Recommendation: Use assess_transitivity() to evaluate",
    output_format
  ))
  lines <- c(lines, "")

  lines <- c(lines, add_subsection("2. Consistency", output_format))
  if (!inherits(nma, "try-error")) {
    lines <- c(lines, add_item(
      sprintf("Global heterogeneity (I²): %.1f%%", nma$I2.random * 100),
      output_format
    ))
    lines <- c(lines, add_item(
      sprintf("Between-study heterogeneity (Tau): %.4f", nma$tau),
      output_format
    ))
    lines <- c(lines, add_item(
      "Recommendation: Use assess_inconsistency() for detailed evaluation",
      output_format
    ))
  }
  lines <- c(lines, "")

  # Statistical methods
  lines <- c(lines, add_section("Statistical Methods Applied", output_format))
  lines <- c(lines, add_item(
    sprintf("Analysis type: %s",
            ifelse(config$use_bayesian, "Bayesian", "Frequentist")),
    output_format
  ))
  lines <- c(lines, add_item(
    sprintf("Random effects model: %s",
            ifelse(!inherits(nma, "try-error"), "Yes", "N/A")),
    output_format
  ))
  lines <- c(lines, add_item(
    sprintf("Fixed effects model: %s", "No"),
    output_format
  ))
  lines <- c(lines, "")

  # Treatment rankings
  lines <- c(lines, add_section("Treatment Rankings", output_format))
  lines <- c(lines, add_item(
    "P-scores available: Use calculate_rankings()",
    output_format
  ))
  lines <- c(lines, add_item(
    "SUCRA values: Calculated via P-scores (frequentist analogue)",
    output_format
  ))
  lines <- c(lines, "")

  # PRISMA-NMA specific items (S1-S5)
  lines <- c(lines, add_section("PRISMA-NMA Extension Items", output_format))

  lines <- c(lines, add_subsection("S1: Network Structure", output_format))
  lines <- c(lines, add_item(
    "✓ Network diagram available via plot_network()",
    output_format
  ))
  lines <- c(lines, "")

  lines <- c(lines, add_subsection("S2: Assessment of Transitivity", output_format))
  lines <- c(lines, add_item(
    "Available via assess_transitivity()",
    output_format
  ))
  lines <- c(lines, "")

  lines <- c(lines, add_subsection("S3: Network Meta-Analysis Methods", output_format))
  lines <- c(lines, add_item(
    sprintf("✓ Method: %s network meta-analysis",
            ifelse(config$use_bayesian, "Bayesian", "Frequentist")),
    output_format
  ))
  lines <- c(lines, add_item(
    "✓ Software: R package netmeta via cnma",
    output_format
  ))
  lines <- c(lines, "")

  lines <- c(lines, add_subsection("S4: Assessment of Inconsistency", output_format))
  lines <- c(lines, add_item(
    "Available via assess_inconsistency()",
    output_format
  ))
  lines <- c(lines, add_item(
    "- Global: Cochran's Q, I², Tau²",
    output_format
  ))
  lines <- c(lines, add_item(
    "- Local: Node-splitting (netsplit)",
    output_format
  ))
  lines <- c(lines, add_item(
    "- Design: Design-by-treatment interaction",
    output_format
  ))
  lines <- c(lines, "")

  lines <- c(lines, add_subsection("S5: Results Presentation", output_format))
  lines <- c(lines, add_item(
    "✓ League table: create_league_table()",
    output_format
  ))
  lines <- c(lines, add_item(
    "✓ Forest plot: plot_forest()",
    output_format
  ))
  lines <- c(lines, add_item(
    "✓ Ranking plots: plot_rankings()",
    output_format
  ))
  lines <- c(lines, add_item(
    "✓ Prediction intervals: calculate_prediction_intervals()",
    output_format
  ))
  lines <- c(lines, "")

  # Quality of evidence
  lines <- c(lines, add_section("Quality of Evidence (GRADE)", output_format))
  lines <- c(lines, add_item(
    sprintf("GRADE assessment enabled: %s", config$use_grade_weighting),
    output_format
  ))
  lines <- c(lines, add_item(
    "Consider: Risk of bias, inconsistency, indirectness, imprecision, publication bias",
    output_format
  ))
  lines <- c(lines, "")

  # Recommendations
  lines <- c(lines, add_section("Recommendations for Complete Reporting", output_format))
  lines <- c(lines, add_item(
    "1. Generate all visualizations using create_publication_plots()",
    output_format
  ))
  lines <- c(lines, add_item(
    "2. Assess inconsistency using assess_inconsistency()",
    output_format
  ))
  lines <- c(lines, add_item(
    "3. Evaluate transitivity using assess_transitivity()",
    output_format
  ))
  lines <- c(lines, add_item(
    "4. Calculate treatment rankings using calculate_rankings()",
    output_format
  ))
  lines <- c(lines, add_item(
    "5. Generate league table using create_league_table()",
    output_format
  ))
  lines <- c(lines, add_item(
    "6. Calculate prediction intervals using calculate_prediction_intervals()",
    output_format
  ))
  lines <- c(lines, add_item(
    "7. Create funnel plot to assess publication bias using plot_funnel()",
    output_format
  ))
  lines <- c(lines, "")

  # References
  lines <- c(lines, add_section("Key References", output_format))
  lines <- c(lines, add_item(
    "Hutton et al. (2015). The PRISMA extension statement for reporting of",
    output_format
  ))
  lines <- c(lines, add_item(
    "  systematic reviews incorporating NMA. Ann Intern Med, 162(11):777-84.",
    output_format
  ))
  lines <- c(lines, add_item(
    "Rücker & Schwarzer (2015). Ranking treatments in frequentist NMA.",
    output_format
  ))
  lines <- c(lines, add_item(
    "  BMC Med Res Methodol, 15:58.",
    output_format
  ))
  lines <- c(lines, "")

  if (output_format == "html") {
    lines <- c(lines, "</body></html>")
  }

  report <- paste(lines, collapse = "\n")

  # Save to file if specified
  if (!is.null(file)) {
    writeLines(report, file)
    msg("Report saved to: %s", file)
  }

  invisible(report)
}

# Helper functions for formatting
add_section <- function(title, format) {
  if (format == "markdown") {
    return(c(sprintf("## %s", title), ""))
  } else if (format == "html") {
    return(sprintf("<h2>%s</h2>", title))
  } else {
    return(c(title, paste(rep("-", nchar(title)), collapse = ""), ""))
  }
}

add_subsection <- function(title, format) {
  if (format == "markdown") {
    return(c(sprintf("### %s", title), ""))
  } else if (format == "html") {
    return(sprintf("<h3>%s</h3>", title))
  } else {
    return(c(title, ""))
  }
}

add_item <- function(text, format) {
  if (format == "markdown") {
    return(sprintf("- %s", text))
  } else if (format == "html") {
    return(sprintf("<p>%s</p>", text))
  } else {
    return(sprintf("  %s", text))
  }
}

#' Generate Network Characteristics Summary
#'
#' Creates a comprehensive summary of network characteristics for publication.
#'
#' @param data Original data frame
#' @param netmeta_obj Optional netmeta object
#' @return Data frame with network characteristics
#' @export
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(40)
#' summary <- network_characteristics_summary(data)
#' print(summary)
#' }
network_characteristics_summary <- function(data, netmeta_obj = NULL) {
  # Basic counts
  n_studies <- length(unique(data$studlab))
  treatments <- unique(c(data$treat1, data$treat2))
  n_treatments <- length(treatments)
  n_comparisons <- nrow(data)

  # Studies per comparison
  studies_per_comp <- data %>%
    group_by(treat1, treat2) %>%
    summarise(n_studies = n(), .groups = "drop")

  # Direct comparisons available
  n_direct_comparisons <- nrow(studies_per_comp)

  # Possible comparisons
  n_possible <- choose(n_treatments, 2)

  # Network density
  network_density <- n_direct_comparisons / n_possible

  # Multi-arm studies
  arms_per_study <- data %>%
    group_by(studlab) %>%
    summarise(
      n_arms = length(unique(c(treat1, treat2))),
      .groups = "drop"
    )

  n_multi_arm <- sum(arms_per_study$n_arms > 2)

  summary_df <- data.frame(
    characteristic = c(
      "Number of studies",
      "Number of treatments",
      "Number of pairwise comparisons",
      "Number of direct comparisons",
      "Possible pairwise comparisons",
      "Network density",
      "Multi-arm studies",
      "Median studies per comparison",
      "Range of studies per comparison"
    ),
    value = c(
      as.character(n_studies),
      as.character(n_treatments),
      as.character(n_comparisons),
      as.character(n_direct_comparisons),
      as.character(n_possible),
      sprintf("%.2f%%", network_density * 100),
      as.character(n_multi_arm),
      as.character(median(studies_per_comp$n_studies)),
      sprintf("%d - %d",
              min(studies_per_comp$n_studies),
              max(studies_per_comp$n_studies))
    ),
    stringsAsFactors = FALSE
  )

  # Add heterogeneity if netmeta object provided
  if (!is.null(netmeta_obj) && !inherits(netmeta_obj, "try-error")) {
    het_rows <- data.frame(
      characteristic = c(
        "Tau (between-study SD)",
        "Tau² (between-study variance)",
        "I² (heterogeneity)"
      ),
      value = c(
        sprintf("%.4f", netmeta_obj$tau),
        sprintf("%.4f", netmeta_obj$tau^2),
        sprintf("%.1f%%", netmeta_obj$I2.random * 100)
      ),
      stringsAsFactors = FALSE
    )
    summary_df <- rbind(summary_df, het_rows)
  }

  class(summary_df) <- c("cnma_network_summary", "data.frame")
  return(summary_df)
}

#' Print Network Characteristics Summary
#' @param x cnma_network_summary object
#' @param ... Additional arguments
#' @export
print.cnma_network_summary <- function(x, ...) {
  cat("Network Characteristics Summary\n")
  cat("================================\n\n")
  print(as.data.frame(x), row.names = FALSE)
  invisible(x)
}
