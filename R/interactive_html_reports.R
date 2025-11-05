# Interactive HTML Report Generation
# Complete interactive reports with all visualizations and results
# Version 1.7.0

# =============================================================================
# Comprehensive Interactive HTML Report
# =============================================================================

#' Generate Comprehensive Interactive HTML Report
#'
#' Creates a complete, interactive HTML report with all network meta-analysis
#' results, visualizations, tables, and statistical summaries. Includes
#' interactive plots, collapsible sections, and professional styling.
#'
#' @param nma_results Network meta-analysis results
#' @param data Study-level data
#' @param title Report title
#' @param author Report author
#' @param output_file Output HTML file path
#' @param include_ml Include machine learning predictions
#' @param include_grade Include GRADE tables
#' @param theme Report theme: "light", "dark", "blue", "green"
#' @param toc_float Floating table of contents
#'
#' @return Path to generated HTML report
#' @export
#'
#' @examples
#' \dontrun{
#' data <- simulate_cnma_data(50)
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#'
#' report <- generate_interactive_html_report(
#'   nma, data,
#'   title = "My Network Meta-Analysis",
#'   author = "Research Team",
#'   output_file = "nma_report.html",
#'   include_ml = TRUE,
#'   include_grade = TRUE
#' )
#' }
generate_interactive_html_report <- function(nma_results,
                                            data,
                                            title = "Network Meta-Analysis Report",
                                            author = "CNMA Analysis",
                                            output_file = "nma_report.html",
                                            include_ml = TRUE,
                                            include_grade = TRUE,
                                            theme = c("light", "dark", "blue", "green"),
                                            toc_float = TRUE) {

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' required. Install with: install.packages('rmarkdown')")
  }

  theme <- match.arg(theme)

  # Create temporary directory for assets
  temp_dir <- tempdir()
  assets_dir <- file.path(temp_dir, "report_assets")
  if (!dir.exists(assets_dir)) dir.create(assets_dir, recursive = TRUE)

  # Generate all visualizations
  message("Generating visualizations...")
  viz_files <- .generate_all_visualizations(nma_results, data, assets_dir)

  # Generate machine learning results if requested
  ml_results <- NULL
  if (include_ml && "year" %in% names(data)) {
    message("Running machine learning predictions...")
    ml_results <- tryCatch({
      ml_predict_rankings(data, nma_results, method = "rf")
    }, error = function(e) NULL)
  }

  # Generate GRADE table if requested
  grade_table <- NULL
  if (include_grade) {
    message("Generating GRADE evidence profile...")
    grade_table <- tryCatch({
      generate_grade_table(nma_results, data, output_file = NULL)
    }, error = function(e) NULL)
  }

  # Create R Markdown content
  message("Creating HTML report...")
  rmd_content <- .create_interactive_report_rmd(
    nma_results, data, title, author, theme, toc_float,
    viz_files, ml_results, grade_table
  )

  # Write R Markdown file
  rmd_file <- file.path(temp_dir, "report.Rmd")
  writeLines(rmd_content, rmd_file)

  # Render to HTML
  rmarkdown::render(
    rmd_file,
    output_format = rmarkdown::html_document(
      toc = TRUE,
      toc_float = toc_float,
      toc_depth = 3,
      theme = switch(theme,
        "light" = "flatly",
        "dark" = "darkly",
        "blue" = "cerulean",
        "green" = "journal"
      ),
      highlight = "tango",
      code_folding = "hide",
      df_print = "paged"
    ),
    output_file = output_file,
    quiet = TRUE
  )

  message("Interactive HTML report generated: ", output_file)

  return(output_file)
}

# =============================================================================
# Helper Functions
# =============================================================================

.generate_all_visualizations <- function(nma_results, data, output_dir) {

  viz_files <- list()

  # Network plot
  viz_files$network <- file.path(output_dir, "network.png")
  png(viz_files$network, width = 1200, height = 900, res = 150)
  netmeta::netgraph(nma_results, plastic = FALSE, thickness = "number.of.studies")
  dev.off()

  # Forest plot
  viz_files$forest <- file.path(output_dir, "forest.png")
  png(viz_files$forest, width = 1400, height = 1000, res = 150)
  netmeta::forest(nma_results, xlim = c(0.5, 2))
  dev.off()

  # Ranking plot
  viz_files$ranking <- file.path(output_dir, "ranking.png")
  rankings <- netmeta::netrank(nma_results)
  png(viz_files$ranking, width = 1200, height = 800, res = 150)
  plot(rankings)
  dev.off()

  # Net heat plot
  viz_files$netheat <- file.path(output_dir, "netheat.png")
  png(viz_files$netheat, width = 1400, height = 1000, res = 150)
  netmeta::netheat(nma_results)
  dev.off()

  # Funnel plot
  viz_files$funnel <- file.path(output_dir, "funnel.png")
  png(viz_files$funnel, width = 1200, height = 900, res = 150)
  netmeta::funnel(nma_results)
  dev.off()

  # Contribution matrix
  viz_files$contribution <- file.path(output_dir, "contribution.png")
  png(viz_files$contribution, width = 1400, height = 1000, res = 150)
  contrib <- netmeta::netcontrib(nma_results)
  plot(contrib)
  dev.off()

  return(viz_files)
}

.create_interactive_report_rmd <- function(nma_results, data, title, author,
                                          theme, toc_float, viz_files,
                                          ml_results, grade_table) {

  # Extract key statistics
  n_studies <- length(unique(data$studlab))
  n_treatments <- length(nma_results$trts)
  n_comparisons <- nrow(data)

  rankings <- netmeta::netrank(nma_results)
  top_treatment <- nma_results$trts[which.max(rankings$Pscore.random)]

  # Create R Markdown content
  rmd <- paste0(
    "---\n",
    "title: \"", title, "\"\n",
    "author: \"", author, "\"\n",
    "date: \"`r Sys.Date()`\"\n",
    "output:\n",
    "  html_document:\n",
    "    toc: true\n",
    "    toc_float: ", tolower(as.character(toc_float)), "\n",
    "    toc_depth: 3\n",
    "    theme: ", switch(theme, "light" = "flatly", "dark" = "darkly",
                          "blue" = "cerulean", "green" = "journal"), "\n",
    "    code_folding: hide\n",
    "---\n\n",

    "```{r setup, include=FALSE}\n",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)\n",
    "```\n\n",

    "# Executive Summary {.tabset}\n\n",
    "## Overview\n\n",
    "This network meta-analysis synthesizes evidence from **", n_studies, " studies** comparing **",
    n_treatments, " treatments** across **", n_comparisons, " comparisons**.\n\n",

    "### Key Findings\n\n",
    "- **Top-ranked treatment**: ", top_treatment, " (P-score: ",
    round(max(rankings$Pscore.random), 3), ")\n",
    "- **Between-study heterogeneity (I²)**: ", round(nma_results$I2 * 100, 1), "%\n",
    "- **Between-study variance (Tau²)**: ", round(nma_results$tau^2, 4), "\n",
    "- **Global inconsistency (Q)**: ", round(nma_results$Q, 2),
    " (p = ", format.pval(nma_results$pval.Q, digits = 3), ")\n\n",

    "## Network Characteristics\n\n",
    "```{r network-table}\n",
    "network_char <- data.frame(\n",
    "  Characteristic = c('Number of Studies', 'Number of Treatments', 'Number of Comparisons',\n",
    "                     'Effect Measure', 'Model', 'Heterogeneity (I²)', 'Tau²'),\n",
    "  Value = c('", n_studies, "', '", n_treatments, "', '", n_comparisons, "',\n",
    "            '", nma_results$sm, "', 'Random Effects', '",
    round(nma_results$I2 * 100, 1), "%', '", round(nma_results$tau^2, 4), "')\n",
    ")\n",
    "knitr::kable(network_char, caption = 'Network Meta-Analysis Characteristics')\n",
    "```\n\n",

    "# Visualizations {.tabset}\n\n",
    "## Network Plot\n\n",
    "![Treatment Network](", viz_files$network, "){width=100%}\n\n",
    "The network plot shows all treatment comparisons with edge thickness representing the number of studies.\n\n",

    "## Forest Plot\n\n",
    "![Forest Plot of Treatment Effects](", viz_files$forest, "){width=100%}\n\n",
    "Forest plot showing treatment effects compared to reference with 95% confidence intervals and prediction intervals.\n\n",

    "## Treatment Rankings\n\n",
    "![Treatment Ranking Plot](", viz_files$ranking, "){width=100%}\n\n",
    "Treatment ranking based on P-scores (frequentist analog of SUCRA).\n\n",

    "## Net Heat Plot\n\n",
    "![Net Heat Plot for Inconsistency](", viz_files$netheat, "){width=100%}\n\n",
    "Net heat plot for detecting inconsistency in the network. Hot colors indicate potential inconsistency.\n\n",

    "## Funnel Plot\n\n",
    "![Comparison-Adjusted Funnel Plot](", viz_files$funnel, "){width=100%}\n\n",
    "Comparison-adjusted funnel plot for assessing publication bias.\n\n",

    "## Contribution Matrix\n\n",
    "![Contribution of Direct Comparisons](", viz_files$contribution, "){width=100%}\n\n",
    "Contribution matrix showing how direct comparisons contribute to network estimates.\n\n",

    "# Results Tables {.tabset}\n\n",
    "## Treatment Rankings\n\n",
    "```{r rankings-table}\n",
    "rankings_df <- data.frame(\n",
    "  Treatment = c('", paste(nma_results$trts, collapse = "', '"), "'),\n",
    "  P_score = c(", paste(round(rankings$Pscore.random, 4), collapse = ", "), "),\n",
    "  Rank = c(", paste(rank(-rankings$Pscore.random), collapse = ", "), ")\n",
    ")\n",
    "rankings_df <- rankings_df[order(-rankings_df$P_score), ]\n",
    "knitr::kable(rankings_df, caption = 'Treatment Rankings', row.names = FALSE)\n",
    "```\n\n",

    "## League Table\n\n",
    "```{r league-table}\n",
    "league <- netmeta::netleague(nma_results, digits = 2)\n",
    "knitr::kable(league$random, caption = 'League Table of Pairwise Comparisons')\n",
    "```\n\n"
  )

  # Add machine learning section if available
  if (!is.null(ml_results)) {
    rmd <- paste0(rmd,
      "# Machine Learning Predictions {.tabset}\n\n",
      "## Model Performance\n\n",
      "```{r ml-performance}\n",
      "ml_perf <- data.frame(\n",
      "  Metric = c('RMSE', 'MAE', 'R²'),\n",
      "  Value = c(", ml_results$performance$RMSE, ", ",
                   ml_results$performance$MAE, ", ",
                   ml_results$performance$R_squared, ")\n",
      ")\n",
      "knitr::kable(ml_perf, caption = 'Machine Learning Model Performance', digits = 4)\n",
      "```\n\n",

      "## Variable Importance\n\n",
      "```{r ml-importance}\n",
      "importance_df <- head(ml_results$importance, 10)\n",
      "knitr::kable(importance_df, caption = 'Top 10 Important Variables', row.names = FALSE)\n",
      "```\n\n"
    )
  }

  # Add GRADE section if available
  if (!is.null(grade_table)) {
    rmd <- paste0(rmd,
      "# GRADE Evidence Profile\n\n",
      "```{r grade-table}\n",
      "knitr::kable(grade_table, caption = 'GRADE Evidence Profile Table')\n",
      "```\n\n"
    )
  }

  # Add heterogeneity and inconsistency section
  rmd <- paste0(rmd,
    "# Statistical Assessment {.tabset}\n\n",
    "## Heterogeneity\n\n",
    "### Summary Statistics\n\n",
    "```{r heterogeneity}\n",
    "het_stats <- data.frame(\n",
    "  Statistic = c('Tau² (between-study variance)', 'Tau (between-study SD)',\n",
    "                'I² (heterogeneity)', 'Q statistic', 'Q p-value'),\n",
    "  Value = c('", round(nma_results$tau^2, 4), "', '", round(nma_results$tau, 4), "',\n",
    "            '", round(nma_results$I2 * 100, 2), "%', '", round(nma_results$Q, 2), "',\n",
    "            '", format.pval(nma_results$pval.Q, digits = 4), "')\n",
    ")\n",
    "knitr::kable(het_stats, caption = 'Heterogeneity Statistics')\n",
    "```\n\n",

    "### Interpretation\n\n",
    ifelse(nma_results$I2 < 0.25,
      "The I² statistic suggests **low heterogeneity** (< 25%).",
    ifelse(nma_results$I2 < 0.50,
      "The I² statistic suggests **moderate heterogeneity** (25-50%).",
    ifelse(nma_results$I2 < 0.75,
      "The I² statistic suggests **substantial heterogeneity** (50-75%).",
      "The I² statistic suggests **considerable heterogeneity** (> 75%)."
    ))), "\n\n",

    "## Inconsistency\n\n",
    "```{r inconsistency}\n",
    "decomp <- netmeta::decomp.design(nma_results)\n",
    "inconsis_stats <- data.frame(\n",
    "  Statistic = c('Q total', 'Q within designs', 'Q between designs', 'p-value (within)', 'p-value (between)'),\n",
    "  Value = c(round(decomp$Q.decomp['Total'], 2),\n",
    "            round(decomp$Q.decomp['Within designs'], 2),\n",
    "            round(decomp$Q.decomp['Between designs'], 2),\n",
    "            format.pval(decomp$pval.Q.decomp['Within designs'], digits = 4),\n",
    "            format.pval(decomp$pval.Q.decomp['Between designs'], digits = 4))\n",
    ")\n",
    "knitr::kable(inconsis_stats, caption = 'Design-Based Inconsistency Decomposition')\n",
    "```\n\n",

    "# Study Characteristics\n\n",
    "```{r study-char}\n",
    "study_summary <- data %>%\n",
    "  dplyr::group_by(studlab) %>%\n",
    "  dplyr::summarise(\n",
    "    Comparisons = dplyr::n(),\n",
    "    Treatments = paste(unique(c(treat1, treat2)), collapse = ', '),\n",
    "    .groups = 'drop'\n",
    "  )\n",
    "knitr::kable(study_summary, caption = 'Study Characteristics')\n",
    "```\n\n",

    "# Methods\n\n",
    "## Statistical Methods\n\n",
    "This network meta-analysis was conducted using frequentist methods with the `netmeta` package in R.\n\n",
    "- **Effect measure**: ", nma_results$sm, "\n",
    "- **Pooling model**: Random effects model\n",
    "- **Heterogeneity estimation**: Restricted maximum likelihood (REML)\n",
    "- **Inconsistency assessment**: Design-by-treatment interaction model\n",
    "- **Treatment ranking**: P-scores (frequentist analog of SUCRA)\n\n",
    "## Software\n\n",
    "- R version: `r R.version.string`\n",
    "- netmeta package version: `r packageVersion('netmeta')`\n",
    "- CNMA package version: 1.7.0\n\n",

    "# References\n\n",
    "1. Rücker G, Schwarzer G. Ranking treatments in frequentist network meta-analysis works without resampling methods. *BMC Med Res Methodol* 2015;15:58.\n\n",
    "2. Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed treatment comparison meta-analysis. *Stat Med* 2010;29(7-8):932-944.\n\n",
    "3. Hutton B, Salanti G, Caldwell DM, et al. The PRISMA extension statement for reporting of systematic reviews incorporating network meta-analyses. *Ann Intern Med* 2015;162(11):777-784.\n\n",

    "---\n\n",
    "**Report generated on**: `r Sys.time()`\n\n",
    "**Generated by**: CNMA Package v1.7.0\n"
  )

  return(rmd)
}

# =============================================================================
# Interactive Dashboard HTML Report
# =============================================================================

#' Generate Interactive Dashboard HTML
#'
#' Creates an interactive dashboard-style HTML report with tabs, widgets,
#' and real-time interactivity using htmlwidgets.
#'
#' @param nma_results Network meta-analysis results
#' @param data Study-level data
#' @param output_file Output HTML file
#'
#' @return Path to generated HTML file
#' @export
#'
#' @examples
#' \dontrun{
#' generate_interactive_dashboard_html(nma, data, "dashboard.html")
#' }
generate_interactive_dashboard_html <- function(nma_results,
                                               data,
                                               output_file = "nma_dashboard.html") {

  if (!requireNamespace("flexdashboard", quietly = TRUE)) {
    message("Package 'flexdashboard' recommended for best experience.")
  }

  # Create dashboard with plotly interactive plots
  temp_dir <- tempdir()
  rmd_file <- file.path(temp_dir, "dashboard.Rmd")

  # Create dashboard content
  dashboard_content <- .create_dashboard_content(nma_results, data)

  writeLines(dashboard_content, rmd_file)

  # Render
  rmarkdown::render(
    rmd_file,
    output_file = output_file,
    quiet = TRUE
  )

  message("Interactive dashboard generated: ", output_file)
  return(output_file)
}

.create_dashboard_content <- function(nma_results, data) {

  rankings <- netmeta::netrank(nma_results)

  content <- paste0(
    "---\n",
    "title: 'Network Meta-Analysis Dashboard'\n",
    "output:\n",
    "  html_document:\n",
    "    theme: cosmo\n",
    "---\n\n",

    "```{r setup, include=FALSE}\n",
    "library(plotly)\n",
    "library(DT)\n",
    "```\n\n",

    "# {.tabset}\n\n",

    "## Summary\n\n",
    "### Key Metrics\n\n",
    "```{r summary-metrics, echo=FALSE}\n",
    "metrics <- data.frame(\n",
    "  Metric = c('Number of Studies', 'Number of Treatments', 'Heterogeneity (I²)', 'Top Treatment'),\n",
    "  Value = c('", length(unique(data$studlab)), "', '",
    length(nma_results$trts), "', '",
    round(nma_results$I2 * 100, 1), "%', '",
    nma_results$trts[which.max(rankings$Pscore.random)], "')\n",
    ")\n",
    "knitr::kable(metrics, format = 'html')\n",
    "```\n\n",

    "## Interactive Rankings\n\n",
    "```{r rankings-plotly, echo=FALSE}\n",
    "rankings_df <- data.frame(\n",
    "  Treatment = c('", paste(nma_results$trts, collapse = "', '"), "'),\n",
    "  P_score = c(", paste(round(rankings$Pscore.random, 4), collapse = ", "), ")\n",
    ")\n",
    "rankings_df <- rankings_df[order(-rankings_df$P_score), ]\n",
    "p <- plot_ly(rankings_df, x = ~reorder(Treatment, P_score), y = ~P_score, type = 'bar',\n",
    "             marker = list(color = 'rgb(158,202,225)', line = list(color = 'rgb(8,48,107)', width = 1.5)))\n",
    "p <- p %>% layout(title = 'Treatment Rankings (P-scores)',\n",
    "                  xaxis = list(title = 'Treatment'),\n",
    "                  yaxis = list(title = 'P-score'))\n",
    "p\n",
    "```\n\n",

    "## Data Table\n\n",
    "```{r data-table, echo=FALSE}\n",
    "DT::datatable(data, filter = 'top', options = list(pageLength = 10, autoWidth = TRUE))\n",
    "```\n\n"
  )

  return(content)
}

# =============================================================================
# Export Complete Report Package
# =============================================================================

#' Export Complete Report Package
#'
#' Creates a complete report package with HTML report, all figures, tables,
#' and data files in a ZIP archive.
#'
#' @param nma_results Network meta-analysis results
#' @param data Study-level data
#' @param project_name Project name for the package
#' @param output_dir Output directory
#'
#' @return Path to ZIP file
#' @export
#'
#' @examples
#' \dontrun{
#' export_complete_report_package(nma, data, "MyNMA", "reports")
#' }
export_complete_report_package <- function(nma_results,
                                          data,
                                          project_name = "NMA_Report",
                                          output_dir = ".") {

  # Create temporary directory structure
  temp_base <- tempfile()
  dir.create(temp_base, recursive = TRUE)

  package_dir <- file.path(temp_base, project_name)
  dir.create(package_dir)

  # Create subdirectories
  figures_dir <- file.path(package_dir, "figures")
  tables_dir <- file.path(package_dir, "tables")
  data_dir <- file.path(package_dir, "data")
  reports_dir <- file.path(package_dir, "reports")

  dir.create(figures_dir)
  dir.create(tables_dir)
  dir.create(data_dir)
  dir.create(reports_dir)

  message("Generating complete report package...")

  # Generate HTML report
  html_report <- file.path(reports_dir, "interactive_report.html")
  generate_interactive_html_report(
    nma_results, data,
    title = paste(project_name, "Analysis"),
    output_file = html_report,
    include_ml = TRUE,
    include_grade = TRUE
  )

  # Generate all figures
  .export_all_figures(nma_results, data, figures_dir)

  # Generate all tables
  .export_all_tables(nma_results, data, tables_dir)

  # Export data
  write.csv(data, file.path(data_dir, "study_data.csv"), row.names = FALSE)
  saveRDS(nma_results, file.path(data_dir, "nma_results.rds"))

  # Create README
  readme_content <- paste0(
    "# ", project_name, " - Network Meta-Analysis Report Package\n\n",
    "Generated: ", Sys.time(), "\n\n",
    "## Contents\n\n",
    "- **reports/**: Interactive HTML report\n",
    "- **figures/**: All visualizations (PNG format)\n",
    "- **tables/**: All results tables (CSV format)\n",
    "- **data/**: Original data and NMA results object\n\n",
    "## How to Use\n\n",
    "1. Open `reports/interactive_report.html` in your web browser\n",
    "2. View figures in the `figures/` directory\n",
    "3. Import tables from `tables/` directory into your manuscript\n",
    "4. Load data and results in R using `data/nma_results.rds`\n\n",
    "Generated with CNMA Package v1.7.0\n"
  )

  writeLines(readme_content, file.path(package_dir, "README.txt"))

  # Create ZIP archive
  zip_file <- file.path(output_dir, paste0(project_name, "_", format(Sys.Date(), "%Y%m%d"), ".zip"))

  current_dir <- getwd()
  setwd(temp_base)
  utils::zip(zip_file, files = project_name)
  setwd(current_dir)

  message("Complete report package exported: ", zip_file)

  return(zip_file)
}

.export_all_figures <- function(nma_results, data, output_dir) {

  # Network plot
  png(file.path(output_dir, "01_network_plot.png"), width = 1200, height = 900, res = 150)
  netmeta::netgraph(nma_results, plastic = FALSE, thickness = "number.of.studies")
  dev.off()

  # Forest plot
  png(file.path(output_dir, "02_forest_plot.png"), width = 1400, height = 1000, res = 150)
  netmeta::forest(nma_results)
  dev.off()

  # Ranking plot
  png(file.path(output_dir, "03_ranking_plot.png"), width = 1200, height = 800, res = 150)
  rankings <- netmeta::netrank(nma_results)
  plot(rankings)
  dev.off()

  # Net heat plot
  png(file.path(output_dir, "04_netheat_plot.png"), width = 1400, height = 1000, res = 150)
  netmeta::netheat(nma_results)
  dev.off()

  # Funnel plot
  png(file.path(output_dir, "05_funnel_plot.png"), width = 1200, height = 900, res = 150)
  netmeta::funnel(nma_results)
  dev.off()

  # Contribution matrix
  png(file.path(output_dir, "06_contribution_matrix.png"), width = 1400, height = 1000, res = 150)
  contrib <- netmeta::netcontrib(nma_results)
  plot(contrib)
  dev.off()

  message("All figures exported to: ", output_dir)
}

.export_all_tables <- function(nma_results, data, output_dir) {

  # Treatment rankings
  rankings <- netmeta::netrank(nma_results)
  rankings_df <- data.frame(
    Treatment = nma_results$trts,
    P_score = rankings$Pscore.random,
    Rank = rank(-rankings$Pscore.random)
  )
  write.csv(rankings_df, file.path(output_dir, "01_treatment_rankings.csv"), row.names = FALSE)

  # League table
  league <- netmeta::netleague(nma_results, digits = 3)
  write.csv(league$random, file.path(output_dir, "02_league_table.csv"))

  # Heterogeneity statistics
  het_stats <- data.frame(
    Statistic = c("Tau2", "Tau", "I2", "Q", "Q_pvalue"),
    Value = c(nma_results$tau^2, nma_results$tau, nma_results$I2,
              nma_results$Q, nma_results$pval.Q)
  )
  write.csv(het_stats, file.path(output_dir, "03_heterogeneity_statistics.csv"), row.names = FALSE)

  # Study summary
  study_summary <- data %>%
    dplyr::group_by(studlab) %>%
    dplyr::summarise(
      n_comparisons = dplyr::n(),
      treatments = paste(unique(c(treat1, treat2)), collapse = ", "),
      .groups = "drop"
    )
  write.csv(study_summary, file.path(output_dir, "04_study_summary.csv"), row.names = FALSE)

  message("All tables exported to: ", output_dir)
}
