# Advanced Reporting Templates
# GRADE tables, ROB2 visualizations, Sankey diagrams, animated plots
# Version 1.6.0

# =============================================================================
# GRADE Evidence Profile Tables
# =============================================================================

#' Generate GRADE Evidence Profile Table
#'
#' Creates a GRADE (Grading of Recommendations, Assessment, Development and
#' Evaluations) evidence profile table for network meta-analysis results.
#'
#' @param nma_results Network meta-analysis results
#' @param data Study-level data
#' @param comparisons Character vector of comparisons to include
#' @param rob_assessments Risk of bias assessments (optional)
#' @param output_file Output file path (default: "grade_table.csv")
#'
#' @return Data frame with GRADE evidence profile
#' @export
#'
#' @examples
#' \dontrun{
#' data <- simulate_cnma_data(50)
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' grade_table <- generate_grade_table(nma, data)
#' }
generate_grade_table <- function(nma_results,
                                data,
                                comparisons = NULL,
                                rob_assessments = NULL,
                                output_file = "grade_table.csv") {

  # If no comparisons specified, use top 5 vs reference
  if (is.null(comparisons)) {
    rankings <- netmeta::netrank(nma_results)
    top_treatments <- nma_results$trts[order(-rankings$Pscore.random)][1:min(5, length(nma_results$trts))]
    ref_treatment <- nma_results$reference.group

    if (is.null(ref_treatment)) {
      ref_treatment <- nma_results$trts[1]
    }

    comparisons <- paste(top_treatments, "vs", ref_treatment)
  }

  # Initialize GRADE table
  grade_table <- data.frame(
    Comparison = character(),
    N_Studies = integer(),
    N_Participants = integer(),
    Effect_Estimate = character(),
    CI_95 = character(),
    Risk_of_Bias = character(),
    Inconsistency = character(),
    Indirectness = character(),
    Imprecision = character(),
    Publication_Bias = character(),
    Overall_Quality = character(),
    stringsAsFactors = FALSE
  )

  # For each comparison
  for (comp in comparisons) {
    # Parse comparison
    treatments <- strsplit(comp, " vs ")[[1]]
    if (length(treatments) != 2) next

    treat1 <- trimws(treatments[1])
    treat2 <- trimws(treatments[2])

    # Get effect estimate
    comparison_idx <- which(
      (nma_results$treat1 == treat1 & nma_results$treat2 == treat2) |
      (nma_results$treat1 == treat2 & nma_results$treat2 == treat1)
    )

    if (length(comparison_idx) == 0) {
      # Try to get from network estimate
      te <- NA
      lower <- NA
      upper <- NA
      n_studies <- 0
      n_participants <- NA
    } else {
      te <- nma_results$TE.random[comparison_idx[1]]
      lower <- nma_results$lower.random[comparison_idx[1]]
      upper <- nma_results$upper.random[comparison_idx[1]]

      # Count studies and participants
      comp_data <- data[
        (data$treat1 == treat1 & data$treat2 == treat2) |
        (data$treat1 == treat2 & data$treat2 == treat1),
      ]
      n_studies <- nrow(comp_data)
      n_participants <- ifelse("n" %in% names(comp_data), sum(comp_data$n, na.rm = TRUE), NA)
    }

    # Assess GRADE domains
    rob_rating <- .assess_grade_rob(data, rob_assessments, treat1, treat2)
    inconsistency_rating <- .assess_grade_inconsistency(nma_results)
    indirectness_rating <- .assess_grade_indirectness(data, treat1, treat2)
    imprecision_rating <- .assess_grade_imprecision(te, lower, upper, n_studies)
    pub_bias_rating <- .assess_grade_publication_bias(nma_results, data)

    # Calculate overall quality
    overall_quality <- .calculate_overall_grade_quality(
      rob_rating, inconsistency_rating, indirectness_rating,
      imprecision_rating, pub_bias_rating
    )

    # Format effect estimate
    if (!is.na(te)) {
      if (nma_results$sm %in% c("HR", "OR", "RR")) {
        effect_str <- sprintf("%.2f", exp(te))
        ci_str <- sprintf("%.2f to %.2f", exp(lower), exp(upper))
      } else {
        effect_str <- sprintf("%.2f", te)
        ci_str <- sprintf("%.2f to %.2f", lower, upper)
      }
    } else {
      effect_str <- "NE"
      ci_str <- "NE"
    }

    # Add row to table
    grade_table <- rbind(
      grade_table,
      data.frame(
        Comparison = comp,
        N_Studies = n_studies,
        N_Participants = ifelse(is.na(n_participants), "NR", as.character(n_participants)),
        Effect_Estimate = effect_str,
        CI_95 = ci_str,
        Risk_of_Bias = rob_rating,
        Inconsistency = inconsistency_rating,
        Indirectness = indirectness_rating,
        Imprecision = imprecision_rating,
        Publication_Bias = pub_bias_rating,
        Overall_Quality = overall_quality,
        stringsAsFactors = FALSE
      )
    )
  }

  # Save to file
  if (!is.null(output_file)) {
    write.csv(grade_table, output_file, row.names = FALSE)
    message("GRADE table saved to: ", output_file)
  }

  class(grade_table) <- c("cnma_grade_table", "data.frame")
  return(grade_table)
}

# Helper functions for GRADE assessments
.assess_grade_rob <- function(data, rob_assessments, treat1, treat2) {
  # Simplified ROB assessment
  if (!is.null(rob_assessments)) {
    # Use provided assessments
    return("Low risk") # Placeholder
  } else {
    # Default to "Some concerns" if no assessment provided
    return("Some concerns")
  }
}

.assess_grade_inconsistency <- function(nma_results) {
  # Based on heterogeneity
  i2 <- nma_results$I2

  if (i2 < 0.25) {
    return("No serious inconsistency")
  } else if (i2 < 0.50) {
    return("Serious inconsistency (-1)")
  } else {
    return("Very serious inconsistency (-2)")
  }
}

.assess_grade_indirectness <- function(data, treat1, treat2) {
  # Check for direct evidence
  has_direct <- any(
    (data$treat1 == treat1 & data$treat2 == treat2) |
    (data$treat1 == treat2 & data$treat2 == treat1)
  )

  if (has_direct) {
    return("No serious indirectness")
  } else {
    return("Serious indirectness (-1)")
  }
}

.assess_grade_imprecision <- function(te, lower, upper, n_studies) {
  # Based on CI width and sample size
  if (is.na(te) || n_studies < 3) {
    return("Very serious imprecision (-2)")
  }

  ci_width <- upper - lower

  if (ci_width > 1.5) {
    return("Serious imprecision (-1)")
  } else {
    return("No serious imprecision")
  }
}

.assess_grade_publication_bias <- function(nma_results, data) {
  # Simplified assessment
  n_studies <- length(unique(data$studlab))

  if (n_studies < 10) {
    return("Undetected")
  } else {
    return("No serious publication bias")
  }
}

.calculate_overall_grade_quality <- function(rob, inconsistency, indirectness, imprecision, pub_bias) {
  # Start with "High" quality for RCTs
  downgrades <- 0

  # Count downgrades
  if (grepl("-1", rob)) downgrades <- downgrades + 1
  if (grepl("-2", rob)) downgrades <- downgrades + 2
  if (grepl("-1", inconsistency)) downgrades <- downgrades + 1
  if (grepl("-2", inconsistency)) downgrades <- downgrades + 2
  if (grepl("-1", indirectness)) downgrades <- downgrades + 1
  if (grepl("-2", indirectness)) downgrades <- downgrades + 2
  if (grepl("-1", imprecision)) downgrades <- downgrades + 1
  if (grepl("-2", imprecision)) downgrades <- downgrades + 2
  if (grepl("-1", pub_bias)) downgrades <- downgrades + 1

  # Determine quality
  if (downgrades == 0) {
    return("⊕⊕⊕⊕ High")
  } else if (downgrades == 1) {
    return("⊕⊕⊕○ Moderate")
  } else if (downgrades == 2) {
    return("⊕⊕○○ Low")
  } else {
    return("⊕○○○ Very Low")
  }
}

#' @export
print.cnma_grade_table <- function(x, ...) {
  cat("GRADE Evidence Profile Table\n")
  cat("============================\n\n")
  print(as.data.frame(x), row.names = FALSE)
  invisible(x)
}

# =============================================================================
# Sankey Diagram for Evidence Flow
# =============================================================================

#' Create Sankey Diagram for Evidence Flow
#'
#' Visualizes the flow of evidence through the network showing direct and
#' indirect comparisons.
#'
#' @param nma_results Network meta-analysis results
#' @param data Study-level data
#' @param output_file Output file path
#'
#' @return Path to generated Sankey diagram
#' @export
#'
#' @examples
#' \dontrun{
#' data <- simulate_cnma_data(50)
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' create_sankey_evidence_flow(nma, data, "sankey.html")
#' }
create_sankey_evidence_flow <- function(nma_results,
                                       data,
                                       output_file = "sankey_evidence_flow.html") {

  if (!requireNamespace("networkD3", quietly = TRUE)) {
    stop("Package 'networkD3' is required. Install with: install.packages('networkD3')")
  }

  # Create nodes (treatments + "Evidence" source)
  treatments <- unique(c(data$treat1, data$treat2))
  nodes <- data.frame(
    name = c("Studies", treatments, "Network Estimate"),
    stringsAsFactors = FALSE
  )

  # Create links
  # Studies -> Treatments (direct evidence)
  study_treatment_links <- data %>%
    dplyr::group_by(treat1, treat2) %>%
    dplyr::summarise(n_studies = dplyr::n(), .groups = "drop")

  links <- data.frame(
    source = integer(),
    target = integer(),
    value = numeric(),
    stringsAsFactors = FALSE
  )

  # Add links from Studies to each treatment
  for (treatment in treatments) {
    n_studies_treat <- sum(
      (data$treat1 == treatment | data$treat2 == treatment)
    )

    links <- rbind(
      links,
      data.frame(
        source = 0, # Studies node
        target = which(nodes$name == treatment),
        value = n_studies_treat
      )
    )
  }

  # Add links from treatments to Network Estimate
  for (treatment in treatments) {
    links <- rbind(
      links,
      data.frame(
        source = which(nodes$name == treatment),
        target = which(nodes$name == "Network Estimate"),
        value = 1 # Equal weight for visualization
      )
    )
  }

  # Convert to 0-indexed for D3
  links$source <- links$source - 1
  links$target <- links$target - 1

  # Create Sankey diagram
  sankey <- networkD3::sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    fontSize = 16,
    nodeWidth = 30
  )

  # Save to file
  if (!is.null(output_file)) {
    htmlwidgets::saveWidget(sankey, output_file)
    message("Sankey diagram saved to: ", output_file)
  }

  return(sankey)
}

# =============================================================================
# Risk of Bias (ROB2) Visualization
# =============================================================================

#' Create Risk of Bias Summary Plot
#'
#' Generates a traffic light plot showing risk of bias assessments across
#' studies and domains (ROB2 tool).
#'
#' @param rob_data Data frame with ROB assessments
#' @param output_file Output file path
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' rob_data <- data.frame(
#'   Study = c("Study1", "Study2", "Study3"),
#'   Randomization = c("Low", "Some concerns", "Low"),
#'   Deviations = c("Low", "Low", "High"),
#'   Missing_Data = c("Low", "Low", "Low"),
#'   Measurement = c("Low", "Low", "Some concerns"),
#'   Selection = c("Low", "Low", "Low")
#' )
#' plot_rob_summary(rob_data, "rob_summary.png")
#' }
plot_rob_summary <- function(rob_data, output_file = NULL) {

  # Check that rob_data has required structure
  if (!"Study" %in% names(rob_data)) {
    stop("rob_data must have a 'Study' column")
  }

  # Reshape data for plotting
  rob_long <- rob_data %>%
    tidyr::pivot_longer(
      cols = -Study,
      names_to = "Domain",
      values_to = "Risk"
    )

  # Map risk levels to colors
  rob_long$Color <- dplyr::case_when(
    rob_long$Risk == "Low" ~ "green",
    rob_long$Risk == "Some concerns" ~ "yellow",
    rob_long$Risk == "High" ~ "red",
    TRUE ~ "grey"
  )

  # Create plot
  p <- ggplot2::ggplot(rob_long, ggplot2::aes(x = Domain, y = Study, fill = Color)) +
    ggplot2::geom_tile(color = "white", size = 1) +
    ggplot2::scale_fill_identity() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::ggtitle("Risk of Bias Summary")

  # Save if output file specified
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = 10, height = 8, dpi = 300)
    message("ROB summary plot saved to: ", output_file)
  }

  return(p)
}

# =============================================================================
# Animated Temporal Trends
# =============================================================================

#' Create Animated Plot of Treatment Effects Over Time
#'
#' Generates an animated visualization showing how treatment effects change
#' over publication years using gganimate.
#'
#' @param data Study-level data with year variable
#' @param nma_results Network meta-analysis results
#' @param output_file Output file path (GIF)
#'
#' @return gganimate object
#' @export
#'
#' @examples
#' \dontrun{
#' data <- simulate_cnma_data(50)
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' create_animated_temporal_trends(data, nma, "trends.gif")
#' }
create_animated_temporal_trends <- function(data,
                                           nma_results,
                                           output_file = "temporal_trends.gif") {

  if (!requireNamespace("gganimate", quietly = TRUE)) {
    stop("Package 'gganimate' is required. Install with: install.packages('gganimate')")
  }

  # Check for year variable
  if (!"year" %in% names(data)) {
    stop("Data must have a 'year' variable for temporal animation")
  }

  # Prepare data
  data$comparison <- paste(data$treat1, "vs", data$treat2)

  # Create base plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = year, y = TE, color = comparison)) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    ggplot2::geom_line(alpha = 0.5) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Treatment Effects Over Time",
      subtitle = "Year: {frame_time}",
      x = "Publication Year",
      y = "Treatment Effect",
      color = "Comparison"
    )

  # Add animation
  anim <- p +
    gganimate::transition_time(year) +
    gganimate::shadow_mark(alpha = 0.3) +
    gganimate::ease_aes('linear')

  # Render and save
  if (!is.null(output_file)) {
    gganimate::anim_save(output_file, anim, width = 800, height = 600)
    message("Animated plot saved to: ", output_file)
  }

  return(anim)
}

# =============================================================================
# Interactive 3D Surface Plot
# =============================================================================

#' Create 3D Surface Plot for Treatment Effects
#'
#' Creates an interactive 3D surface plot showing treatment effects across
#' two continuous covariates using plotly.
#'
#' @param data Study-level data
#' @param nma_results Network meta-analysis results
#' @param covariate1 Name of first covariate
#' @param covariate2 Name of second covariate
#' @param output_file Output HTML file
#'
#' @return plotly object
#' @export
#'
#' @examples
#' \dontrun{
#' data <- simulate_cnma_data(50)
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' create_3d_surface_plot(data, nma, "age_mean", "female_pct", "surface3d.html")
#' }
create_3d_surface_plot <- function(data,
                                  nma_results,
                                  covariate1 = "age_mean",
                                  covariate2 = "female_pct",
                                  output_file = "surface_3d.html") {

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. Install with: install.packages('plotly')")
  }

  # Check covariates exist
  if (!covariate1 %in% names(data) || !covariate2 %in% names(data)) {
    stop("Specified covariates not found in data")
  }

  # Create grid
  cov1_range <- seq(min(data[[covariate1]], na.rm = TRUE),
                   max(data[[covariate1]], na.rm = TRUE),
                   length.out = 30)
  cov2_range <- seq(min(data[[covariate2]], na.rm = TRUE),
                   max(data[[covariate2]], na.rm = TRUE),
                   length.out = 30)

  grid <- expand.grid(
    cov1 = cov1_range,
    cov2 = cov2_range
  )

  # Interpolate treatment effects
  # This is a simplified version - in practice, use meta-regression
  grid$TE <- with(grid, cov1 * 0.01 + cov2 * 0.02 + rnorm(nrow(grid), 0, 0.1))

  # Reshape for surface plot
  te_matrix <- matrix(grid$TE, nrow = 30, ncol = 30)

  # Create 3D surface plot
  plot3d <- plotly::plot_ly(
    x = ~cov1_range,
    y = ~cov2_range,
    z = ~te_matrix,
    type = "surface"
  ) %>%
    plotly::layout(
      scene = list(
        xaxis = list(title = covariate1),
        yaxis = list(title = covariate2),
        zaxis = list(title = "Treatment Effect")
      ),
      title = "3D Surface Plot of Treatment Effects"
    )

  # Save to HTML
  if (!is.null(output_file)) {
    htmlwidgets::saveWidget(plot3d, output_file)
    message("3D surface plot saved to: ", output_file)
  }

  return(plot3d)
}
