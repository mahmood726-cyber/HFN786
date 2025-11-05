# =========================================================
# Interactive Visualizations
# Modern interactive plots with plotly for exploration
# =========================================================

#' Create Interactive Network Plot
#'
#' Generates an interactive network plot using plotly, allowing users to
#' hover over nodes and edges to see details, zoom, and pan.
#'
#' @param netmeta_obj A netmeta object
#' @param node_size_var Variable to scale node size ("n_studies", "constant")
#' @param edge_width_var Variable to scale edge width ("n_studies", "precision")
#' @param layout Layout algorithm ("fr", "circle", "star")
#' @param title Plot title
#' @return plotly object (interactive plot)
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   data <- simulate_cnma_data(30)
#'   config <- setup_cnma(use_bayesian = FALSE)
#'   results <- run_cnma_analysis(data, config = config)
#'   plot_interactive_network(results$results$main_nma)
#' }
#' }
plot_interactive_network <- function(netmeta_obj,
                                     node_size_var = "n_studies",
                                     edge_width_var = "n_studies",
                                     layout = "fr",
                                     title = "Interactive Network Plot") {

  if (!requireNamespace("plotly", quietly = TRUE)) {
    msg("Package 'plotly' is required for interactive plots.")
    msg("Install it with: install.packages('plotly')")
    return(invisible(NULL))
  }

  if (!requireNamespace("igraph", quietly = TRUE)) {
    msg("Package 'igraph' is required for network layouts.")
    msg("Install it with: install.packages('igraph')")
    return(invisible(NULL))
  }

  # Extract network structure
  treatments <- netmeta_obj$trts
  comparisons <- data.frame(
    treat1 = netmeta_obj$treat1,
    treat2 = netmeta_obj$treat2,
    studlab = netmeta_obj$studlab,
    stringsAsFactors = FALSE
  )

  # Count studies per comparison
  edge_counts <- comparisons %>%
    group_by(treat1, treat2) %>%
    summarise(n_studies = n(), .groups = "drop")

  # Create igraph object
  edges_for_graph <- edge_counts %>%
    select(treat1, treat2, n_studies)

  g <- igraph::graph_from_data_frame(edges_for_graph, directed = FALSE,
                                     vertices = treatments)

  # Layout
  layout_mat <- switch(layout,
                      "fr" = igraph::layout_with_fr(g),
                      "circle" = igraph::layout_in_circle(g),
                      "star" = igraph::layout_as_star(g),
                      igraph::layout_with_fr(g))

  # Node coordinates
  node_coords <- data.frame(
    treatment = treatments,
    x = layout_mat[, 1],
    y = layout_mat[, 2],
    degree = igraph::degree(g),
    stringsAsFactors = FALSE
  )

  # Edge coordinates
  edge_list <- igraph::get.edgelist(g)
  edge_coords <- list()
  for(i in 1:nrow(edge_list)) {
    from_node <- node_coords[node_coords$treatment == edge_list[i, 1], ]
    to_node <- node_coords[node_coords$treatment == edge_list[i, 2], ]

    n_st <- edge_counts$n_studies[edge_counts$treat1 == edge_list[i, 1] &
                                   edge_counts$treat2 == edge_list[i, 2]]
    if(length(n_st) == 0) {
      n_st <- edge_counts$n_studies[edge_counts$treat1 == edge_list[i, 2] &
                                     edge_counts$treat2 == edge_list[i, 1]]
    }

    edge_coords[[i]] <- data.frame(
      x = c(from_node$x, to_node$x, NA),
      y = c(from_node$y, to_node$y, NA),
      comparison = paste(edge_list[i, 1], "vs", edge_list[i, 2]),
      n_studies = n_st,
      stringsAsFactors = FALSE
    )
  }
  edge_df <- do.call(rbind, edge_coords)

  # Create plotly plot
  p <- plotly::plot_ly()

  # Add edges
  p <- p %>%
    plotly::add_trace(
      data = edge_df,
      x = ~x, y = ~y,
      type = "scatter",
      mode = "lines",
      line = list(color = "lightgray", width = ~n_studies * 2),
      hoverinfo = "text",
      text = ~paste0(comparison, "<br>Studies: ", n_studies),
      showlegend = FALSE
    )

  # Add nodes
  p <- p %>%
    plotly::add_trace(
      data = node_coords,
      x = ~x, y = ~y,
      type = "scatter",
      mode = "markers+text",
      marker = list(
        size = ~degree * 8 + 15,
        color = "steelblue",
        line = list(color = "white", width = 2)
      ),
      text = ~treatment,
      textposition = "top center",
      hoverinfo = "text",
      hovertext = ~paste0("Treatment: ", treatment,
                          "<br>Connections: ", degree),
      showlegend = FALSE
    )

  # Layout
  p <- p %>%
    plotly::layout(
      title = title,
      xaxis = list(showgrid = FALSE, showticklabels = FALSE, title = ""),
      yaxis = list(showgrid = FALSE, showticklabels = FALSE, title = ""),
      hovermode = "closest",
      plot_bgcolor = "white"
    )

  return(p)
}

#' Create Interactive Forest Plot
#'
#' Generates an interactive forest plot showing treatment effects with
#' confidence and prediction intervals.
#'
#' @param netmeta_obj A netmeta object
#' @param reference Reference treatment
#' @param prediction_interval Include prediction intervals
#' @param title Plot title
#' @return plotly object
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   data <- simulate_cnma_data(25)
#'   config <- setup_cnma(use_bayesian = FALSE)
#'   results <- run_cnma_analysis(data, config = config)
#'   plot_interactive_forest(results$results$main_nma)
#' }
#' }
plot_interactive_forest <- function(netmeta_obj,
                                   reference = NULL,
                                   prediction_interval = TRUE,
                                   title = "Interactive Forest Plot") {

  if (!requireNamespace("plotly", quietly = TRUE)) {
    msg("Package 'plotly' required. Install with: install.packages('plotly')")
    return(invisible(NULL))
  }

  if (is.null(reference)) {
    reference <- netmeta_obj$reference.group
  }

  # Extract effects
  treatments <- netmeta_obj$trts
  treatments <- treatments[treatments != reference]

  effects <- netmeta_obj$TE.random[reference, treatments]
  lower_ci <- netmeta_obj$lower.random[reference, treatments]
  upper_ci <- netmeta_obj$upper.random[reference, treatments]

  # Prediction intervals if requested
  if (prediction_interval) {
    pred_int <- calculate_prediction_intervals(netmeta_obj, level = 0.95)
    pred_int_ref <- pred_int[pred_int$treatment_1 == reference |
                             pred_int$treatment_2 == reference, ]
  }

  # Create data frame
  forest_data <- data.frame(
    treatment = treatments,
    effect = as.numeric(effects),
    ci_lower = as.numeric(lower_ci),
    ci_upper = as.numeric(upper_ci),
    stringsAsFactors = FALSE
  )

  # Sort by effect
  forest_data <- forest_data[order(forest_data$effect), ]
  forest_data$y <- nrow(forest_data):1

  # Create plot
  p <- plotly::plot_ly()

  # Add reference line
  p <- p %>%
    plotly::add_trace(
      x = c(0, 0),
      y = c(0, nrow(forest_data) + 1),
      type = "scatter",
      mode = "lines",
      line = list(color = "red", dash = "dash", width = 1),
      showlegend = FALSE,
      hoverinfo = "skip"
    )

  # Add CI bars
  for(i in 1:nrow(forest_data)) {
    p <- p %>%
      plotly::add_trace(
        x = c(forest_data$ci_lower[i], forest_data$ci_upper[i]),
        y = c(forest_data$y[i], forest_data$y[i]),
        type = "scatter",
        mode = "lines",
        line = list(color = "black", width = 2),
        showlegend = FALSE,
        hoverinfo = "skip"
      )
  }

  # Add point estimates
  p <- p %>%
    plotly::add_trace(
      data = forest_data,
      x = ~effect,
      y = ~y,
      type = "scatter",
      mode = "markers",
      marker = list(
        size = 10,
        color = "navy",
        symbol = "square",
        line = list(color = "white", width = 1)
      ),
      text = ~paste0(treatment, "<br>Effect: ", round(effect, 3),
                     "<br>95% CI: [", round(ci_lower, 3), ", ",
                     round(ci_upper, 3), "]"),
      hoverinfo = "text",
      showlegend = FALSE
    )

  # Layout
  p <- p %>%
    plotly::layout(
      title = title,
      xaxis = list(title = paste("Effect vs", reference), zeroline = FALSE),
      yaxis = list(
        title = "",
        tickvals = forest_data$y,
        ticktext = forest_data$treatment
      ),
      hovermode = "closest"
    )

  return(p)
}

#' Create Interactive Ranking Plot
#'
#' Generates an interactive visualization of treatment rankings with
#' P-scores displayed on hover.
#'
#' @param netmeta_obj A netmeta object
#' @param type Type of plot ("bar", "lollipop")
#' @param title Plot title
#' @return plotly object
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   data <- simulate_cnma_data(25)
#'   config <- setup_cnma(use_bayesian = FALSE)
#'   results <- run_cnma_analysis(data, config = config)
#'   plot_interactive_rankings(results$results$main_nma)
#' }
#' }
plot_interactive_rankings <- function(netmeta_obj,
                                     type = c("bar", "lollipop"),
                                     title = "Interactive Treatment Rankings") {

  if (!requireNamespace("plotly", quietly = TRUE)) {
    msg("Package 'plotly' required. Install with: install.packages('plotly')")
    return(invisible(NULL))
  }

  type <- match.arg(type)

  # Calculate rankings
  rankings <- calculate_rankings(netmeta_obj)
  rankings <- rankings[order(-rankings$p_score), ]
  rankings$rank_order <- 1:nrow(rankings)

  # Create plot based on type
  if (type == "bar") {
    p <- plotly::plot_ly(
      data = rankings,
      x = ~reorder(treatment, p_score),
      y = ~p_score,
      type = "bar",
      marker = list(
        color = ~p_score,
        colorscale = list(c(0, "lightblue"), c(1, "darkblue")),
        showscale = TRUE,
        colorbar = list(title = "P-score")
      ),
      text = ~paste0("Treatment: ", treatment,
                     "<br>P-score: ", round(p_score, 3),
                     "<br>Rank: ", rank,
                     "<br>", interpretation),
      hoverinfo = "text"
    )

    p <- p %>%
      plotly::layout(
        title = title,
        xaxis = list(title = "Treatment"),
        yaxis = list(title = "P-score", range = c(0, 1)),
        showlegend = FALSE
      )

  } else {  # lollipop
    p <- plotly::plot_ly(data = rankings)

    # Add stems
    for(i in 1:nrow(rankings)) {
      p <- p %>%
        plotly::add_trace(
          x = c(rankings$p_score[i], rankings$p_score[i]),
          y = c(0, rankings$rank_order[i]),
          type = "scatter",
          mode = "lines",
          line = list(color = "gray", width = 2),
          showlegend = FALSE,
          hoverinfo = "skip"
        )
    }

    # Add points
    p <- p %>%
      plotly::add_trace(
        data = rankings,
        x = ~p_score,
        y = ~rank_order,
        type = "scatter",
        mode = "markers+text",
        marker = list(
          size = 15,
          color = ~p_score,
          colorscale = list(c(0, "lightcoral"), c(1, "darkgreen")),
          showscale = TRUE,
          colorbar = list(title = "P-score"),
          line = list(color = "white", width = 2)
        ),
        text = ~treatment,
        textposition = "right",
        hovertext = ~paste0("Treatment: ", treatment,
                           "<br>P-score: ", round(p_score, 3),
                           "<br>Rank: ", rank),
        hoverinfo = "text",
        showlegend = FALSE
      )

    p <- p %>%
      plotly::layout(
        title = title,
        xaxis = list(title = "P-score", range = c(0, 1)),
        yaxis = list(title = "Rank", tickvals = rankings$rank_order,
                    ticktext = rankings$rank)
      )
  }

  return(p)
}

#' Create Interactive Dashboard
#'
#' Generates a complete interactive dashboard with multiple linked plots
#' for comprehensive exploration of NMA results.
#'
#' @param cnma_results A cnma results object
#' @param save_html Save to HTML file
#' @param filename HTML filename if save_html = TRUE
#' @return List of plotly objects or HTML widget
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   data <- simulate_cnma_data(30)
#'   config <- setup_cnma(use_bayesian = FALSE)
#'   results <- run_cnma_analysis(data, config = config)
#'   dashboard <- create_interactive_dashboard(results)
#' }
#' }
create_interactive_dashboard <- function(cnma_results,
                                        save_html = FALSE,
                                        filename = "nma_dashboard.html") {

  if (!inherits(cnma_results, "cnma")) {
    .stop_hint("Input must be a cnma results object.")
  }

  nma <- cnma_results$results$main_nma

  msg("Creating interactive dashboard...")

  # Create individual plots
  plots <- list(
    network = plot_interactive_network(nma, title = "Network Structure"),
    forest = plot_interactive_forest(nma, title = "Treatment Effects"),
    rankings = plot_interactive_rankings(nma, type = "bar", title = "Treatment Rankings")
  )

  if (save_html) {
    if (!requireNamespace("htmltools", quietly = TRUE)) {
      msg("Package 'htmltools' required for HTML saving.")
      msg("Install with: install.packages('htmltools')")
    } else {
      msg("Saving dashboard to: %s", filename)
      # Create HTML with all plots
      html_content <- htmltools::tags$html(
        htmltools::tags$head(
          htmltools::tags$title("CNMA Interactive Dashboard")
        ),
        htmltools::tags$body(
          htmltools::tags$h1("Network Meta-Analysis Dashboard"),
          htmltools::tags$div(plots$network),
          htmltools::tags$div(plots$forest),
          htmltools::tags$div(plots$rankings)
        )
      )
      htmltools::save_html(html_content, file = filename)
    }
  }

  return(plots)
}
