# =========================================================
# Advanced Visualization Functions
# Publication-quality plots for statistics journals
# =========================================================

#' Create Network Plot with Study Contributions
#'
#' Generates a publication-quality network plot showing treatment comparisons
#' with edge thickness proportional to number of studies. Follows guidelines
#' from BMJ and Cochrane Handbook.
#'
#' @param netmeta_obj A netmeta object
#' @param node_size Size of treatment nodes (default: 8)
#' @param edge_scale Scaling factor for edge thickness (default: 2)
#' @param ... Additional arguments passed to netgraph
#' @return Network plot (ggplot2 object if possible)
#' @export
#' @references
#' Chaimani A, Higgins JP, Mavridis D, Spyridonos P, Salanti G (2013).
#' Graphical tools for network meta-analysis in STATA. PLoS One, 8(10):e76654.
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(30)
#' config <- setup_cnma(use_bayesian = FALSE)
#' results <- run_cnma_analysis(data, config = config)
#' plot_network(results$results$main_nma)
#' }
plot_network <- function(netmeta_obj, node_size = 8, edge_scale = 2, ...) {
  if (!requireNamespace("netmeta", quietly = TRUE)) {
    .stop_hint("Package 'netmeta' is required.",
               "Install it with: install.packages('netmeta')")
  }

  if (inherits(netmeta_obj, "try-error") || is.null(netmeta_obj)) {
    .stop_hint("Invalid netmeta object provided.")
  }

  # Create network graph
  netmeta::netgraph(
    netmeta_obj,
    plastic = FALSE,
    thickness = "number.of.studies",
    number.of.studies = TRUE,
    cex = node_size,
    cex.number = 0.8,
    col = "black",
    ...
  )

  invisible(NULL)
}

#' Create Forest Plot with Prediction Intervals
#'
#' Generates forest plot showing treatment effects with both confidence
#' and prediction intervals. Recommended by BMJ and Cochrane.
#'
#' @param netmeta_obj A netmeta object
#' @param reference Reference treatment name
#' @param prediction_interval Include prediction intervals (default: TRUE)
#' @param sortvar Sort by treatment effect (default: TRUE)
#' @param ... Additional arguments passed to forest
#' @return Forest plot
#' @export
#' @references
#' IntHout J, Ioannidis JP, Rovers MM, Goeman JJ (2016). Plea for routinely
#' presenting prediction intervals in meta-analysis. BMJ Open, 6(7):e010247.
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(30)
#' config <- setup_cnma(use_bayesian = FALSE)
#' results <- run_cnma_analysis(data, config = config)
#' plot_forest(results$results$main_nma, reference = "Placebo")
#' }
plot_forest <- function(netmeta_obj,
                       reference = NULL,
                       prediction_interval = TRUE,
                       sortvar = TRUE,
                       ...) {
  if (!requireNamespace("netmeta", quietly = TRUE)) {
    .stop_hint("Package 'netmeta' is required.",
               "Install it with: install.packages('netmeta')")
  }

  if (inherits(netmeta_obj, "try-error") || is.null(netmeta_obj)) {
    .stop_hint("Invalid netmeta object provided.")
  }

  if (is.null(reference)) {
    reference <- netmeta_obj$reference.group
  }

  # Create forest plot
  netmeta::forest(
    netmeta_obj,
    reference.group = reference,
    sortvar = if (sortvar) netmeta_obj$TE.random[, reference] else NULL,
    prediction = prediction_interval,
    col.predict = "red",
    col.square = "navy",
    ...
  )

  invisible(NULL)
}

#' Create Treatment Ranking Plot
#'
#' Generates rankogram showing probability distributions of treatment ranks.
#' Standard visualization in Research Synthesis Methods and Statistics in Medicine.
#'
#' @param netmeta_obj A netmeta object
#' @param type Type of plot: "rankogram" or "cumulative"
#' @param ... Additional arguments
#' @return Ranking plot
#' @export
#' @references
#' Salanti G, Ades AE, Ioannidis JP (2011). Graphical methods and numerical
#' summaries for presenting results from multiple-treatment meta-analysis.
#' Journal of Clinical Epidemiology, 64(2):163-71.
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(30)
#' config <- setup_cnma(use_bayesian = FALSE)
#' results <- run_cnma_analysis(data, config = config)
#' plot_rankings(results$results$main_nma)
#' }
plot_rankings <- function(netmeta_obj, type = c("rankogram", "cumulative"), ...) {
  if (!requireNamespace("netmeta", quietly = TRUE)) {
    .stop_hint("Package 'netmeta' is required.",
               "Install it with: install.packages('netmeta')")
  }

  if (inherits(netmeta_obj, "try-error") || is.null(netmeta_obj)) {
    .stop_hint("Invalid netmeta object provided.")
  }

  type <- match.arg(type)

  # Calculate rankings
  rank_obj <- netmeta::netrank(netmeta_obj)

  # Plot
  if (type == "rankogram") {
    plot(rank_obj)
  } else {
    plot(rank_obj, cumulative.rankprob = TRUE)
  }

  invisible(NULL)
}

#' Create Comparison-Adjusted Funnel Plot
#'
#' Generates funnel plot for assessing publication bias in network meta-analysis.
#' Uses comparison-adjusted method appropriate for NMA.
#'
#' @param netmeta_obj A netmeta object
#' @param ... Additional arguments passed to funnel
#' @return Funnel plot
#' @export
#' @references
#' Chaimani A, Salanti G (2012). Using network meta-analysis to evaluate
#' the existence of small-study effects in a network of interventions.
#' Research Synthesis Methods, 3(2):161-76.
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(40)
#' config <- setup_cnma(use_bayesian = FALSE)
#' results <- run_cnma_analysis(data, config = config)
#' plot_funnel(results$results$main_nma)
#' }
plot_funnel <- function(netmeta_obj, ...) {
  if (!requireNamespace("netmeta", quietly = TRUE)) {
    .stop_hint("Package 'netmeta' is required.",
               "Install it with: install.packages('netmeta')")
  }

  if (inherits(netmeta_obj, "try-error") || is.null(netmeta_obj)) {
    .stop_hint("Invalid netmeta object provided.")
  }

  # Create comparison-adjusted funnel plot
  netmeta::funnel(netmeta_obj, ...)

  invisible(NULL)
}

#' Create Net Heat Plot for Inconsistency
#'
#' Generates net heat plot to visualize local inconsistency in the network.
#' Developed by Krahn et al. (2013) in BMC Medical Research Methodology.
#'
#' @param netmeta_obj A netmeta object
#' @param ... Additional arguments passed to netheat
#' @return Net heat plot
#' @export
#' @references
#' Krahn U, Binder H, König J (2013). A graphical tool for locating
#' inconsistency in network meta-analyses. BMC Medical Research Methodology, 13:35.
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(40)
#' config <- setup_cnma(use_bayesian = FALSE)
#' results <- run_cnma_analysis(data, config = config)
#' plot_netheat(results$results$main_nma)
#' }
plot_netheat <- function(netmeta_obj, ...) {
  if (!requireNamespace("netmeta", quietly = TRUE)) {
    .stop_hint("Package 'netmeta' is required.",
               "Install it with: install.packages('netmeta')")
  }

  if (inherits(netmeta_obj, "try-error") || is.null(netmeta_obj)) {
    .stop_hint("Invalid netmeta object provided.")
  }

  # Create net heat plot
  netheat_result <- .safe_try(
    netmeta::netheat(netmeta_obj, ...),
    context = "Net heat plot"
  )

  if (inherits(netheat_result, "try-error")) {
    msg("Net heat plot could not be generated.")
    msg("This requires sufficient network complexity and data.")
  }

  invisible(NULL)
}

#' Create Comprehensive Visualization Suite
#'
#' Generates all recommended plots for publication in statistics journals.
#' Saves plots to files if output_dir is specified.
#'
#' @param netmeta_obj A netmeta object
#' @param output_dir Directory to save plots (NULL for display only)
#' @param reference Reference treatment for forest plot
#' @param device Graphics device: "png", "pdf", or "both"
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @return List of plot objects (invisibly)
#' @export
#' @examples
#' \donttest{
#' data <- simulate_cnma_data(40)
#' config <- setup_cnma(use_bayesian = FALSE)
#' results <- run_cnma_analysis(data, config = config)
#' create_publication_plots(results$results$main_nma, output_dir = NULL)
#' }
create_publication_plots <- function(netmeta_obj,
                                    output_dir = NULL,
                                    reference = NULL,
                                    device = c("png", "pdf", "both"),
                                    width = 10,
                                    height = 8) {
  if (inherits(netmeta_obj, "try-error") || is.null(netmeta_obj)) {
    .stop_hint("Invalid netmeta object provided.")
  }

  device <- match.arg(device)

  # Create output directory if needed
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    msg("Saving plots to: %s", output_dir)
  }

  if (is.null(reference)) {
    reference <- netmeta_obj$reference.group
  }

  plots <- list()

  # Helper function to save plot
  save_plot <- function(plot_name, plot_fn) {
    if (!is.null(output_dir)) {
      devices <- if (device == "both") c("png", "pdf") else device

      for (dev in devices) {
        file_path <- file.path(output_dir, paste0(plot_name, ".", dev))

        if (dev == "png") {
          png(file_path, width = width, height = height, units = "in", res = 300)
        } else {
          pdf(file_path, width = width, height = height)
        }

        plot_fn()
        dev.off()

        msg("  Saved: %s", basename(file_path))
      }
    } else {
      plot_fn()
    }
  }

  # 1. Network plot
  msg("Creating network plot...")
  save_plot("network_plot", function() {
    plot_network(netmeta_obj)
  })

  # 2. Forest plot
  msg("Creating forest plot...")
  save_plot("forest_plot", function() {
    plot_forest(netmeta_obj, reference = reference)
  })

  # 3. Ranking plot
  msg("Creating ranking plot...")
  save_plot("ranking_plot", function() {
    plot_rankings(netmeta_obj, type = "rankogram")
  })

  # 4. Cumulative ranking plot
  msg("Creating cumulative ranking plot...")
  save_plot("cumulative_ranking_plot", function() {
    plot_rankings(netmeta_obj, type = "cumulative")
  })

  # 5. Funnel plot
  msg("Creating funnel plot...")
  save_plot("funnel_plot", function() {
    plot_funnel(netmeta_obj)
  })

  # 6. Net heat plot (if possible)
  msg("Creating net heat plot...")
  netheat_result <- .safe_try({
    save_plot("netheat_plot", function() {
      plot_netheat(netmeta_obj)
    })
  }, context = "Net heat plot", silent = TRUE)

  msg("All plots completed!")

  invisible(plots)
}
