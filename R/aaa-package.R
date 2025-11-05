# =========================================================
# CNMA: Comprehensive Network Meta-Analysis
# Package Documentation and Lifecycle Hooks
# =========================================================

#' CNMA: Comprehensive Network Meta-Analysis
#'
#' A comprehensive toolkit for conducting network meta-analysis with both
#' frequentist and Bayesian approaches, including transportability analysis,
#' GRADE weighting, bias assessment, and extensive diagnostics.
#'
#' @docType package
#' @name cnma-package
#' @aliases cnma
#' @author CNMA Development Team
#' @import netmeta ggplot2 dplyr tidyr purrr tibble
#' @importFrom stats lm rnorm runif plogis qlogis cov median quantile sd na.omit
#' @importFrom stats AIC coef predict ks.test binomial
#' @importFrom utils capture.output write.csv installed.packages sessionInfo
#' @importFrom grDevices dev.off pdf png svg
#' @importFrom graphics plot
NULL

# ---- Package Environment ----
.cnma_env <- new.env(parent = emptyenv())

#' Initialize CNMA Environment
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  .cnma_env$cache <- new.env(parent = emptyenv())
  .cnma_env$parallel_enabled <- FALSE
}

#' Clean up CNMA Environment
#' @keywords internal
.onUnload <- function(libpath) {
  if (exists("parallel_enabled", envir = .cnma_env)) {
    if (.cnma_env$parallel_enabled && requireNamespace("future", quietly = TRUE)) {
      future::plan(future::sequential)
    }
  }
}

#' @keywords internal
"_PACKAGE"
