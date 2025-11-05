# =========================================================
# S3 Methods
# =========================================================

#' Print CNMA Results
#'
#' @param x CNMA object
#' @param ... Additional arguments
#' @export
print.cnma <- function(x, ...) {
  cat("<CNMA Results>\n")
  cat(sprintf("  Summary measure: %s\n", x$config$sm))
  cat(sprintf("  Reference: %s\n", x$ref_treatment))

  if (!is.null(x$results$main_nma) &&
      !inherits(x$results$main_nma, "try-error")) {
    cat(sprintf("  Tau: %.4f\n", x$results$main_nma$tau))
    cat(sprintf("  I²: %.1f%%\n", x$results$main_nma$I2.random * 100))
  }

  invisible(x)
}

#' Summary of CNMA Results
#'
#' @param object CNMA object
#' @param ... Additional arguments
#' @export
summary.cnma <- function(object, ...) {
  cat("CNMA Analysis Summary\n")
  cat("=====================\n\n")

  print(object)

  if (!is.null(object$results$main_nma) &&
      !inherits(object$results$main_nma, "try-error")) {
    cat("\nNetwork characteristics:\n")

    n_studies <- length(unique(object$data$studlab))
    n_treatments <- length(unique(c(object$data$treat1, object$data$treat2)))
    n_comparisons <- nrow(object$data)

    cat(sprintf("  Studies: %d\n", n_studies))
    cat(sprintf("  Treatments: %d\n", n_treatments))
    cat(sprintf("  Comparisons: %d\n", n_comparisons))

    # Additional network info
    cat("\nHeterogeneity:\n")
    cat(sprintf("  Tau²: %.4f\n", object$results$main_nma$tau^2))
    cat(sprintf("  Tau: %.4f\n", object$results$main_nma$tau))
    cat(sprintf("  I²: %.1f%%\n", object$results$main_nma$I2.random * 100))
  }

  invisible(object)
}
