# =========================================================
# Parallel Processing Functions
# =========================================================

#' Enable parallel processing
#'
#' @param strategy Parallel strategy: "sequential", "multisession", or "multicore"
#' @param workers Number of workers (NULL for default)
#' @return Logical indicating success
#' @export
#' @examples
#' \donttest{
#' cnma_parallel_on("sequential")
#' }
cnma_parallel_on <- function(strategy = c("sequential", "multisession", "multicore"),
                            workers = NULL) {
  strategy <- match.arg(strategy)

  if (is_cran() || !has_pkg("future")) {
    .cnma_env$parallel_enabled <- FALSE
    return(invisible(FALSE))
  }

  future <- getNamespace("future")
  if (!is.null(workers)) {
    future$plan(strategy, workers = workers)
  } else {
    future$plan(strategy)
  }

  .cnma_env$parallel_enabled <- TRUE
  invisible(TRUE)
}

#' Disable parallel processing
#' @export
#' @examples
#' cnma_parallel_off()
cnma_parallel_off <- function() {
  if (has_pkg("future")) {
    future <- getNamespace("future")
    future$plan(future$sequential)
  }
  .cnma_env$parallel_enabled <- FALSE
}

#' Parallel apply
#' @param x List
#' @param FUN Function to apply
#' @return List of results
#' @keywords internal
.papply <- function(x, FUN) {
  if (.cnma_env$parallel_enabled && has_pkg("future.apply")) {
    future.apply <- getNamespace("future.apply")
    future.apply$future_lapply(x, FUN)
  } else {
    lapply(x, FUN)
  }
}
