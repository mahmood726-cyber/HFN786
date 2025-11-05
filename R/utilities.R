# =========================================================
# Utility Functions
# =========================================================

#' Null coalescing operator
#' @param a First value
#' @param b Alternative value if a is NULL
#' @return a if not NULL, otherwise b
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Safe value clipping
#' @param x Numeric vector
#' @param lo Lower bound
#' @param hi Upper bound
#' @return Clipped values
#' @keywords internal
safe_clip <- function(x, lo, hi) pmin(hi, pmax(lo, x))

#' Replace NA with value
#' @param x Vector
#' @param val Replacement value
#' @return Vector with NAs replaced
#' @keywords internal
vcoalesce <- function(x, val = 0) {
  x[is.na(x)] <- val
  x
}

#' Check if package is available
#' @param p Package name
#' @return Logical
#' @keywords internal
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

#' Print message with timestamp
#' @param ... Message components
#' @keywords internal
msg <- function(...) {
  if (!getOption("cnma.quiet", FALSE)) {
    cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
        sprintf(...), "\n")
  }
}

#' Check if running on CRAN
#' @return Logical
#' @keywords internal
is_cran <- function() {
  !identical(Sys.getenv("NOT_CRAN"), "true") &&
    (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")) ||
     nzchar(Sys.getenv("_R_CHECK_SIZE_OF_TARBALL_")))
}

#' Check if JAGS is available
#' @return Logical
#' @keywords internal
has_jags <- function() {
  if (is_cran()) return(FALSE)
  out1 <- suppressWarnings(Sys.which("jags"))
  out2 <- suppressWarnings(Sys.which("JAGS"))
  nzchar(out1) || nzchar(out2)
}

#' Stop with hint
#' @param msg0 Error message
#' @param hint Optional hint
#' @keywords internal
.stop_hint <- function(msg0, hint = NULL) {
  if (!is.null(hint)) msg0 <- paste0(msg0, "\nHint: ", hint)
  stop(msg0, call. = FALSE)
}

#' Safe try wrapper
#' @param expr Expression to evaluate
#' @param context Context string
#' @param silent Silent mode
#' @return Result or try-error
#' @keywords internal
.safe_try <- function(expr, context = "", silent = TRUE) {
  out <- try(expr, silent = silent)
  if (inherits(out, "try-error")) {
    msg("ERROR in %s: %s", context,
        as.character(attr(out, "condition")$message %||% out))
  }
  out
}

# ---- Cache Functions ----

#' Generate cache key
#' @param tag Tag identifier
#' @param ... Additional parameters
#' @return Cache key string
#' @keywords internal
cache_key <- function(tag, ...) {
  if (!has_pkg("digest")) {
    # Improved fallback: use timestamp and more random characters
    return(paste0(tag, "_",
                  format(Sys.time(), "%Y%m%d_%H%M%S_"),
                  paste(sample(c(letters, LETTERS, 0:9), 12, replace = TRUE),
                        collapse = "")))
  }
  digest::digest(list(tag = tag, ...), algo = "xxhash64")
}

#' Memoize function results
#' @param key Cache key
#' @param expr Expression to evaluate
#' @param enable_cache Enable caching
#' @return Cached or computed result
#' @keywords internal
memoize <- function(key, expr, enable_cache = FALSE) {
  if (!enable_cache) return(eval.parent(substitute(expr)))

  cache <- .cnma_env$cache
  if (exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache))
  }

  val <- eval.parent(substitute(expr))
  assign(key, val, envir = cache)
  val
}
