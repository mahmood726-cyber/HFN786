# =========================================================
# Demo Data Generation
# =========================================================

#' Simulate CNMA Data
#'
#' Generates simulated network meta-analysis data for testing and examples.
#'
#' @param n_studies Number of studies to simulate
#' @param seed Random seed for reproducibility
#' @return Data frame with simulated NMA data
#' @export
#' @examples
#' data <- simulate_cnma_data(20, seed = 123)
#' head(data)
simulate_cnma_data <- function(n_studies = 40, seed = 42) {
  # Input validation
  if (!is.numeric(n_studies) || n_studies < 1) {
    .stop_hint("n_studies must be a positive integer.",
               "Provide a value >= 1.")
  }

  set.seed(seed)

  tr <- c("Placebo", "DrugA", "DrugB", "DrugC", "DrugD")
  true <- c(Placebo = 0, DrugA = log(0.85), DrugB = log(0.75),
           DrugC = log(0.90), DrugD = log(0.78))

  designs <- list(
    c("Placebo", "DrugA"),
    c("Placebo", "DrugB"),
    c("Placebo", "DrugC"),
    c("DrugA", "DrugB"),
    c("DrugA", "DrugC"),
    c("DrugB", "DrugD"),
    c("Placebo", "DrugA", "DrugB"),
    c("Placebo", "DrugC", "DrugD")
  )

  L <- sample(designs, n_studies, replace = TRUE)

  pw <- purrr::map_dfr(seq_along(L), function(i) {
    s <- sprintf("Study_%02d", i)
    prs <- combn(L[[i]], 2, simplify = FALSE)

    purrr::map_dfr(prs, function(p) {
      t1 <- p[1]
      t2 <- p[2]
      true_diff <- true[t1] - true[t2]
      re <- rnorm(1, 0, 0.10)
      se <- runif(1, 0.08, 0.25)
      TE <- rnorm(1, true_diff + re, se)

      tibble(
        studlab = s,
        treat1 = t1,
        treat2 = t2,
        TE = TE,
        seTE = se
      )
    })
  })

  info <- pw %>%
    distinct(studlab) %>%
    mutate(
      year = sample(2010:2025, n(), TRUE),
      is_rct = rbinom(n(), 1, 0.8),
      study_design = factor(ifelse(is_rct == 1, "RCT", "Non-RCT")),
      grade = factor(
        sample(c("High", "Moderate", "Low", "Very low"), n(), TRUE,
               prob = c(0.4, 0.4, 0.15, 0.05)),
        levels = c("High", "Moderate", "Low", "Very low")
      ),
      rob2 = sample(c("low", "some concerns", "high"), n(), TRUE,
                   prob = c(0.5, 0.35, 0.15)),
      age_mean = round(rnorm(n(), 65, 5), 1),
      female_pct = pmin(0.8, pmax(0.2, rnorm(n(), 0.45, 0.10))),
      bmi_mean = round(rnorm(n(), 28, 2), 1),
      charlson = round(pmax(0, rnorm(n(), 1.5, 0.5)), 1),
      baseline_risk = plogis(rnorm(n(), qlogis(0.25), 0.4))
    )

  left_join(pw, info, by = "studlab")
}
