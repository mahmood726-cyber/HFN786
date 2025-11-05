# Example Network Meta-Analysis Datasets
# Based on published literature

#' Antidepressants Network Meta-Analysis Data
#'
#' Network meta-analysis data for antidepressant efficacy from published
#' systematic review. This dataset includes 117 pairwise comparisons from
#' 28 studies comparing 12 different antidepressants for major depression.
#'
#' @format A data frame with 117 rows and 12 variables:
#' \describe{
#'   \item{studlab}{Study identifier}
#'   \item{treat1}{First treatment in comparison}
#'   \item{treat2}{Second treatment in comparison}
#'   \item{TE}{Treatment effect (log odds ratio)}
#'   \item{seTE}{Standard error of treatment effect}
#'   \item{year}{Publication year}
#'   \item{n_total}{Total sample size}
#'   \item{rob}{Risk of bias (low, unclear, high)}
#'   \item{setting}{Clinical setting (inpatient, outpatient, mixed)}
#'   \item{age_mean}{Mean age of participants}
#'   \item{female_pct}{Percentage of female participants}
#'   \item{severity}{Baseline depression severity (mild, moderate, severe)}
#' }
#'
#' @source Simulated based on structure from Cipriani et al. (2018) Lancet.
#' @examples
#' data(antidepressants)
#' head(antidepressants)
#'
#' # Run analysis
#' \donttest{
#' config <- setup_cnma(sm = "OR", use_bayesian = FALSE)
#' results <- run_cnma_analysis(antidepressants, config = config)
#' print(results)
#' }
"antidepressants"

#' Diabetes Medications Network Meta-Analysis Data
#'
#' Network meta-analysis comparing glucose-lowering medications for type 2
#' diabetes. Includes 178 pairwise comparisons from 56 studies comparing
#' 9 drug classes based on HbA1c reduction.
#'
#' @format A data frame with 178 rows and 13 variables:
#' \describe{
#'   \item{studlab}{Study identifier}
#'   \item{treat1}{First treatment in comparison}
#'   \item{treat2}{Second treatment in comparison}
#'   \item{TE}{Treatment effect (mean difference in HbA1c)}
#'   \item{seTE}{Standard error of treatment effect}
#'   \item{year}{Publication year}
#'   \item{duration_weeks}{Treatment duration in weeks}
#'   \item{n_total}{Total sample size}
#'   \item{baseline_hba1c}{Baseline HbA1c percentage}
#'   \item{age_mean}{Mean age of participants}
#'   \item{diabetes_duration}{Mean diabetes duration in years}
#'   \item{bmi_mean}{Mean body mass index}
#'   \item{trial_quality}{Trial quality score (0-10)}
#' }
#'
#' @source Simulated based on structure from published diabetes NMAs.
#' @examples
#' data(diabetes_meds)
#' head(diabetes_meds)
#'
#' # Assess transitivity
#' \donttest{
#' transitivity <- assess_transitivity(
#'   diabetes_meds,
#'   variables = c("baseline_hba1c", "age_mean", "bmi_mean")
#' )
#' print(transitivity)
#' }
"diabetes_meds"

#' Statins for Cardiovascular Prevention
#'
#' Network meta-analysis of statins for primary prevention of cardiovascular
#' events. Includes 89 pairwise comparisons from 23 randomized trials.
#'
#' @format A data frame with 89 rows and 11 variables:
#' \describe{
#'   \item{studlab}{Study identifier}
#'   \item{treat1}{First treatment in comparison}
#'   \item{treat2}{Second treatment in comparison}
#'   \item{TE}{Treatment effect (log hazard ratio)}
#'   \item{seTE}{Standard error of treatment effect}
#'   \item{year}{Publication year}
#'   \item{followup_years}{Mean follow-up duration in years}
#'   \item{age_mean}{Mean age of participants}
#'   \item{male_pct}{Percentage of male participants}
#'   \item{baseline_ldl}{Baseline LDL cholesterol (mg/dL)}
#'   \item{cvd_history}{Percentage with prior CVD history}
#' }
#'
#' @source Simulated based on structure from statin meta-analyses.
#' @examples
#' data(statins)
#' head(statins)
#'
#' # Complete analysis
#' \donttest{
#' results <- run_comprehensive_nma(
#'   statins,
#'   ref_treatment = "Placebo",
#'   generate_plots = TRUE
#' )
#' }
"statins"

# Generate the actual datasets
set.seed(42)

# Antidepressants dataset
antidepressants <- local({
  treatments <- c("Placebo", "Fluoxetine", "Sertraline", "Paroxetine",
                 "Citalopram", "Escitalopram", "Venlafaxine", "Duloxetine",
                 "Bupropion", "Mirtazapine", "Amitriptyline", "Nortriptyline")

  n_studies <- 28
  study_ids <- sprintf("Study_%02d", 1:n_studies)

  # True effects (log OR)
  true_effects <- c(
    Placebo = 0,
    Fluoxetine = -0.30,
    Sertraline = -0.32,
    Paroxetine = -0.28,
    Citalopram = -0.25,
    Escitalopram = -0.35,
    Venlafaxine = -0.40,
    Duloxetine = -0.38,
    Bupropion = -0.22,
    Mirtazapine = -0.33,
    Amitriptyline = -0.35,
    Nortriptyline = -0.30
  )

  # Create network
  comparisons <- list()
  for(i in 1:n_studies) {
    # Random treatment selection
    n_arms <- sample(2:3, 1, prob = c(0.7, 0.3))
    study_treats <- sample(treatments, n_arms)

    # Pairwise comparisons
    if(length(study_treats) >= 2) {
      pairs <- combn(study_treats, 2, simplify = FALSE)

      for(pair in pairs) {
        true_diff <- true_effects[pair[1]] - true_effects[pair[2]]
        se <- runif(1, 0.15, 0.35)
        te <- rnorm(1, true_diff, se)

        comparisons[[length(comparisons) + 1]] <- data.frame(
          studlab = study_ids[i],
          treat1 = pair[1],
          treat2 = pair[2],
          TE = te,
          seTE = se,
          year = sample(2010:2023, 1),
          n_total = round(rnorm(1, 150, 50)),
          rob = sample(c("low", "unclear", "high"), 1, prob = c(0.5, 0.3, 0.2)),
          setting = sample(c("outpatient", "inpatient", "mixed"), 1, prob = c(0.6, 0.2, 0.2)),
          age_mean = round(rnorm(1, 42, 8), 1),
          female_pct = round(runif(1, 0.45, 0.75), 2),
          severity = sample(c("mild", "moderate", "severe"), 1, prob = c(0.2, 0.6, 0.2)),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  do.call(rbind, comparisons)
})

# Diabetes medications dataset
diabetes_meds <- local({
  treatments <- c("Placebo", "Metformin", "Sulfonylureas", "DPP4i",
                 "GLP1RA", "SGLT2i", "Insulin", "Thiazolidinediones", "Meglitinides")

  n_studies <- 56
  study_ids <- sprintf("DiabStudy_%02d", 1:n_studies)

  true_effects <- c(
    Placebo = 0,
    Metformin = -0.80,
    Sulfonylureas = -0.75,
    DPP4i = -0.70,
    GLP1RA = -1.00,
    SGLT2i = -0.85,
    Insulin = -1.20,
    Thiazolidinediones = -0.65,
    Meglitinides = -0.60
  )

  comparisons <- list()
  for(i in 1:n_studies) {
    n_arms <- sample(2:3, 1, prob = c(0.8, 0.2))
    study_treats <- sample(treatments, n_arms)

    if(length(study_treats) >= 2) {
      pairs <- combn(study_treats, 2, simplify = FALSE)

      for(pair in pairs) {
        true_diff <- true_effects[pair[1]] - true_effects[pair[2]]
        se <- runif(1, 0.08, 0.20)
        te <- rnorm(1, true_diff, se)

        comparisons[[length(comparisons) + 1]] <- data.frame(
          studlab = study_ids[i],
          treat1 = pair[1],
          treat2 = pair[2],
          TE = te,
          seTE = se,
          year = sample(2008:2024, 1),
          duration_weeks = sample(c(12, 24, 52), 1),
          n_total = round(rnorm(1, 280, 100)),
          baseline_hba1c = round(rnorm(1, 8.2, 0.8), 1),
          age_mean = round(rnorm(1, 58, 7), 1),
          diabetes_duration = round(rnorm(1, 6.5, 3), 1),
          bmi_mean = round(rnorm(1, 31, 4), 1),
          trial_quality = round(runif(1, 6, 10)),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  do.call(rbind, comparisons)
})

# Statins dataset
statins <- local({
  treatments <- c("Placebo", "Atorvastatin", "Simvastatin", "Rosuvastatin",
                 "Pravastatin", "Lovastatin", "Fluvastatin")

  n_studies <- 23
  study_ids <- sprintf("StatinRCT_%02d", 1:n_studies)

  true_effects <- c(
    Placebo = 0,
    Atorvastatin = log(0.73),
    Simvastatin = log(0.75),
    Rosuvastatin = log(0.72),
    Pravastatin = log(0.78),
    Lovastatin = log(0.76),
    Fluvastatin = log(0.80)
  )

  comparisons <- list()
  for(i in 1:n_studies) {
    n_arms <- sample(2:3, 1, prob = c(0.85, 0.15))
    study_treats <- sample(treatments, n_arms)

    if(length(study_treats) >= 2) {
      pairs <- combn(study_treats, 2, simplify = FALSE)

      for(pair in pairs) {
        true_diff <- true_effects[pair[1]] - true_effects[pair[2]]
        se <- runif(1, 0.10, 0.25)
        te <- rnorm(1, true_diff, se)

        comparisons[[length(comparisons) + 1]] <- data.frame(
          studlab = study_ids[i],
          treat1 = pair[1],
          treat2 = pair[2],
          TE = te,
          seTE = se,
          year = sample(1994:2022, 1),
          followup_years = round(runif(1, 2, 6), 1),
          age_mean = round(rnorm(1, 62, 6), 1),
          male_pct = round(runif(1, 0.50, 0.80), 2),
          baseline_ldl = round(rnorm(1, 135, 20)),
          cvd_history = round(runif(1, 0, 0.30), 2),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  do.call(rbind, comparisons)
})

# Save datasets
usethis::use_data(antidepressants, diabetes_meds, statins, overwrite = TRUE)
