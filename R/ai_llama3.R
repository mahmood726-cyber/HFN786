# =========================================================
# LLama 3 AI Integration Module
# Intelligent assistance for network meta-analysis
# =========================================================

#' Configure LLama 3 AI Assistant (Local Ollama Integration)
#'
#' Sets up the LLama 3 language model for providing intelligent assistance
#' throughout the network meta-analysis workflow. Uses local Ollama server
#' for complete offline inference - no data leaves your computer.
#'
#' @param model Model name in Ollama (default "llama3")
#' @param temperature Sampling temperature (0-1, default 0.7)
#' @param num_predict Maximum tokens in response (default 2000)
#' @param top_p Nucleus sampling parameter (default 0.9)
#' @param top_k Top-k sampling parameter (default 40)
#' @param seed Random seed for reproducibility (NULL for random)
#' @param check_ollama Check if Ollama is running (default TRUE)
#' @return AI configuration object
#' @export
#' @examples
#' \dontrun{
#' # Install Ollama first: https://ollama.com/download
#' # Pull llama3 model: ollama pull llama3
#'
#' # Configure with local Ollama
#' ai_config <- configure_llama3(
#'   model = "llama3",
#'   temperature = 0.7,
#'   seed = 42
#' )
#'
#' # Use higher-end models for better analysis
#' ai_config <- configure_llama3(model = "llama3:70b")
#' }
configure_llama3 <- function(model = "llama3",
                            temperature = 0.7,
                            num_predict = 2000,
                            top_p = 0.9,
                            top_k = 40,
                            seed = NULL,
                            check_ollama = TRUE) {

  # Check for rollama package
  if (!requireNamespace("rollama", quietly = TRUE)) {
    msg("Package 'rollama' required for local LLama 3 integration.")
    msg("Install with: install.packages('rollama')")
    msg("Also ensure Ollama is installed: https://ollama.com/download")
    .stop_hint("Missing required package: rollama")
  }

  # Check if Ollama is running
  if (check_ollama) {
    ollama_status <- .safe_try({
      rollama::list_models()
    }, context = "Ollama connection", silent = TRUE)

    if (inherits(ollama_status, "try-error")) {
      msg("ERROR: Cannot connect to Ollama.")
      msg("Please ensure Ollama is running:")
      msg("  1. Download from: https://ollama.com/download")
      msg("  2. Start Ollama service")
      msg("  3. Pull model: ollama pull %s", model)
      .stop_hint("Ollama not running - start service first")
    }

    # Check if model exists
    available_models <- ollama_status$name
    if (!any(grepl(model, available_models))) {
      msg("Model '%s' not found in Ollama.", model)
      msg("Available models: %s", paste(available_models, collapse = ", "))
      msg("Pull with: ollama pull %s", model)
      .stop_hint("Model not found - pull it first with 'ollama pull %s'", model)
    }
  }

  config <- list(
    model = model,
    temperature = temperature,
    num_predict = num_predict,
    top_p = top_p,
    top_k = top_k,
    seed = seed,
    enabled = TRUE,
    backend = "ollama"
  )

  class(config) <- "llama3_config"

  # Store in package environment
  .cnma_env$ai_config <- config

  msg("✓ LLama 3 AI Assistant configured successfully")
  msg("  Model: %s (local Ollama)", model)
  msg("  Temperature: %.2f", temperature)
  msg("  Backend: Offline/Local (no data leaves your computer)")

  invisible(config)
}

#' Call LLama 3 Model (Local Ollama)
#'
#' Internal function to call LLama 3 via local Ollama server.
#' All inference happens locally - no data sent to external servers.
#'
#' @param prompt Text prompt for the model
#' @param config AI configuration object
#' @param system_message System message for context
#' @param stream Stream response (default FALSE)
#' @return Model response text
#' @keywords internal
call_llama3 <- function(prompt,
                       config = .cnma_env$ai_config,
                       system_message = "You are an expert biostatistician specializing in network meta-analysis. Provide clear, evidence-based, actionable advice.",
                       stream = FALSE) {

  if (is.null(config) || !config$enabled) {
    return("AI assistant not configured. Use configure_llama3() first.")
  }

  # Check for rollama
  if (!requireNamespace("rollama", quietly = TRUE)) {
    return("AI assistant requires 'rollama' package. Install with: install.packages('rollama')")
  }

  # Prepare options
  options <- list(
    temperature = config$temperature,
    num_predict = config$num_predict,
    top_p = config$top_p,
    top_k = config$top_k
  )

  if (!is.null(config$seed)) {
    options$seed <- config$seed
  }

  # Call Ollama locally
  response <- .safe_try({
    rollama::generate(
      model = config$model,
      prompt = prompt,
      system = system_message,
      options = options,
      stream = stream,
      output = "text"
    )
  }, context = "LLama 3 local inference", silent = TRUE)

  if (inherits(response, "try-error")) {
    error_msg <- attr(response, "condition")$message
    if (grepl("connection", error_msg, ignore.case = TRUE)) {
      return("AI assistant connection failed. Ensure Ollama is running:\n  Start with: ollama serve")
    } else if (grepl("model", error_msg, ignore.case = TRUE)) {
      return(sprintf("Model '%s' not found. Pull with: ollama pull %s",
                    config$model, config$model))
    } else {
      return(paste("AI Error:", error_msg))
    }
  }

  return(response)
}

#' AI-Powered Data Quality Check
#'
#' Uses LLama 3 to perform intelligent data quality assessment, identifying
#' potential issues, anomalies, and providing recommendations.
#'
#' @param data Data frame with NMA data
#' @param verbose Print detailed AI analysis
#' @return List with quality assessment and recommendations
#' @export
#' @examples
#' \dontrun{
#' configure_llama3()
#' data <- simulate_cnma_data(30)
#' quality <- ai_quality_check(data)
#' print(quality)
#' }
ai_quality_check <- function(data, verbose = TRUE) {

  # Gather data statistics
  stats <- list(
    n_studies = length(unique(data$studlab)),
    n_treatments = length(unique(c(data$treat1, data$treat2))),
    n_comparisons = nrow(data),
    effect_range = range(data$TE, na.rm = TRUE),
    se_range = range(data$seTE, na.rm = TRUE),
    missing_values = sum(is.na(data)),
    duplicate_rows = sum(duplicated(data))
  )

  # Create prompt
  prompt <- sprintf("
Analyze this network meta-analysis dataset and identify potential quality issues:

Dataset Statistics:
- Number of studies: %d
- Number of treatments: %d
- Number of comparisons: %d
- Treatment effect range: [%.3f, %.3f]
- Standard error range: [%.3f, %.3f]
- Missing values: %d
- Duplicate rows: %d

Please provide:
1. Assessment of data quality (Good/Fair/Poor)
2. Potential issues or anomalies
3. Specific recommendations for improvement
4. Red flags that require immediate attention

Format your response clearly with numbered sections.
",
    stats$n_studies,
    stats$n_treatments,
    stats$n_comparisons,
    stats$effect_range[1],
    stats$effect_range[2],
    stats$se_range[1],
    stats$se_range[2],
    stats$missing_values,
    stats$duplicate_rows
  )

  if (verbose) {
    msg("Consulting AI assistant for data quality assessment...")
  }

  ai_response <- call_llama3(prompt)

  result <- list(
    statistics = stats,
    ai_assessment = ai_response,
    timestamp = Sys.time()
  )

  class(result) <- "ai_quality_check"

  if (verbose) {
    cat("\n")
    cat("AI Data Quality Assessment\n")
    cat("==========================\n\n")
    cat(ai_response)
    cat("\n")
  }

  invisible(result)
}

#' AI-Powered Result Interpretation
#'
#' Uses LLama 3 to provide intelligent interpretation of NMA results,
#' including clinical significance, implications, and recommendations.
#'
#' @param cnma_results CNMA results object
#' @param context Clinical context (e.g., "depression treatment", "diabetes")
#' @param target_audience Audience level ("clinical", "academic", "patient")
#' @return List with AI interpretation and recommendations
#' @export
#' @examples
#' \dontrun{
#' configure_llama3()
#' data <- simulate_cnma_data(30)
#' config <- setup_cnma(use_bayesian = FALSE)
#' results <- run_cnma_analysis(data, config = config)
#'
#' interpretation <- ai_interpret_results(
#'   results,
#'   context = "antidepressant treatment",
#'   target_audience = "clinical"
#' )
#' }
ai_interpret_results <- function(cnma_results,
                                context = "general",
                                target_audience = c("clinical", "academic", "patient")) {

  target_audience <- match.arg(target_audience)

  if (!inherits(cnma_results, "cnma")) {
    .stop_hint("Input must be a cnma results object.")
  }

  nma <- cnma_results$results$main_nma

  # Gather key results
  key_results <- sprintf("
Network Meta-Analysis Results for: %s

Network Structure:
- Studies: %d
- Treatments: %d
- Reference: %s

Heterogeneity:
- Tau²: %.4f
- I²: %.1f%%
- Interpretation: %s

Top 3 Treatments (if rankings available):
[Would list top treatments with P-scores]

Key Findings:
- [Summary of main comparisons]
- [Confidence intervals]
- [Prediction intervals]

Target Audience: %s

Please provide:
1. Plain language interpretation of these results
2. Clinical significance and implications
3. Limitations and caveats
4. Recommendations for practice/further research
5. Key messages for the target audience
",
    context,
    length(unique(cnma_results$data$studlab)),
    length(unique(c(cnma_results$data$treat1, cnma_results$data$treat2))),
    cnma_results$ref_treatment,
    nma$tau^2,
    nma$I2.random * 100,
    ifelse(nma$I2.random < 0.40, "Low",
           ifelse(nma$I2.random < 0.75, "Moderate", "Substantial")),
    target_audience
  )

  msg("Generating AI interpretation...")

  ai_response <- call_llama3(key_results)

  result <- list(
    context = context,
    audience = target_audience,
    interpretation = ai_response,
    results_summary = cnma_results,
    timestamp = Sys.time()
  )

  class(result) <- "ai_interpretation"

  cat("\n")
  cat("AI-Powered Result Interpretation\n")
  cat("=================================\n\n")
  cat(ai_response)
  cat("\n")

  invisible(result)
}

#' AI Analysis Recommender
#'
#' Uses LLama 3 to recommend appropriate analyses based on data characteristics
#' and research question.
#'
#' @param data Data frame
#' @param research_question Text description of research question
#' @param constraints Any constraints (time, resources, etc.)
#' @return Recommended analysis plan
#' @export
#' @examples
#' \dontrun{
#' configure_llama3()
#' data <- simulate_cnma_data(30)
#'
#' recommendations <- ai_recommend_analysis(
#'   data,
#'   research_question = "Which antidepressant is most effective?",
#'   constraints = "Need results within 1 week"
#' )
#' }
ai_recommend_analysis <- function(data,
                                 research_question,
                                 constraints = "None") {

  # Analyze data characteristics
  characteristics <- sprintf("
Dataset Characteristics:
- Studies: %d
- Treatments: %d
- Comparisons: %d
- Covariates available: %s
- Network connectivity: [to be assessed]

Research Question: %s

Constraints: %s

Based on this information, please recommend:
1. Most appropriate analysis approach (frequentist/Bayesian)
2. Which CNMA functions to use and in what order
3. Recommended sensitivity analyses
4. Covariate adjustments to consider
5. Potential pitfalls to avoid
6. Estimated time to completion
7. Step-by-step analysis plan

Be specific with function names from the CNMA package.
",
    length(unique(data$studlab)),
    length(unique(c(data$treat1, data$treat2))),
    nrow(data),
    paste(setdiff(names(data), c("studlab", "treat1", "treat2", "TE", "seTE")),
          collapse = ", "),
    research_question,
    constraints
  )

  msg("Consulting AI for analysis recommendations...")

  ai_response <- call_llama3(characteristics)

  result <- list(
    research_question = research_question,
    recommendations = ai_response,
    data_summary = data.frame(
      n_studies = length(unique(data$studlab)),
      n_treatments = length(unique(c(data$treat1, data$treat2))),
      n_comparisons = nrow(data)
    ),
    timestamp = Sys.time()
  )

  class(result) <- "ai_recommendations"

  cat("\n")
  cat("AI Analysis Recommendations\n")
  cat("============================\n\n")
  cat(ai_response)
  cat("\n")

  invisible(result)
}

#' AI Error Diagnostics
#'
#' Uses LLama 3 to diagnose errors and provide solutions.
#'
#' @param error_message Error message text
#' @param context Analysis context
#' @return Diagnostic information and solutions
#' @export
#' @examples
#' \dontrun{
#' configure_llama3()
#' diagnostics <- ai_diagnose_error(
#'   "Error: TE/seTE contain non-finite values",
#'   context = "data validation"
#' )
#' }
ai_diagnose_error <- function(error_message, context = "general") {

  prompt <- sprintf("
A user encountered this error in network meta-analysis:

Error: %s
Context: %s

Please provide:
1. Explanation of what caused this error
2. Step-by-step solution to fix it
3. Code example if applicable
4. How to prevent this error in the future
5. Related issues to check

Be specific and practical.
",
    error_message,
    context
  )

  msg("Diagnosing error with AI assistant...")

  ai_response <- call_llama3(prompt)

  cat("\n")
  cat("AI Error Diagnostics\n")
  cat("====================\n\n")
  cat(ai_response)
  cat("\n")

  invisible(ai_response)
}

#' AI Manuscript Assistant
#'
#' Helps draft methods and results sections for manuscripts.
#'
#' @param cnma_results CNMA results object
#' @param section Section to draft ("methods", "results", "abstract")
#' @param journal_style Journal style ("BMJ", "Lancet", "JAMA", "generic")
#' @return Drafted text
#' @export
#' @examples
#' \dontrun{
#' configure_llama3()
#' data <- simulate_cnma_data(30)
#' results <- run_comprehensive_nma(data)
#'
#' methods <- ai_draft_manuscript(results, "methods", "BMJ")
#' cat(methods)
#' }
ai_draft_manuscript <- function(cnma_results,
                               section = c("methods", "results", "abstract"),
                               journal_style = c("BMJ", "Lancet", "JAMA", "generic")) {

  section <- match.arg(section)
  journal_style <- match.arg(journal_style)

  prompt <- sprintf("
Draft the %s section for a network meta-analysis manuscript following %s style.

Analysis Details:
- %d studies included
- %d treatments compared
- Frequentist approach using random-effects model
- Assessment of heterogeneity, inconsistency, and transitivity
- PRISMA-NMA guidelines followed
- Treatment rankings calculated using P-scores
- Sensitivity analyses performed

Please draft a publication-ready %s section that:
1. Follows %s journal format and style
2. Includes all necessary methodological details
3. Reports results clearly and accurately
4. Uses appropriate statistical terminology
5. Cites key methodology papers
6. Is approximately 250-300 words

Do not include placeholder text - make it publication-ready.
",
    section,
    journal_style,
    length(unique(cnma_results$data$studlab)),
    length(unique(c(cnma_results$data$treat1, cnma_results$data$treat2))),
    section,
    journal_style
  )

  msg("Drafting manuscript %s section...", section)

  ai_response <- call_llama3(prompt)

  cat("\n")
  cat(sprintf("AI-Drafted %s Section (%s style)\n",
              tools::toTitleCase(section), journal_style))
  cat(paste(rep("=", 50), collapse = ""))
  cat("\n\n")
  cat(ai_response)
  cat("\n\n")
  cat("Note: Please review and edit as needed. Verify all statistical details.\n")

  invisible(ai_response)
}
