# Database Backend for Analysis Storage
# SQLite database for storing and retrieving NMA analyses
# Version 1.6.0

# =============================================================================
# Database Initialization
# =============================================================================

#' Initialize CNMA Analysis Database
#'
#' Creates a SQLite database for storing network meta-analysis projects,
#' datasets, results, and reports.
#'
#' @param db_path Path to database file (default: "cnma_analyses.db")
#'
#' @return Database connection
#' @export
#'
#' @examples
#' \dontrun{
#' db <- initialize_cnma_database("my_analyses.db")
#' }
initialize_cnma_database <- function(db_path = "cnma_analyses.db") {

  if (!requireNamespace("RSQLite", quietly = TRUE)) {
    stop("Package 'RSQLite' is required. Install with: install.packages('RSQLite')")
  }

  if (!requireNamespace("DBI", quietly = TRUE)) {
    stop("Package 'DBI' is required. Install with: install.packages('DBI')")
  }

  # Connect to database
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Create projects table
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS projects (
      project_id INTEGER PRIMARY KEY AUTOINCREMENT,
      project_name TEXT NOT NULL UNIQUE,
      description TEXT,
      created_date TEXT NOT NULL,
      last_modified TEXT NOT NULL,
      status TEXT DEFAULT 'active'
    )
  ")

  # Create datasets table
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS datasets (
      dataset_id INTEGER PRIMARY KEY AUTOINCREMENT,
      project_id INTEGER,
      dataset_name TEXT NOT NULL,
      n_studies INTEGER,
      n_treatments INTEGER,
      n_comparisons INTEGER,
      data_format TEXT,
      uploaded_date TEXT NOT NULL,
      FOREIGN KEY (project_id) REFERENCES projects(project_id)
    )
  ")

  # Create analyses table
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS analyses (
      analysis_id INTEGER PRIMARY KEY AUTOINCREMENT,
      dataset_id INTEGER,
      analysis_type TEXT NOT NULL,
      model_type TEXT,
      summary_measure TEXT,
      analysis_date TEXT NOT NULL,
      status TEXT DEFAULT 'completed',
      FOREIGN KEY (dataset_id) REFERENCES datasets(dataset_id)
    )
  ")

  # Create results table
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS results (
      result_id INTEGER PRIMARY KEY AUTOINCREMENT,
      analysis_id INTEGER,
      tau2 REAL,
      i2 REAL,
      q_statistic REAL,
      q_pvalue REAL,
      n_treatments INTEGER,
      top_treatment TEXT,
      top_pscore REAL,
      FOREIGN KEY (analysis_id) REFERENCES analyses(analysis_id)
    )
  ")

  # Create treatment rankings table
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS treatment_rankings (
      ranking_id INTEGER PRIMARY KEY AUTOINCREMENT,
      analysis_id INTEGER,
      treatment TEXT NOT NULL,
      p_score REAL,
      rank_position INTEGER,
      FOREIGN KEY (analysis_id) REFERENCES analyses(analysis_id)
    )
  ")

  # Create reports table
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS reports (
      report_id INTEGER PRIMARY KEY AUTOINCREMENT,
      analysis_id INTEGER,
      report_type TEXT NOT NULL,
      report_format TEXT,
      file_path TEXT,
      generated_date TEXT NOT NULL,
      FOREIGN KEY (analysis_id) REFERENCES analyses(analysis_id)
    )
  ")

  message("Database initialized: ", db_path)

  return(con)
}

# =============================================================================
# Project Management
# =============================================================================

#' Create New Project
#'
#' Creates a new project entry in the database.
#'
#' @param con Database connection
#' @param project_name Project name (must be unique)
#' @param description Project description
#'
#' @return Project ID
#' @export
#'
#' @examples
#' \dontrun{
#' db <- initialize_cnma_database()
#' project_id <- create_project(db, "Antidepressants NMA", "Comparing antidepressants")
#' }
create_project <- function(con, project_name, description = "") {

  # Check if project already exists
  existing <- DBI::dbGetQuery(
    con,
    "SELECT project_id FROM projects WHERE project_name = ?",
    params = list(project_name)
  )

  if (nrow(existing) > 0) {
    stop("Project '", project_name, "' already exists")
  }

  # Insert new project
  DBI::dbExecute(
    con,
    "INSERT INTO projects (project_name, description, created_date, last_modified)
     VALUES (?, ?, ?, ?)",
    params = list(
      project_name,
      description,
      as.character(Sys.time()),
      as.character(Sys.time())
    )
  )

  # Get project ID
  project_id <- DBI::dbGetQuery(con, "SELECT last_insert_rowid() as id")$id

  message("Project created with ID: ", project_id)

  return(project_id)
}

#' List All Projects
#'
#' Retrieves all projects from the database.
#'
#' @param con Database connection
#' @param status Filter by status (default: all)
#'
#' @return Data frame with project information
#' @export
#'
#' @examples
#' \dontrun{
#' db <- initialize_cnma_database()
#' projects <- list_projects(db)
#' }
list_projects <- function(con, status = c("all", "active", "archived")) {

  status <- match.arg(status)

  if (status == "all") {
    query <- "SELECT * FROM projects ORDER BY last_modified DESC"
    projects <- DBI::dbGetQuery(con, query)
  } else {
    query <- "SELECT * FROM projects WHERE status = ? ORDER BY last_modified DESC"
    projects <- DBI::dbGetQuery(con, query, params = list(status))
  }

  return(projects)
}

# =============================================================================
# Dataset Storage
# =============================================================================

#' Save Dataset to Database
#'
#' Stores a dataset in the database associated with a project.
#'
#' @param con Database connection
#' @param project_id Project ID
#' @param dataset_name Dataset name
#' @param data Data frame with study-level data
#'
#' @return Dataset ID
#' @export
#'
#' @examples
#' \dontrun{
#' db <- initialize_cnma_database()
#' project_id <- create_project(db, "My NMA")
#' data <- simulate_cnma_data(50)
#' dataset_id <- save_dataset(db, project_id, "Main Dataset", data)
#' }
save_dataset <- function(con, project_id, dataset_name, data) {

  # Check if project exists
  project_exists <- DBI::dbGetQuery(
    con,
    "SELECT COUNT(*) as n FROM projects WHERE project_id = ?",
    params = list(project_id)
  )$n > 0

  if (!project_exists) {
    stop("Project ID ", project_id, " does not exist")
  }

  # Calculate dataset statistics
  n_studies <- length(unique(data$studlab))
  n_treatments <- length(unique(c(data$treat1, data$treat2)))
  n_comparisons <- nrow(data)

  # Determine data format
  data_format <- ifelse(all(c("TE", "seTE") %in% names(data)), "contrast-level", "arm-level")

  # Insert dataset metadata
  DBI::dbExecute(
    con,
    "INSERT INTO datasets (project_id, dataset_name, n_studies, n_treatments, n_comparisons, data_format, uploaded_date)
     VALUES (?, ?, ?, ?, ?, ?, ?)",
    params = list(
      project_id,
      dataset_name,
      n_studies,
      n_treatments,
      n_comparisons,
      data_format,
      as.character(Sys.time())
    )
  )

  dataset_id <- DBI::dbGetQuery(con, "SELECT last_insert_rowid() as id")$id

  # Save actual data as serialized blob
  data_table_name <- paste0("dataset_", dataset_id)

  # Store data in separate table
  DBI::dbWriteTable(con, data_table_name, data, overwrite = TRUE)

  # Update project last_modified
  DBI::dbExecute(
    con,
    "UPDATE projects SET last_modified = ? WHERE project_id = ?",
    params = list(as.character(Sys.time()), project_id)
  )

  message("Dataset saved with ID: ", dataset_id)

  return(dataset_id)
}

#' Load Dataset from Database
#'
#' Retrieves a dataset from the database.
#'
#' @param con Database connection
#' @param dataset_id Dataset ID
#'
#' @return Data frame with study-level data
#' @export
#'
#' @examples
#' \dontrun{
#' db <- initialize_cnma_database()
#' data <- load_dataset(db, dataset_id = 1)
#' }
load_dataset <- function(con, dataset_id) {

  # Check if dataset exists
  dataset_exists <- DBI::dbGetQuery(
    con,
    "SELECT COUNT(*) as n FROM datasets WHERE dataset_id = ?",
    params = list(dataset_id)
  )$n > 0

  if (!dataset_exists) {
    stop("Dataset ID ", dataset_id, " does not exist")
  }

  # Load data
  data_table_name <- paste0("dataset_", dataset_id)

  if (!DBI::dbExistsTable(con, data_table_name)) {
    stop("Dataset table '", data_table_name, "' not found")
  }

  data <- DBI::dbReadTable(con, data_table_name)

  message("Dataset loaded: ", nrow(data), " rows")

  return(data)
}

# =============================================================================
# Analysis Storage
# =============================================================================

#' Save Analysis Results to Database
#'
#' Stores network meta-analysis results in the database.
#'
#' @param con Database connection
#' @param dataset_id Dataset ID
#' @param nma_results Network meta-analysis results
#' @param analysis_type Type of analysis
#'
#' @return Analysis ID
#' @export
#'
#' @examples
#' \dontrun{
#' db <- initialize_cnma_database()
#' data <- load_dataset(db, 1)
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' analysis_id <- save_analysis(db, dataset_id = 1, nma, "Frequentist NMA")
#' }
save_analysis <- function(con,
                         dataset_id,
                         nma_results,
                         analysis_type = "Frequentist NMA") {

  # Extract key results
  tau2 <- nma_results$tau^2
  i2 <- nma_results$I2
  q_stat <- nma_results$Q
  q_pval <- nma_results$pval.Q

  # Get treatment rankings
  rankings <- netmeta::netrank(nma_results)
  top_idx <- which.max(rankings$Pscore.random)
  top_treatment <- nma_results$trts[top_idx]
  top_pscore <- rankings$Pscore.random[top_idx]

  n_treatments <- length(nma_results$trts)

  # Insert analysis metadata
  DBI::dbExecute(
    con,
    "INSERT INTO analyses (dataset_id, analysis_type, model_type, summary_measure, analysis_date, status)
     VALUES (?, ?, ?, ?, ?, ?)",
    params = list(
      dataset_id,
      analysis_type,
      "Random effects",
      nma_results$sm,
      as.character(Sys.time()),
      "completed"
    )
  )

  analysis_id <- DBI::dbGetQuery(con, "SELECT last_insert_rowid() as id")$id

  # Insert results summary
  DBI::dbExecute(
    con,
    "INSERT INTO results (analysis_id, tau2, i2, q_statistic, q_pvalue, n_treatments, top_treatment, top_pscore)
     VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
    params = list(
      analysis_id,
      tau2,
      i2,
      q_stat,
      q_pval,
      n_treatments,
      top_treatment,
      top_pscore
    )
  )

  # Insert treatment rankings
  for (i in seq_along(nma_results$trts)) {
    DBI::dbExecute(
      con,
      "INSERT INTO treatment_rankings (analysis_id, treatment, p_score, rank_position)
       VALUES (?, ?, ?, ?)",
      params = list(
        analysis_id,
        nma_results$trts[i],
        rankings$Pscore.random[i],
        rank(-rankings$Pscore.random)[i]
      )
    )
  }

  # Save complete NMA object as serialized file
  nma_object_name <- paste0("nma_object_", analysis_id, ".rds")
  nma_dir <- file.path(dirname(DBI::dbGetInfo(con)$dbname), "nma_objects")

  if (!dir.exists(nma_dir)) {
    dir.create(nma_dir, recursive = TRUE)
  }

  nma_file_path <- file.path(nma_dir, nma_object_name)
  saveRDS(nma_results, nma_file_path)

  message("Analysis saved with ID: ", analysis_id)

  return(analysis_id)
}

#' Load Analysis from Database
#'
#' Retrieves a complete NMA analysis from the database.
#'
#' @param con Database connection
#' @param analysis_id Analysis ID
#'
#' @return NMA results object
#' @export
#'
#' @examples
#' \dontrun{
#' db <- initialize_cnma_database()
#' nma <- load_analysis(db, analysis_id = 1)
#' }
load_analysis <- function(con, analysis_id) {

  # Check if analysis exists
  analysis_exists <- DBI::dbGetQuery(
    con,
    "SELECT COUNT(*) as n FROM analyses WHERE analysis_id = ?",
    params = list(analysis_id)
  )$n > 0

  if (!analysis_exists) {
    stop("Analysis ID ", analysis_id, " does not exist")
  }

  # Load NMA object
  nma_object_name <- paste0("nma_object_", analysis_id, ".rds")
  nma_dir <- file.path(dirname(DBI::dbGetInfo(con)$dbname), "nma_objects")
  nma_file_path <- file.path(nma_dir, nma_object_name)

  if (!file.exists(nma_file_path)) {
    stop("NMA object file not found: ", nma_file_path)
  }

  nma_results <- readRDS(nma_file_path)

  message("Analysis loaded: ", analysis_id)

  return(nma_results)
}

# =============================================================================
# Query and Search Functions
# =============================================================================

#' Search Analyses by Criteria
#'
#' Searches for analyses matching specified criteria.
#'
#' @param con Database connection
#' @param project_id Filter by project ID (optional)
#' @param analysis_type Filter by analysis type (optional)
#' @param min_studies Minimum number of studies (optional)
#' @param date_from Filter by date from (optional)
#'
#' @return Data frame with matching analyses
#' @export
#'
#' @examples
#' \dontrun{
#' db <- initialize_cnma_database()
#' analyses <- search_analyses(db, min_studies = 10)
#' }
search_analyses <- function(con,
                           project_id = NULL,
                           analysis_type = NULL,
                           min_studies = NULL,
                           date_from = NULL) {

  query <- "
    SELECT
      a.analysis_id,
      p.project_name,
      d.dataset_name,
      d.n_studies,
      d.n_treatments,
      a.analysis_type,
      a.analysis_date,
      r.top_treatment,
      r.top_pscore
    FROM analyses a
    JOIN datasets d ON a.dataset_id = d.dataset_id
    JOIN projects p ON d.project_id = p.project_id
    LEFT JOIN results r ON a.analysis_id = r.analysis_id
    WHERE 1=1
  "

  params <- list()

  if (!is.null(project_id)) {
    query <- paste(query, "AND p.project_id = ?")
    params <- c(params, project_id)
  }

  if (!is.null(analysis_type)) {
    query <- paste(query, "AND a.analysis_type = ?")
    params <- c(params, analysis_type)
  }

  if (!is.null(min_studies)) {
    query <- paste(query, "AND d.n_studies >= ?")
    params <- c(params, min_studies)
  }

  if (!is.null(date_from)) {
    query <- paste(query, "AND a.analysis_date >= ?")
    params <- c(params, date_from)
  }

  query <- paste(query, "ORDER BY a.analysis_date DESC")

  results <- DBI::dbGetQuery(con, query, params = params)

  return(results)
}

#' Get Treatment Rankings Across Analyses
#'
#' Compares treatment rankings across multiple analyses.
#'
#' @param con Database connection
#' @param analysis_ids Vector of analysis IDs to compare
#'
#' @return Data frame with treatment rankings across analyses
#' @export
#'
#' @examples
#' \dontrun{
#' db <- initialize_cnma_database()
#' rankings <- get_rankings_comparison(db, c(1, 2, 3))
#' }
get_rankings_comparison <- function(con, analysis_ids) {

  rankings_list <- list()

  for (aid in analysis_ids) {
    query <- "
      SELECT treatment, p_score, rank_position
      FROM treatment_rankings
      WHERE analysis_id = ?
      ORDER BY rank_position
    "

    rankings <- DBI::dbGetQuery(con, query, params = list(aid))
    rankings$analysis_id <- aid
    rankings_list[[length(rankings_list) + 1]] <- rankings
  }

  all_rankings <- do.call(rbind, rankings_list)

  # Reshape for comparison
  comparison <- all_rankings %>%
    tidyr::pivot_wider(
      names_from = analysis_id,
      values_from = c(p_score, rank_position),
      names_prefix = "Analysis_"
    )

  return(comparison)
}

# =============================================================================
# Database Utilities
# =============================================================================

#' Close Database Connection
#'
#' Properly closes the database connection.
#'
#' @param con Database connection
#'
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' db <- initialize_cnma_database()
#' # ... do work ...
#' close_database(db)
#' }
close_database <- function(con) {
  DBI::dbDisconnect(con)
  message("Database connection closed")
  invisible(NULL)
}

#' Get Database Summary
#'
#' Retrieves summary statistics about the database.
#'
#' @param con Database connection
#'
#' @return List with database statistics
#' @export
#'
#' @examples
#' \dontrun{
#' db <- initialize_cnma_database()
#' summary <- get_database_summary(db)
#' print(summary)
#' }
get_database_summary <- function(con) {

  n_projects <- DBI::dbGetQuery(con, "SELECT COUNT(*) as n FROM projects")$n
  n_datasets <- DBI::dbGetQuery(con, "SELECT COUNT(*) as n FROM datasets")$n
  n_analyses <- DBI::dbGetQuery(con, "SELECT COUNT(*) as n FROM analyses")$n
  n_reports <- DBI::dbGetQuery(con, "SELECT COUNT(*) as n FROM reports")$n

  # Get database size
  db_path <- DBI::dbGetInfo(con)$dbname
  db_size_mb <- file.size(db_path) / 1024 / 1024

  summary <- list(
    n_projects = n_projects,
    n_datasets = n_datasets,
    n_analyses = n_analyses,
    n_reports = n_reports,
    database_size_mb = round(db_size_mb, 2),
    database_path = db_path
  )

  class(summary) <- "cnma_db_summary"
  return(summary)
}

#' @export
print.cnma_db_summary <- function(x, ...) {
  cat("CNMA Database Summary\n")
  cat("=====================\n\n")
  cat("Database path:", x$database_path, "\n")
  cat("Database size:", x$database_size_mb, "MB\n\n")
  cat("Contents:\n")
  cat("  Projects:", x$n_projects, "\n")
  cat("  Datasets:", x$n_datasets, "\n")
  cat("  Analyses:", x$n_analyses, "\n")
  cat("  Reports:", x$n_reports, "\n")
  invisible(x)
}
