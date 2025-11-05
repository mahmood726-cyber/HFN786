# Machine Learning Predictions for Network Meta-Analysis
# Advanced pattern detection, treatment ranking prediction, and outcome forecasting
# Version 1.6.0

# =============================================================================
# Machine Learning Treatment Ranking Prediction
# =============================================================================

#' Predict Treatment Rankings Using Machine Learning
#'
#' Uses random forest or gradient boosting to predict treatment rankings based
#' on study characteristics and network structure.
#'
#' @param data Data frame with study-level data
#' @param nma_results Network meta-analysis results
#' @param predictors Character vector of predictor variables
#' @param method Machine learning method: "rf" (random forest) or "gbm" (gradient boosting)
#' @param train_proportion Proportion of data for training (default: 0.8)
#' @param cv_folds Number of cross-validation folds (default: 5)
#' @param seed Random seed for reproducibility
#'
#' @return List with predictions, model performance, variable importance
#' @export
#'
#' @examples
#' \dontrun{
#' data <- simulate_cnma_data(50)
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' ml_pred <- ml_predict_rankings(data, nma, predictors = c("year", "age_mean"))
#' print(ml_pred)
#' plot(ml_pred)
#' }
ml_predict_rankings <- function(data,
                               nma_results,
                               predictors = NULL,
                               method = c("rf", "gbm"),
                               train_proportion = 0.8,
                               cv_folds = 5,
                               seed = 42) {

  method <- match.arg(method)
  set.seed(seed)

  # Extract rankings from NMA
  rankings <- netmeta::netrank(nma_results)
  treatment_scores <- data.frame(
    Treatment = nma_results$trts,
    P_score = rankings$Pscore.random,
    Rank = rank(-rankings$Pscore.random)
  )

  # If no predictors specified, use available numeric columns
  if (is.null(predictors)) {
    numeric_cols <- names(data)[sapply(data, is.numeric)]
    exclude_cols <- c("TE", "seTE", "n")
    predictors <- setdiff(numeric_cols, exclude_cols)
  }

  # Aggregate study data by treatment
  treatment_data <- data %>%
    dplyr::group_by(treat1) %>%
    dplyr::summarise(
      across(all_of(predictors), ~mean(.x, na.rm = TRUE)),
      n_studies = dplyr::n()
    ) %>%
    dplyr::rename(Treatment = treat1)

  # Merge with rankings
  ml_data <- dplyr::left_join(treatment_scores, treatment_data, by = "Treatment")

  # Remove rows with missing data
  ml_data <- ml_data[complete.cases(ml_data), ]

  if (nrow(ml_data) < 5) {
    stop("Insufficient data for machine learning (need at least 5 treatments with complete data)")
  }

  # Split into train/test
  n_train <- floor(nrow(ml_data) * train_proportion)
  train_idx <- sample(1:nrow(ml_data), n_train)
  train_data <- ml_data[train_idx, ]
  test_data <- ml_data[-train_idx, ]

  # Prepare formula
  formula_str <- paste("P_score ~", paste(c(predictors, "n_studies"), collapse = " + "))
  formula_obj <- as.formula(formula_str)

  # Train model
  if (method == "rf") {
    # Random Forest
    if (!requireNamespace("randomForest", quietly = TRUE)) {
      stop("Package 'randomForest' is required. Install with: install.packages('randomForest')")
    }

    model <- randomForest::randomForest(
      formula_obj,
      data = train_data,
      ntree = 500,
      importance = TRUE
    )

    var_importance <- randomForest::importance(model)

  } else {
    # Gradient Boosting
    if (!requireNamespace("gbm", quietly = TRUE)) {
      stop("Package 'gbm' is required. Install with: install.packages('gbm')")
    }

    model <- gbm::gbm(
      formula_obj,
      data = train_data,
      distribution = "gaussian",
      n.trees = 500,
      interaction.depth = 3,
      shrinkage = 0.01,
      cv.folds = cv_folds
    )

    var_importance <- summary(model, plotit = FALSE)
  }

  # Predictions
  if (nrow(test_data) > 0) {
    if (method == "rf") {
      predictions <- predict(model, newdata = test_data)
    } else {
      predictions <- predict(model, newdata = test_data, n.trees = 500)
    }

    # Calculate performance metrics
    rmse <- sqrt(mean((test_data$P_score - predictions)^2))
    mae <- mean(abs(test_data$P_score - predictions))
    r_squared <- cor(test_data$P_score, predictions)^2

    performance <- list(
      RMSE = rmse,
      MAE = mae,
      R_squared = r_squared,
      predictions = data.frame(
        Treatment = test_data$Treatment,
        Actual = test_data$P_score,
        Predicted = predictions,
        Error = test_data$P_score - predictions
      )
    )
  } else {
    performance <- list(
      RMSE = NA,
      MAE = NA,
      R_squared = NA,
      predictions = NULL,
      note = "No test data available (insufficient sample size)"
    )
  }

  # Feature importance
  if (method == "rf") {
    importance_df <- data.frame(
      Variable = rownames(var_importance),
      Importance = var_importance[, "%IncMSE"],
      stringsAsFactors = FALSE
    )
  } else {
    importance_df <- var_importance
    names(importance_df) <- c("Variable", "Importance")
  }

  importance_df <- importance_df[order(-importance_df$Importance), ]

  result <- list(
    model = model,
    method = method,
    performance = performance,
    importance = importance_df,
    predictors = c(predictors, "n_studies"),
    train_data = train_data,
    test_data = test_data,
    all_data = ml_data
  )

  class(result) <- "cnma_ml_predictions"
  return(result)
}

#' @export
print.cnma_ml_predictions <- function(x, ...) {
  cat("CNMA Machine Learning Treatment Ranking Predictions\n")
  cat("====================================================\n\n")

  cat("Method:", toupper(x$method), "\n")
  cat("Predictors:", paste(x$predictors, collapse = ", "), "\n")
  cat("Training samples:", nrow(x$train_data), "\n")
  cat("Test samples:", nrow(x$test_data), "\n\n")

  if (!is.na(x$performance$RMSE)) {
    cat("Model Performance (on test set):\n")
    cat("  RMSE:", round(x$performance$RMSE, 4), "\n")
    cat("  MAE:", round(x$performance$MAE, 4), "\n")
    cat("  R²:", round(x$performance$R_squared, 4), "\n\n")
  }

  cat("Variable Importance (Top 5):\n")
  print(head(x$importance, 5), row.names = FALSE)

  invisible(x)
}

#' @export
plot.cnma_ml_predictions <- function(x, type = c("importance", "predictions", "residuals"), ...) {
  type <- match.arg(type)

  if (type == "importance") {
    # Variable importance plot
    top_vars <- head(x$importance, min(10, nrow(x$importance)))

    barplot(
      top_vars$Importance,
      names.arg = top_vars$Variable,
      las = 2,
      col = "#3498db",
      main = "Variable Importance for Treatment Ranking Prediction",
      ylab = "Importance",
      cex.names = 0.8
    )

  } else if (type == "predictions" && !is.null(x$performance$predictions)) {
    # Actual vs Predicted
    pred_df <- x$performance$predictions

    plot(
      pred_df$Actual,
      pred_df$Predicted,
      pch = 19,
      col = "#3498db",
      xlab = "Actual P-score",
      ylab = "Predicted P-score",
      main = "Actual vs Predicted Treatment Rankings"
    )
    abline(0, 1, col = "red", lty = 2, lwd = 2)
    text(pred_df$Actual, pred_df$Predicted, labels = pred_df$Treatment, pos = 3, cex = 0.7)

    # Add R² to plot
    r2_text <- paste("R² =", round(x$performance$R_squared, 3))
    legend("topleft", legend = r2_text, bty = "n")

  } else if (type == "residuals" && !is.null(x$performance$predictions)) {
    # Residual plot
    pred_df <- x$performance$predictions

    plot(
      pred_df$Predicted,
      pred_df$Error,
      pch = 19,
      col = "#e74c3c",
      xlab = "Predicted P-score",
      ylab = "Residual (Actual - Predicted)",
      main = "Residual Plot"
    )
    abline(h = 0, col = "red", lty = 2, lwd = 2)
    text(pred_df$Predicted, pred_df$Error, labels = pred_df$Treatment, pos = 3, cex = 0.7)
  } else {
    cat("No predictions available for plotting.\n")
  }
}

# =============================================================================
# Pattern Detection in Network Meta-Analysis
# =============================================================================

#' Detect Patterns in Network Meta-Analysis Data
#'
#' Uses clustering and anomaly detection to identify patterns in study characteristics,
#' treatment effects, and network structure.
#'
#' @param data Data frame with study-level data
#' @param nma_results Network meta-analysis results
#' @param method Clustering method: "kmeans", "hierarchical", or "dbscan"
#' @param n_clusters Number of clusters (for kmeans and hierarchical)
#' @param variables Variables to use for clustering
#'
#' @return List with cluster assignments, centroids, and visualizations
#' @export
#'
#' @examples
#' \dontrun{
#' data <- simulate_cnma_data(50)
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' patterns <- ml_detect_patterns(data, nma, method = "kmeans", n_clusters = 3)
#' print(patterns)
#' plot(patterns)
#' }
ml_detect_patterns <- function(data,
                               nma_results,
                               method = c("kmeans", "hierarchical", "dbscan"),
                               n_clusters = 3,
                               variables = NULL) {

  method <- match.arg(method)

  # If no variables specified, use available numeric columns
  if (is.null(variables)) {
    numeric_cols <- names(data)[sapply(data, is.numeric)]
    exclude_cols <- c("seTE")
    variables <- setdiff(numeric_cols, exclude_cols)
  }

  # Prepare clustering data
  cluster_data <- data[, variables, drop = FALSE]
  cluster_data <- cluster_data[complete.cases(cluster_data), ]

  # Scale data
  scaled_data <- scale(cluster_data)

  # Perform clustering
  if (method == "kmeans") {
    cluster_result <- stats::kmeans(scaled_data, centers = n_clusters, nstart = 25)
    clusters <- cluster_result$cluster
    centroids <- cluster_result$centers

  } else if (method == "hierarchical") {
    dist_matrix <- dist(scaled_data)
    hc <- stats::hclust(dist_matrix, method = "ward.D2")
    clusters <- stats::cutree(hc, k = n_clusters)

    # Calculate centroids
    centroids <- aggregate(scaled_data, by = list(clusters), FUN = mean)
    rownames(centroids) <- centroids[, 1]
    centroids <- as.matrix(centroids[, -1])

  } else if (method == "dbscan") {
    if (!requireNamespace("dbscan", quietly = TRUE)) {
      stop("Package 'dbscan' is required. Install with: install.packages('dbscan')")
    }

    db_result <- dbscan::dbscan(scaled_data, eps = 0.5, minPts = 5)
    clusters <- db_result$cluster

    # Calculate centroids for each cluster
    unique_clusters <- unique(clusters[clusters != 0])
    n_clusters <- length(unique_clusters)

    centroids <- matrix(NA, nrow = n_clusters, ncol = ncol(scaled_data))
    for (i in seq_along(unique_clusters)) {
      cluster_id <- unique_clusters[i]
      centroids[i, ] <- colMeans(scaled_data[clusters == cluster_id, , drop = FALSE])
    }
    rownames(centroids) <- paste0("Cluster_", unique_clusters)
    colnames(centroids) <- colnames(scaled_data)
  }

  # Cluster statistics
  cluster_summary <- data.frame(
    Cluster = 1:n_clusters,
    Size = as.numeric(table(clusters)),
    stringsAsFactors = FALSE
  )

  # Add original data with cluster assignments
  result_data <- cbind(data[complete.cases(data[, variables]), ], Cluster = clusters)

  # Silhouette analysis (for kmeans and hierarchical)
  silhouette_score <- NA
  if (method %in% c("kmeans", "hierarchical") && n_clusters > 1) {
    if (requireNamespace("cluster", quietly = TRUE)) {
      dist_matrix <- dist(scaled_data)
      sil <- cluster::silhouette(clusters, dist_matrix)
      silhouette_score <- mean(sil[, 3])
    }
  }

  result <- list(
    method = method,
    clusters = clusters,
    centroids = centroids,
    cluster_summary = cluster_summary,
    data = result_data,
    scaled_data = scaled_data,
    variables = variables,
    silhouette_score = silhouette_score
  )

  class(result) <- "cnma_ml_patterns"
  return(result)
}

#' @export
print.cnma_ml_patterns <- function(x, ...) {
  cat("CNMA Pattern Detection Results\n")
  cat("==============================\n\n")

  cat("Method:", x$method, "\n")
  cat("Variables:", paste(x$variables, collapse = ", "), "\n")
  cat("Number of clusters:", nrow(x$cluster_summary), "\n\n")

  cat("Cluster Summary:\n")
  print(x$cluster_summary, row.names = FALSE)

  if (!is.na(x$silhouette_score)) {
    cat("\nSilhouette Score:", round(x$silhouette_score, 3), "\n")
    cat("(0.71-1.0: Strong structure, 0.51-0.70: Reasonable, 0.26-0.50: Weak, <0.25: No structure)\n")
  }

  invisible(x)
}

#' @export
plot.cnma_ml_patterns <- function(x, type = c("scatter", "heatmap", "dendrogram"), ...) {
  type <- match.arg(type)

  if (type == "scatter") {
    # PCA for 2D visualization
    if (ncol(x$scaled_data) >= 2) {
      pca <- stats::prcomp(x$scaled_data)

      # Determine number of colors needed
      n_clusters <- length(unique(x$clusters))
      colors <- grDevices::rainbow(n_clusters)

      plot(
        pca$x[, 1],
        pca$x[, 2],
        col = colors[x$clusters],
        pch = 19,
        xlab = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
        ylab = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)"),
        main = "Pattern Detection: PCA Visualization"
      )
      legend("topright", legend = paste("Cluster", 1:n_clusters), col = colors, pch = 19)
    }

  } else if (type == "heatmap") {
    # Heatmap of cluster centroids
    if (requireNamespace("pheatmap", quietly = TRUE)) {
      pheatmap::pheatmap(
        x$centroids,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        main = "Cluster Centroids Heatmap"
      )
    } else {
      heatmap(
        x$centroids,
        main = "Cluster Centroids Heatmap",
        col = grDevices::colorRampPalette(c("blue", "white", "red"))(100)
      )
    }

  } else if (type == "dendrogram" && x$method == "hierarchical") {
    # Dendrogram for hierarchical clustering
    dist_matrix <- dist(x$scaled_data)
    hc <- stats::hclust(dist_matrix, method = "ward.D2")
    plot(hc, main = "Hierarchical Clustering Dendrogram", xlab = "", sub = "")
  }
}

# =============================================================================
# Anomaly Detection
# =============================================================================

#' Detect Anomalies in Network Meta-Analysis Data
#'
#' Identifies outlier studies based on treatment effects, standard errors,
#' and study characteristics using isolation forest or local outlier factor.
#'
#' @param data Data frame with study-level data
#' @param nma_results Network meta-analysis results
#' @param method Anomaly detection method: "isolation_forest" or "lof"
#' @param contamination Expected proportion of outliers (default: 0.1)
#'
#' @return List with anomaly scores and classifications
#' @export
#'
#' @examples
#' \dontrun{
#' data <- simulate_cnma_data(50)
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' anomalies <- ml_detect_anomalies(data, nma)
#' print(anomalies)
#' }
ml_detect_anomalies <- function(data,
                                nma_results,
                                method = c("isolation_forest", "lof"),
                                contamination = 0.1) {

  method <- match.arg(method)

  # Prepare data for anomaly detection
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  anomaly_data <- data[, numeric_cols, drop = FALSE]
  anomaly_data <- anomaly_data[complete.cases(anomaly_data), ]

  # Scale data
  scaled_data <- scale(anomaly_data)

  # Detect anomalies
  if (method == "isolation_forest") {
    if (!requireNamespace("solitude", quietly = TRUE)) {
      stop("Package 'solitude' is required. Install with: install.packages('solitude')")
    }

    iso_forest <- solitude::isolationForest$new()
    iso_forest$fit(as.data.frame(scaled_data))
    anomaly_scores <- iso_forest$predict(as.data.frame(scaled_data))

    # Classify as anomaly if score exceeds threshold
    threshold <- stats::quantile(anomaly_scores$anomaly_score, 1 - contamination)
    is_anomaly <- anomaly_scores$anomaly_score > threshold

  } else if (method == "lof") {
    if (!requireNamespace("dbscan", quietly = TRUE)) {
      stop("Package 'dbscan' is required. Install with: install.packages('dbscan')")
    }

    lof_scores <- dbscan::lof(scaled_data, k = 5)

    # Classify as anomaly if LOF score exceeds threshold
    threshold <- stats::quantile(lof_scores, 1 - contamination)
    is_anomaly <- lof_scores > threshold

    anomaly_scores <- data.frame(
      anomaly_score = lof_scores,
      stringsAsFactors = FALSE
    )
  }

  # Create result dataframe
  result_data <- cbind(
    data[complete.cases(data[, numeric_cols]), ],
    anomaly_score = anomaly_scores$anomaly_score,
    is_anomaly = is_anomaly
  )

  result <- list(
    method = method,
    data = result_data,
    anomalies = result_data[is_anomaly, ],
    n_anomalies = sum(is_anomaly),
    contamination = contamination,
    threshold = threshold
  )

  class(result) <- "cnma_ml_anomalies"
  return(result)
}

#' @export
print.cnma_ml_anomalies <- function(x, ...) {
  cat("CNMA Anomaly Detection Results\n")
  cat("===============================\n\n")

  cat("Method:", x$method, "\n")
  cat("Contamination level:", x$contamination, "\n")
  cat("Threshold:", round(x$threshold, 4), "\n")
  cat("Number of anomalies detected:", x$n_anomalies, "\n\n")

  if (x$n_anomalies > 0) {
    cat("Anomalous Studies:\n")
    print(x$anomalies[, c("studlab", "treat1", "treat2", "TE", "anomaly_score")], row.names = FALSE)
  } else {
    cat("No anomalies detected.\n")
  }

  invisible(x)
}

# =============================================================================
# Treatment Effect Prediction for New Studies
# =============================================================================

#' Predict Treatment Effects for New Studies
#'
#' Uses trained models to predict expected treatment effects for new studies
#' based on their characteristics.
#'
#' @param ml_model Trained ML model from ml_predict_rankings()
#' @param new_data Data frame with characteristics of new studies
#'
#' @return Data frame with predicted treatment effects
#' @export
#'
#' @examples
#' \dontrun{
#' # Train model
#' ml_model <- ml_predict_rankings(data, nma, predictors = c("year", "age_mean"))
#'
#' # Predict for new studies
#' new_studies <- data.frame(year = 2025, age_mean = 65, n_studies = 1)
#' predictions <- ml_predict_new_studies(ml_model, new_studies)
#' }
ml_predict_new_studies <- function(ml_model, new_data) {

  if (!inherits(ml_model, "cnma_ml_predictions")) {
    stop("ml_model must be output from ml_predict_rankings()")
  }

  # Check that new_data has all required predictors
  missing_predictors <- setdiff(ml_model$predictors, names(new_data))
  if (length(missing_predictors) > 0) {
    stop("New data is missing required predictors: ", paste(missing_predictors, collapse = ", "))
  }

  # Make predictions
  if (ml_model$method == "rf") {
    predictions <- predict(ml_model$model, newdata = new_data)
  } else {
    predictions <- predict(ml_model$model, newdata = new_data, n.trees = 500)
  }

  result <- data.frame(
    new_data,
    Predicted_P_score = predictions,
    Predicted_Rank = rank(-predictions)
  )

  return(result)
}
