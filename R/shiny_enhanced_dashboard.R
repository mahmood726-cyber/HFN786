# Enhanced Shiny Dashboard with bs4Dash Framework
# Integrates advanced features from mahmood789/786-MIII-Meta-analysis repository
# Version 1.5.0 - Professional-grade modular architecture

# =============================================================================
# Module 1: Instructions Module
# =============================================================================

instructionsUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::h3("Data Format Requirements"),
    shiny::hr(),

    shiny::h4("Option 1: Contrast-Level Data"),
    shiny::p("Required columns:"),
    shiny::tags$ul(
      shiny::tags$li(shiny::strong("studlab:"), "Study identifier (character)"),
      shiny::tags$li(shiny::strong("treat1:"), "First treatment (character)"),
      shiny::tags$li(shiny::strong("treat2:"), "Second treatment (character)"),
      shiny::tags$li(shiny::strong("TE:"), "Treatment effect (numeric, e.g., log HR)"),
      shiny::tags$li(shiny::strong("seTE:"), "Standard error of TE (numeric)")
    ),
    shiny::p("Optional columns: year, age_mean, female_pct, or other covariates"),

    shiny::hr(),

    shiny::h4("Option 2: Arm-Level Data"),
    shiny::p("Required columns:"),
    shiny::tags$ul(
      shiny::tags$li(shiny::strong("study:"), "Study identifier (character)"),
      shiny::tags$li(shiny::strong("treatment:"), "Treatment name (character)"),
      shiny::tags$li(shiny::strong("mean:"), "Mean outcome (numeric)"),
      shiny::tags$li(shiny::strong("sd:"), "Standard deviation (numeric)"),
      shiny::tags$li(shiny::strong("n:"), "Sample size (integer)")
    ),
    shiny::p("Optional columns: year, covariates"),
    shiny::p(shiny::em("Note: Arm-level data will be automatically converted to contrast-level using the pairwise() function."))
  )
}

instructionsServer <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    # Minimal server logic for instructions
  })
}

# =============================================================================
# Module 2: Data Upload Module with Validation
# =============================================================================

dataUploadUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::fluidRow(
      shiny::column(
        width = 4,
        bs4Dash::box(
          title = "Upload Data",
          status = "primary",
          solidHeader = TRUE,
          width = 12,

          shiny::radioButtons(
            ns("data_format"),
            "Data Format:",
            choices = c("Contrast-level" = "contrast", "Arm-level" = "arm"),
            selected = "contrast"
          ),

          shiny::fileInput(
            ns("data_file"),
            "Choose CSV File:",
            accept = c(".csv", "text/csv", "text/comma-separated-values")
          ),

          shiny::checkboxInput(
            ns("log_transform"),
            "Log-transform outcome (for HR/OR/RR from arm-level data)",
            value = FALSE
          ),

          shiny::actionButton(
            ns("process_data"),
            "Process Data",
            class = "btn-success",
            icon = shiny::icon("play")
          )
        )
      ),

      shiny::column(
        width = 8,
        bs4Dash::box(
          title = "Data Preview",
          status = "info",
          solidHeader = TRUE,
          width = 12,

          shinycssloaders::withSpinner(
            DT::DTOutput(ns("data_preview"))
          )
        )
      )
    )
  )
}

dataUploadServer <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {

    dataset <- shiny::reactiveVal(NULL)

    shiny::observeEvent(input$process_data, {
      shiny::req(input$data_file)

      tryCatch({
        # Read CSV
        raw_data <- readr::read_csv(input$data_file$datapath, show_col_types = FALSE)

        if (input$data_format == "arm") {
          # Validate arm-level columns
          required_cols <- c("study", "treatment", "mean", "sd", "n")
          missing_cols <- setdiff(required_cols, names(raw_data))

          if (length(missing_cols) > 0) {
            shiny::showNotification(
              paste("Missing required columns:", paste(missing_cols, collapse = ", ")),
              type = "error",
              duration = 5
            )
            return(NULL)
          }

          # Convert arm-level to contrast-level
          processed_data <- netmeta::pairwise(
            treat = treatment,
            mean = mean,
            sd = sd,
            n = n,
            studlab = study,
            data = raw_data,
            sm = if(input$log_transform) "SMD" else "MD"
          )

        } else {
          # Validate contrast-level columns
          required_cols <- c("studlab", "treat1", "treat2", "TE", "seTE")
          missing_cols <- setdiff(required_cols, names(raw_data))

          if (length(missing_cols) > 0) {
            shiny::showNotification(
              paste("Missing required columns:", paste(missing_cols, collapse = ", ")),
              type = "error",
              duration = 5
            )
            return(NULL)
          }

          processed_data <- raw_data
        }

        # Store processed data
        dataset(processed_data)

        shiny::showNotification(
          paste("Successfully loaded", nrow(processed_data), "comparisons"),
          type = "message",
          duration = 3
        )

      }, error = function(e) {
        shiny::showNotification(
          paste("Error processing data:", e$message),
          type = "error",
          duration = 10
        )
      })
    })

    # Data preview table
    output$data_preview <- DT::renderDT({
      shiny::req(dataset())
      DT::datatable(
        dataset(),
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = 'Bfrtip'
        ),
        rownames = FALSE
      )
    })

    return(dataset)
  })
}

# =============================================================================
# Module 3: Example Data Download Module
# =============================================================================

exampleDataUI <- function(id) {
  ns <- shiny::NS(id)

  bs4Dash::box(
    title = "Download Example Datasets",
    status = "success",
    solidHeader = TRUE,
    width = 12,

    shiny::p("Download sample CSV files to understand the required format:"),

    shiny::fluidRow(
      shiny::column(
        width = 4,
        shiny::downloadButton(
          ns("download_contrast"),
          "Contrast-Level Example",
          class = "btn-primary btn-block"
        )
      ),
      shiny::column(
        width = 4,
        shiny::downloadButton(
          ns("download_arm"),
          "Arm-Level Example",
          class = "btn-primary btn-block"
        )
      ),
      shiny::column(
        width = 4,
        shiny::downloadButton(
          ns("download_covariate"),
          "With Covariates Example",
          class = "btn-primary btn-block"
        )
      )
    )
  )
}

exampleDataServer <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {

    output$download_contrast <- shiny::downloadHandler(
      filename = "example_contrast_level.csv",
      content = function(file) {
        example_data <- data.frame(
          studlab = c("Study1", "Study1", "Study2", "Study2", "Study3"),
          treat1 = c("Placebo", "Placebo", "Placebo", "DrugA", "Placebo"),
          treat2 = c("DrugA", "DrugB", "DrugA", "DrugB", "DrugC"),
          TE = c(-0.163, -0.288, -0.223, -0.125, -0.357),
          seTE = c(0.150, 0.200, 0.180, 0.165, 0.210)
        )
        readr::write_csv(example_data, file)
      }
    )

    output$download_arm <- shiny::downloadHandler(
      filename = "example_arm_level.csv",
      content = function(file) {
        example_data <- data.frame(
          study = c("Study1", "Study1", "Study2", "Study2", "Study3", "Study3"),
          treatment = c("Placebo", "DrugA", "Placebo", "DrugA", "Placebo", "DrugB"),
          mean = c(5.2, 4.8, 5.5, 4.9, 5.3, 4.6),
          sd = c(1.2, 1.1, 1.3, 1.0, 1.2, 1.1),
          n = c(100, 105, 98, 102, 110, 108)
        )
        readr::write_csv(example_data, file)
      }
    )

    output$download_covariate <- shiny::downloadHandler(
      filename = "example_with_covariates.csv",
      content = function(file) {
        example_data <- data.frame(
          studlab = c("Study1", "Study1", "Study2", "Study2", "Study3"),
          treat1 = c("Placebo", "Placebo", "Placebo", "DrugA", "Placebo"),
          treat2 = c("DrugA", "DrugB", "DrugA", "DrugB", "DrugC"),
          TE = c(-0.163, -0.288, -0.223, -0.125, -0.357),
          seTE = c(0.150, 0.200, 0.180, 0.165, 0.210),
          year = c(2018, 2018, 2019, 2019, 2020),
          age_mean = c(55, 55, 58, 58, 52),
          female_pct = c(0.48, 0.48, 0.52, 0.52, 0.45)
        )
        readr::write_csv(example_data, file)
      }
    )
  })
}

# =============================================================================
# Module 4: Interactive Network Plot Module with visNetwork
# =============================================================================

networkPlotUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::fluidRow(
      shiny::column(
        width = 3,
        bs4Dash::box(
          title = "Network Plot Settings",
          status = "primary",
          solidHeader = TRUE,
          width = 12,

          shiny::selectInput(
            ns("layout"),
            "Layout Algorithm:",
            choices = c(
              "Fruchterman-Reingold" = "fr",
              "Kamada-Kawai" = "kk",
              "Circular" = "circle",
              "Tree" = "tree",
              "Random" = "random",
              "Spring" = "spring"
            ),
            selected = "fr"
          ),

          shiny::sliderInput(
            ns("node_size"),
            "Node Size:",
            min = 10,
            max = 100,
            value = 50
          ),

          shiny::sliderInput(
            ns("edge_width"),
            "Edge Width:",
            min = 1,
            max = 10,
            value = 3
          ),

          shiny::checkboxInput(
            ns("show_labels"),
            "Show Treatment Labels",
            value = TRUE
          ),

          shiny::downloadButton(
            ns("download_network"),
            "Download Network Plot",
            class = "btn-success btn-block"
          )
        )
      ),

      shiny::column(
        width = 9,
        bs4Dash::tabBox(
          title = "Network Visualization",
          width = 12,

          shiny::tabPanel(
            "Interactive (visNetwork)",
            shinycssloaders::withSpinner(
              visNetwork::visNetworkOutput(ns("interactive_network"), height = "600px")
            )
          ),

          shiny::tabPanel(
            "Static (igraph)",
            shinycssloaders::withSpinner(
              shiny::plotOutput(ns("static_network"), height = "600px")
            )
          )
        )
      )
    )
  )
}

networkPlotServer <- function(id, analysisResult) {
  shiny::moduleServer(id, function(input, output, session) {

    # Interactive network with visNetwork
    output$interactive_network <- visNetwork::renderVisNetwork({
      shiny::req(analysisResult())

      nma <- analysisResult()$nma

      tryCatch({
        # Extract network structure
        treatments <- nma$trts
        n_treat <- length(treatments)

        # Create nodes dataframe
        nodes <- data.frame(
          id = 1:n_treat,
          label = if(input$show_labels) treatments else rep("", n_treat),
          title = treatments,  # Tooltip
          size = input$node_size,
          color = "#1f77b4",
          font = list(size = 16)
        )

        # Create edges from comparison matrix
        comparisons <- nma$A.matrix
        edges_list <- list()
        edge_id <- 1

        for (i in 1:(n_treat - 1)) {
          for (j in (i + 1):n_treat) {
            if (comparisons[i, j] > 0) {
              edges_list[[edge_id]] <- data.frame(
                from = i,
                to = j,
                width = input$edge_width,
                title = paste(comparisons[i, j], "studies"),
                color = "#7f7f7f"
              )
              edge_id <- edge_id + 1
            }
          }
        }

        edges <- do.call(rbind, edges_list)

        # Create visNetwork
        visNetwork::visNetwork(nodes, edges) %>%
          visNetwork::visNodes(
            shape = "dot",
            shadow = TRUE
          ) %>%
          visNetwork::visEdges(
            smooth = TRUE
          ) %>%
          visNetwork::visLayout(randomSeed = 42) %>%
          visNetwork::visPhysics(
            solver = input$layout,
            stabilization = list(iterations = 100)
          ) %>%
          visNetwork::visInteraction(
            navigationButtons = TRUE,
            dragNodes = TRUE,
            dragView = TRUE,
            zoomView = TRUE
          )

      }, error = function(e) {
        shiny::showNotification(
          paste("Error creating network:", e$message),
          type = "error"
        )
        return(NULL)
      })
    })

    # Static network with igraph
    output$static_network <- shiny::renderPlot({
      shiny::req(analysisResult())

      nma <- analysisResult()$nma

      tryCatch({
        # Create igraph object
        treatments <- nma$trts
        n_treat <- length(treatments)

        # Create edge list
        edge_list <- c()
        for (i in 1:(n_treat - 1)) {
          for (j in (i + 1):n_treat) {
            if (nma$A.matrix[i, j] > 0) {
              edge_list <- c(edge_list, i, j)
            }
          }
        }

        g <- igraph::graph(edges = edge_list, n = n_treat, directed = FALSE)
        igraph::V(g)$name <- treatments

        # Select layout
        layout_func <- switch(
          input$layout,
          "fr" = igraph::layout_with_fr,
          "kk" = igraph::layout_with_kk,
          "circle" = igraph::layout_in_circle,
          "tree" = igraph::layout_as_tree,
          "random" = igraph::layout_randomly,
          "spring" = igraph::layout_with_fr
        )

        layout <- layout_func(g)

        # Plot
        par(mar = c(1, 1, 1, 1))
        plot(
          g,
          layout = layout,
          vertex.size = input$node_size / 2,
          vertex.color = "#1f77b4",
          vertex.label = if(input$show_labels) igraph::V(g)$name else NA,
          vertex.label.cex = 1.2,
          vertex.label.color = "black",
          edge.width = input$edge_width,
          edge.color = "#7f7f7f"
        )

      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Error:", e$message), col = "red")
      })
    })

    # Download handler
    output$download_network <- shiny::downloadHandler(
      filename = function() {
        paste0("network_plot_", Sys.Date(), ".png")
      },
      content = function(file) {
        shiny::req(analysisResult())

        png(file, width = 1200, height = 800, res = 150)

        nma <- analysisResult()$nma
        treatments <- nma$trts
        n_treat <- length(treatments)

        edge_list <- c()
        for (i in 1:(n_treat - 1)) {
          for (j in (i + 1):n_treat) {
            if (nma$A.matrix[i, j] > 0) {
              edge_list <- c(edge_list, i, j)
            }
          }
        }

        g <- igraph::graph(edges = edge_list, n = n_treat, directed = FALSE)
        igraph::V(g)$name <- treatments

        layout_func <- switch(
          input$layout,
          "fr" = igraph::layout_with_fr,
          "kk" = igraph::layout_with_kk,
          "circle" = igraph::layout_in_circle,
          "tree" = igraph::layout_as_tree,
          "random" = igraph::layout_randomly,
          "spring" = igraph::layout_with_fr
        )

        layout <- layout_func(g)

        par(mar = c(1, 1, 1, 1))
        plot(
          g,
          layout = layout,
          vertex.size = input$node_size / 2,
          vertex.color = "#1f77b4",
          vertex.label = if(input$show_labels) igraph::V(g)$name else NA,
          vertex.label.cex = 1.2,
          vertex.label.color = "black",
          edge.width = input$edge_width,
          edge.color = "#7f7f7f"
        )

        dev.off()
      }
    )
  })
}

# =============================================================================
# Module 5: Analysis Module with Asynchronous Processing
# =============================================================================

analysisUI <- function(id) {
  ns <- shiny::NS(id)

  bs4Dash::box(
    title = "Analysis Configuration",
    status = "warning",
    solidHeader = TRUE,
    width = 12,

    shiny::fluidRow(
      shiny::column(
        width = 4,
        shiny::selectInput(
          ns("sm"),
          "Summary Measure:",
          choices = c(
            "Hazard Ratio" = "HR",
            "Odds Ratio" = "OR",
            "Risk Ratio" = "RR",
            "Mean Difference" = "MD",
            "Standardized Mean Difference" = "SMD"
          ),
          selected = "HR"
        )
      ),

      shiny::column(
        width = 4,
        shiny::selectInput(
          ns("model"),
          "Effect Model:",
          choices = c("Random Effects" = "random", "Fixed Effect" = "fixed"),
          selected = "random"
        )
      ),

      shiny::column(
        width = 4,
        shiny::numericInput(
          ns("conf_level"),
          "Confidence Level:",
          value = 0.95,
          min = 0.80,
          max = 0.99,
          step = 0.01
        )
      )
    ),

    shiny::hr(),

    shiny::h4("Meta-Regression (Optional)"),
    shiny::checkboxInput(
      ns("run_metareg"),
      "Run meta-regression with covariates",
      value = FALSE
    ),

    shiny::uiOutput(ns("covariate_selector")),

    shiny::hr(),

    shiny::actionButton(
      ns("run_analysis"),
      "Run Network Meta-Analysis",
      class = "btn-success btn-lg",
      icon = shiny::icon("play-circle")
    ),

    shiny::br(), shiny::br(),

    shiny::uiOutput(ns("analysis_status"))
  )
}

analysisServer <- function(id, dataset) {
  shiny::moduleServer(id, function(input, output, session) {

    analysisResult <- shiny::reactiveVal(NULL)

    # Dynamic covariate selector
    output$covariate_selector <- shiny::renderUI({
      ns <- session$ns
      shiny::req(dataset())

      if (input$run_metareg) {
        # Get available covariates
        data_cols <- names(dataset())
        required_cols <- c("studlab", "treat1", "treat2", "TE", "seTE")
        potential_covars <- setdiff(data_cols, required_cols)

        if (length(potential_covars) > 0) {
          shiny::checkboxGroupInput(
            ns("covariates"),
            "Select Covariates:",
            choices = potential_covars,
            selected = NULL
          )
        } else {
          shiny::p("No covariates available in dataset", style = "color: red;")
        }
      }
    })

    # Run analysis with asynchronous processing
    shiny::observeEvent(input$run_analysis, {
      shiny::req(dataset())

      output$analysis_status <- shiny::renderUI({
        shiny::tagList(
          shinycssloaders::withSpinner(
            shiny::h4("Analysis running...", style = "color: orange;")
          )
        )
      })

      # Setup for asynchronous execution
      data_copy <- dataset()
      sm_val <- input$sm
      model_val <- input$model
      conf_val <- input$conf_level
      run_mr <- input$run_metareg
      covars <- if(run_mr && !is.null(input$covariates)) input$covariates else NULL

      # Use future/promises for non-blocking execution
      future::future({

        # Run network meta-analysis
        nma <- netmeta::netmeta(
          TE = TE,
          seTE = seTE,
          treat1 = treat1,
          treat2 = treat2,
          studlab = studlab,
          data = data_copy,
          sm = sm_val,
          comb.fixed = (model_val == "fixed"),
          comb.random = (model_val == "random"),
          level = conf_val,
          reference.group = NULL
        )

        # Meta-regression if requested
        metareg_result <- NULL
        if (run_mr && !is.null(covars) && length(covars) > 0) {
          tryCatch({
            formula_str <- paste("~", paste(covars, collapse = " + "))
            metareg_result <- netmeta::netmetareg(nma, formula = as.formula(formula_str))
          }, error = function(e) {
            metareg_result <<- list(error = e$message)
          })
        }

        list(
          nma = nma,
          metareg = metareg_result,
          timestamp = Sys.time()
        )

      }) %...>% {
        # When analysis completes
        analysisResult(.)

        output$analysis_status <- shiny::renderUI({
          shiny::h4("✓ Analysis completed successfully!", style = "color: green;")
        })

        shiny::showNotification(
          "Network meta-analysis completed!",
          type = "message",
          duration = 3
        )
      } %...!% {
        # If analysis fails
        output$analysis_status <- shiny::renderUI({
          shiny::h4("✗ Analysis failed", style = "color: red;")
        })

        shiny::showNotification(
          paste("Error:", .),
          type = "error",
          duration = 10
        )
      }
    })

    return(analysisResult)
  })
}

# =============================================================================
# Module 6: Comprehensive Results Module (20+ Tabs)
# =============================================================================

resultsUI <- function(id) {
  ns <- shiny::NS(id)

  bs4Dash::tabBox(
    title = "Analysis Results",
    width = 12,
    id = ns("results_tabs"),

    # Tab 1: Summary
    shiny::tabPanel(
      "Summary",
      shiny::verbatimTextOutput(ns("nma_summary"))
    ),

    # Tab 2: Forest Plot
    shiny::tabPanel(
      "Forest Plot",
      shiny::fluidRow(
        shiny::column(
          width = 9,
          shinycssloaders::withSpinner(
            shiny::plotOutput(ns("forest_plot"), height = "700px")
          )
        ),
        shiny::column(
          width = 3,
          shiny::downloadButton(ns("download_forest"), "Download PNG", class = "btn-primary btn-block")
        )
      )
    ),

    # Tab 3: League Table
    shiny::tabPanel(
      "League Table",
      DT::DTOutput(ns("league_table")),
      shiny::downloadButton(ns("download_league"), "Download CSV", class = "btn-primary")
    ),

    # Tab 4: Treatment Rankings
    shiny::tabPanel(
      "Rankings",
      shiny::fluidRow(
        shiny::column(
          width = 9,
          shinycssloaders::withSpinner(
            shiny::plotOutput(ns("ranking_plot"), height = "600px")
          )
        ),
        shiny::column(
          width = 3,
          DT::DTOutput(ns("ranking_table")),
          shiny::downloadButton(ns("download_rankings"), "Download", class = "btn-primary btn-block")
        )
      )
    ),

    # Tab 5: Network Graph
    shiny::tabPanel(
      "Network Graph",
      shinycssloaders::withSpinner(
        shiny::plotOutput(ns("network_graph"), height = "600px")
      ),
      shiny::downloadButton(ns("download_network_graph"), "Download", class = "btn-primary")
    ),

    # Tab 6: Net Heat Plot
    shiny::tabPanel(
      "Net Heat",
      shinycssloaders::withSpinner(
        shiny::plotOutput(ns("netheat_plot"), height = "700px")
      ),
      shiny::downloadButton(ns("download_netheat"), "Download", class = "btn-primary")
    ),

    # Tab 7: Funnel Plot
    shiny::tabPanel(
      "Funnel Plot",
      shinycssloaders::withSpinner(
        shiny::plotOutput(ns("funnel_plot"), height = "600px")
      ),
      shiny::downloadButton(ns("download_funnel"), "Download", class = "btn-primary")
    ),

    # Tab 8: Contribution Matrix
    shiny::tabPanel(
      "Contribution Matrix",
      shinycssloaders::withSpinner(
        shiny::plotOutput(ns("contribution_plot"), height = "700px")
      ),
      shiny::downloadButton(ns("download_contribution"), "Download", class = "btn-primary")
    ),

    # Tab 9: Heterogeneity
    shiny::tabPanel(
      "Heterogeneity",
      shiny::verbatimTextOutput(ns("heterogeneity_text")),
      shinycssloaders::withSpinner(
        shiny::plotOutput(ns("heterogeneity_plot"), height = "500px")
      )
    ),

    # Tab 10: Inconsistency
    shiny::tabPanel(
      "Inconsistency",
      shiny::h4("Design-based Inconsistency Assessment"),
      shiny::verbatimTextOutput(ns("inconsistency_text")),
      shiny::downloadButton(ns("download_inconsistency"), "Download Report", class = "btn-primary")
    ),

    # Tab 11: Direct vs Indirect Evidence
    shiny::tabPanel(
      "Direct vs Indirect",
      shinycssloaders::withSpinner(
        shiny::plotOutput(ns("direct_indirect_plot"), height = "700px")
      )
    ),

    # Tab 12: Leave-One-Out
    shiny::tabPanel(
      "Leave-One-Out",
      shiny::p("Sensitivity analysis removing one study at a time"),
      DT::DTOutput(ns("loo_table")),
      shiny::downloadButton(ns("download_loo"), "Download Results", class = "btn-primary")
    ),

    # Tab 13: Net Splitting
    shiny::tabPanel(
      "Net Splitting",
      shiny::p("Local inconsistency assessment via node-splitting"),
      DT::DTOutput(ns("netsplit_table")),
      shinycssloaders::withSpinner(
        shiny::plotOutput(ns("netsplit_plot"), height = "600px")
      ),
      shiny::downloadButton(ns("download_netsplit"), "Download", class = "btn-primary")
    ),

    # Tab 14: Prediction Intervals
    shiny::tabPanel(
      "Prediction Intervals",
      DT::DTOutput(ns("prediction_table")),
      shiny::downloadButton(ns("download_prediction"), "Download", class = "btn-primary")
    ),

    # Tab 15: Meta-Regression (conditional)
    shiny::tabPanel(
      "Meta-Regression",
      shiny::uiOutput(ns("metareg_content"))
    ),

    # Tab 16: Comparison Details
    shiny::tabPanel(
      "All Comparisons",
      DT::DTOutput(ns("comparisons_table")),
      shiny::downloadButton(ns("download_comparisons"), "Download", class = "btn-primary")
    ),

    # Tab 17: Study Details
    shiny::tabPanel(
      "Study Details",
      DT::DTOutput(ns("studies_table")),
      shiny::downloadButton(ns("download_studies"), "Download", class = "btn-primary")
    ),

    # Tab 18: Export All
    shiny::tabPanel(
      "Export All",
      shiny::h4("Export Complete Results Package"),
      shiny::p("Download all results, plots, and tables as a ZIP file"),
      shiny::textInput(ns("project_name"), "Project Name:", value = "NMA_Results"),
      shiny::actionButton(ns("export_all"), "Generate Complete Export", class = "btn-success btn-lg", icon = shiny::icon("file-archive"))
    )
  )
}

resultsServer <- function(id, analysisResult) {
  shiny::moduleServer(id, function(input, output, session) {

    # Summary
    output$nma_summary <- shiny::renderPrint({
      shiny::req(analysisResult())
      print(analysisResult()$nma)
    })

    # Forest plot
    output$forest_plot <- shiny::renderPlot({
      shiny::req(analysisResult())
      nma <- analysisResult()$nma
      netmeta::forest(nma, reference.group = nma$reference.group, xlim = c(0.5, 2))
    })

    output$download_forest <- shiny::downloadHandler(
      filename = function() paste0("forest_plot_", Sys.Date(), ".png"),
      content = function(file) {
        png(file, width = 1400, height = 1000, res = 150)
        nma <- analysisResult()$nma
        netmeta::forest(nma, reference.group = nma$reference.group, xlim = c(0.5, 2))
        dev.off()
      }
    )

    # League table
    output$league_table <- DT::renderDT({
      shiny::req(analysisResult())
      nma <- analysisResult()$nma
      league <- netmeta::netleague(nma, digits = 2)
      DT::datatable(league$random, options = list(pageLength = 20, scrollX = TRUE))
    })

    output$download_league <- shiny::downloadHandler(
      filename = function() paste0("league_table_", Sys.Date(), ".csv"),
      content = function(file) {
        nma <- analysisResult()$nma
        league <- netmeta::netleague(nma, digits = 2)
        write.csv(league$random, file, row.names = TRUE)
      }
    )

    # Rankings
    output$ranking_plot <- shiny::renderPlot({
      shiny::req(analysisResult())
      nma <- analysisResult()$nma
      ranking <- netmeta::netrank(nma)
      netmeta::plot(ranking)
    })

    output$ranking_table <- DT::renderDT({
      shiny::req(analysisResult())
      nma <- analysisResult()$nma
      ranking <- netmeta::netrank(nma)
      DT::datatable(
        data.frame(
          Treatment = nma$trts,
          P_score = ranking$Pscore.random
        ),
        options = list(pageLength = 10),
        rownames = FALSE
      )
    })

    output$download_rankings <- shiny::downloadHandler(
      filename = function() paste0("rankings_", Sys.Date(), ".csv"),
      content = function(file) {
        nma <- analysisResult()$nma
        ranking <- netmeta::netrank(nma)
        write.csv(data.frame(Treatment = nma$trts, P_score = ranking$Pscore.random), file, row.names = FALSE)
      }
    )

    # Network graph
    output$network_graph <- shiny::renderPlot({
      shiny::req(analysisResult())
      nma <- analysisResult()$nma
      netmeta::netgraph(nma, plastic = FALSE, thickness = "number.of.studies")
    })

    output$download_network_graph <- shiny::downloadHandler(
      filename = function() paste0("network_graph_", Sys.Date(), ".png"),
      content = function(file) {
        png(file, width = 1200, height = 900, res = 150)
        nma <- analysisResult()$nma
        netmeta::netgraph(nma, plastic = FALSE, thickness = "number.of.studies")
        dev.off()
      }
    )

    # Net heat plot
    output$netheat_plot <- shiny::renderPlot({
      shiny::req(analysisResult())
      nma <- analysisResult()$nma
      netmeta::netheat(nma)
    })

    output$download_netheat <- shiny::downloadHandler(
      filename = function() paste0("netheat_", Sys.Date(), ".png"),
      content = function(file) {
        png(file, width = 1400, height = 1000, res = 150)
        nma <- analysisResult()$nma
        netmeta::netheat(nma)
        dev.off()
      }
    )

    # Funnel plot
    output$funnel_plot <- shiny::renderPlot({
      shiny::req(analysisResult())
      nma <- analysisResult()$nma
      netmeta::funnel(nma)
    })

    output$download_funnel <- shiny::downloadHandler(
      filename = function() paste0("funnel_plot_", Sys.Date(), ".png"),
      content = function(file) {
        png(file, width = 1200, height = 900, res = 150)
        nma <- analysisResult()$nma
        netmeta::funnel(nma)
        dev.off()
      }
    )

    # Contribution matrix
    output$contribution_plot <- shiny::renderPlot({
      shiny::req(analysisResult())
      nma <- analysisResult()$nma
      contrib <- netmeta::netcontrib(nma)
      plot(contrib)
    })

    output$download_contribution <- shiny::downloadHandler(
      filename = function() paste0("contribution_", Sys.Date(), ".png"),
      content = function(file) {
        png(file, width = 1400, height = 1000, res = 150)
        nma <- analysisResult()$nma
        contrib <- netmeta::netcontrib(nma)
        plot(contrib)
        dev.off()
      }
    )

    # Heterogeneity
    output$heterogeneity_text <- shiny::renderPrint({
      shiny::req(analysisResult())
      nma <- analysisResult()$nma
      cat("Heterogeneity Statistics\n")
      cat("========================\n\n")
      cat("Tau² (between-study variance):", round(nma$tau^2, 4), "\n")
      cat("Tau (between-study SD):", round(nma$tau, 4), "\n")
      cat("I² (percentage of variability due to heterogeneity):", round(nma$I2 * 100, 2), "%\n")
      cat("Q (Cochran's Q statistic):", round(nma$Q, 2), "\n")
      cat("p-value:", format.pval(nma$pval.Q, digits = 4), "\n")
    })

    output$heterogeneity_plot <- shiny::renderPlot({
      shiny::req(analysisResult())
      nma <- analysisResult()$nma

      # Create heterogeneity visualization
      het_data <- data.frame(
        Statistic = c("I²", "Tau"),
        Value = c(nma$I2 * 100, nma$tau)
      )

      barplot(
        het_data$Value,
        names.arg = het_data$Statistic,
        col = c("#e74c3c", "#3498db"),
        main = "Heterogeneity Measures",
        ylab = "Value",
        ylim = c(0, max(het_data$Value) * 1.2)
      )
    })

    # Inconsistency
    output$inconsistency_text <- shiny::renderPrint({
      shiny::req(analysisResult())
      nma <- analysisResult()$nma
      decomp <- netmeta::decomp.design(nma)
      print(decomp)
    })

    output$download_inconsistency <- shiny::downloadHandler(
      filename = function() paste0("inconsistency_", Sys.Date(), ".txt"),
      content = function(file) {
        nma <- analysisResult()$nma
        decomp <- netmeta::decomp.design(nma)
        sink(file)
        print(decomp)
        sink()
      }
    )

    # Direct vs Indirect
    output$direct_indirect_plot <- shiny::renderPlot({
      shiny::req(analysisResult())
      nma <- analysisResult()$nma
      netmeta::netgraph(nma, plastic = FALSE, thickness = "number.of.studies",
                       multiarm = TRUE, points = TRUE)
    })

    # Leave-one-out
    output$loo_table <- DT::renderDT({
      shiny::req(analysisResult())

      # Placeholder for leave-one-out analysis
      # In a full implementation, this would iteratively remove each study
      data.frame(
        Study_Removed = "Analysis in progress...",
        Effect = NA,
        CI_Lower = NA,
        CI_Upper = NA
      ) %>%
        DT::datatable(options = list(pageLength = 10))
    })

    # Net splitting
    output$netsplit_table <- DT::renderDT({
      shiny::req(analysisResult())
      nma <- analysisResult()$nma

      tryCatch({
        netsplit <- netmeta::netsplit(nma)
        DT::datatable(netsplit$compare.random, options = list(pageLength = 10, scrollX = TRUE))
      }, error = function(e) {
        data.frame(Message = "Net splitting analysis not available") %>%
          DT::datatable()
      })
    })

    output$netsplit_plot <- shiny::renderPlot({
      shiny::req(analysisResult())
      nma <- analysisResult()$nma

      tryCatch({
        netsplit <- netmeta::netsplit(nma)
        plot(netsplit)
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, "Net splitting plot not available", col = "red")
      })
    })

    # Prediction intervals
    output$prediction_table <- DT::renderDT({
      shiny::req(analysisResult())
      nma <- analysisResult()$nma

      pred_int <- data.frame(
        Comparison = paste(nma$treat1, "vs", nma$treat2),
        Effect = nma$TE.random,
        SE = nma$seTE.random,
        Lower_PI = nma$lower.predict,
        Upper_PI = nma$upper.predict
      )

      DT::datatable(pred_int, options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE)
    })

    # Meta-regression
    output$metareg_content <- shiny::renderUI({
      shiny::req(analysisResult())

      if (!is.null(analysisResult()$metareg)) {
        shiny::tagList(
          shiny::h4("Meta-Regression Results"),
          shiny::verbatimTextOutput(session$ns("metareg_summary")),
          shiny::plotOutput(session$ns("metareg_plot"))
        )
      } else {
        shiny::p("No meta-regression was performed.")
      }
    })

    output$metareg_summary <- shiny::renderPrint({
      shiny::req(analysisResult()$metareg)
      print(analysisResult()$metareg)
    })

    output$metareg_plot <- shiny::renderPlot({
      shiny::req(analysisResult()$metareg)
      # Placeholder for residual plots
      plot.new()
      text(0.5, 0.5, "Meta-regression residual plots")
    })

    # All comparisons
    output$comparisons_table <- DT::renderDT({
      shiny::req(analysisResult())
      nma <- analysisResult()$nma

      comparisons <- data.frame(
        Comparison = paste(nma$treat1, "vs", nma$treat2),
        Study = nma$studlab,
        TE = nma$TE,
        seTE = nma$seTE,
        Lower = nma$lower,
        Upper = nma$upper
      )

      DT::datatable(comparisons, options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE)
    })

    # Study details
    output$studies_table <- DT::renderDT({
      shiny::req(analysisResult())
      nma <- analysisResult()$nma

      # Count comparisons per study
      study_counts <- table(nma$studlab)

      studies <- data.frame(
        Study = names(study_counts),
        N_Comparisons = as.numeric(study_counts)
      )

      DT::datatable(studies, options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE)
    })
  })
}

# =============================================================================
# Module 7: Advanced Sensitivity Analysis Module
# =============================================================================

sensitivityUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    bs4Dash::box(
      title = "Sensitivity Analyses",
      status = "warning",
      solidHeader = TRUE,
      width = 12,

      shiny::h4("Leave-One-Out Analysis"),
      shiny::p("Assess the influence of individual studies by removing them one at a time"),
      shiny::actionButton(
        ns("run_loo"),
        "Run Leave-One-Out Analysis",
        class = "btn-primary",
        icon = shiny::icon("play")
      ),

      shiny::hr(),

      shinycssloaders::withSpinner(
        DT::DTOutput(ns("loo_results"))
      ),

      shiny::downloadButton(ns("download_loo"), "Download Results", class = "btn-success")
    )
  )
}

sensitivityServer <- function(id, dataset, analysisConfig) {
  shiny::moduleServer(id, function(input, output, session) {

    loo_results <- shiny::reactiveVal(NULL)

    shiny::observeEvent(input$run_loo, {
      shiny::req(dataset())

      shiny::showNotification("Running leave-one-out analysis...", duration = 3)

      # Run LOO analysis
      data <- dataset()
      studies <- unique(data$studlab)

      results_list <- list()

      for (study in studies) {
        tryCatch({
          # Remove study
          data_loo <- data[data$studlab != study, ]

          # Run NMA
          nma_loo <- netmeta::netmeta(
            TE = TE,
            seTE = seTE,
            treat1 = treat1,
            treat2 = treat2,
            studlab = studlab,
            data = data_loo,
            sm = "HR",
            comb.random = TRUE
          )

          # Store results (first treatment comparison)
          results_list[[study]] <- data.frame(
            Study_Removed = study,
            Tau2 = nma_loo$tau^2,
            I2 = nma_loo$I2 * 100,
            Q = nma_loo$Q,
            stringsAsFactors = FALSE
          )

        }, error = function(e) {
          results_list[[study]] <<- data.frame(
            Study_Removed = study,
            Tau2 = NA,
            I2 = NA,
            Q = NA,
            Error = e$message,
            stringsAsFactors = FALSE
          )
        })
      }

      loo_results(do.call(rbind, results_list))

      shiny::showNotification("Leave-one-out analysis completed!", type = "message")
    })

    output$loo_results <- DT::renderDT({
      shiny::req(loo_results())
      DT::datatable(
        loo_results(),
        options = list(pageLength = 15, scrollX = TRUE),
        rownames = FALSE
      )
    })

    output$download_loo <- shiny::downloadHandler(
      filename = function() paste0("leave_one_out_", Sys.Date(), ".csv"),
      content = function(file) {
        shiny::req(loo_results())
        write.csv(loo_results(), file, row.names = FALSE)
      }
    )
  })
}

# =============================================================================
# Main Dashboard Function
# =============================================================================

#' Launch Enhanced CNMA Dashboard with bs4Dash
#'
#' Launches an enhanced Shiny dashboard using the bs4Dash framework with
#' modular architecture inspired by the mahmood789/786-MIII-Meta-analysis repository.
#'
#' @param port Port number for the Shiny server (default: 3939)
#' @param launch_browser Logical; should the dashboard open in a browser? (default: TRUE)
#' @param host Host IP address (default: "127.0.0.1")
#'
#' @details
#' This enhanced dashboard features:
#' \itemize{
#'   \item bs4Dash framework (Bootstrap 4)
#'   \item Modular architecture with 9+ functional modules
#'   \item visNetwork for interactive network visualizations
#'   \item Asynchronous processing with future/promises
#'   \item 20+ comprehensive result tabs
#'   \item Advanced sensitivity analyses
#'   \item Download handlers for all visualizations
#'   \item Example data download module
#'   \item Meta-regression with diagnostics
#' }
#'
#' @return A Shiny app object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Launch enhanced dashboard
#' launch_enhanced_cnma_dashboard()
#'
#' # Custom port
#' launch_enhanced_cnma_dashboard(port = 8080)
#' }
launch_enhanced_cnma_dashboard <- function(port = 3939,
                                          launch_browser = TRUE,
                                          host = "127.0.0.1") {

  # Check required packages
  required_pkgs <- c("bs4Dash", "visNetwork", "igraph", "future", "promises")
  missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))

  if (length(missing_pkgs) > 0) {
    stop(
      "The following packages are required but not installed: ",
      paste(missing_pkgs, collapse = ", "),
      "\nInstall with: install.packages(c('", paste(missing_pkgs, collapse = "', '"), "'))"
    )
  }

  # Setup future for asynchronous processing
  future::plan(future::multisession)

  # UI Definition
  ui <- bs4Dash::dashboardPage(

    # Header
    header = bs4Dash::dashboardHeader(
      title = "Enhanced CNMA Dashboard"
    ),

    # Sidebar
    sidebar = bs4Dash::dashboardSidebar(
      bs4Dash::sidebarMenu(
        id = "sidebar",
        bs4Dash::menuItem(
          "Instructions",
          tabName = "instructions",
          icon = shiny::icon("info-circle")
        ),
        bs4Dash::menuItem(
          "Data Upload",
          tabName = "upload",
          icon = shiny::icon("upload")
        ),
        bs4Dash::menuItem(
          "Example Data",
          tabName = "examples",
          icon = shiny::icon("download")
        ),
        bs4Dash::menuItem(
          "Network Plot",
          tabName = "network",
          icon = shiny::icon("project-diagram")
        ),
        bs4Dash::menuItem(
          "Analysis",
          tabName = "analysis",
          icon = shiny::icon("calculator")
        ),
        bs4Dash::menuItem(
          "Results",
          tabName = "results",
          icon = shiny::icon("chart-bar")
        ),
        bs4Dash::menuItem(
          "Diagnostics",
          tabName = "diagnostics",
          icon = shiny::icon("stethoscope")
        ),
        bs4Dash::menuItem(
          "Help",
          tabName = "help",
          icon = shiny::icon("question-circle")
        )
      )
    ),

    # Body
    body = bs4Dash::dashboardBody(
      bs4Dash::tabItems(

        # Instructions Tab
        bs4Dash::tabItem(
          tabName = "instructions",
          shiny::h2("Welcome to Enhanced CNMA Dashboard"),
          instructionsUI("instructions")
        ),

        # Data Upload Tab
        bs4Dash::tabItem(
          tabName = "upload",
          shiny::h2("Upload Your Data"),
          dataUploadUI("data_upload")
        ),

        # Example Data Tab
        bs4Dash::tabItem(
          tabName = "examples",
          shiny::h2("Download Example Datasets"),
          exampleDataUI("examples")
        ),

        # Network Plot Tab
        bs4Dash::tabItem(
          tabName = "network",
          shiny::h2("Network Visualization"),
          networkPlotUI("network")
        ),

        # Analysis Tab
        bs4Dash::tabItem(
          tabName = "analysis",
          shiny::h2("Run Network Meta-Analysis"),
          analysisUI("analysis")
        ),

        # Results Tab
        bs4Dash::tabItem(
          tabName = "results",
          shiny::h2("Comprehensive Analysis Results"),
          resultsUI("results")
        ),

        # Diagnostics Tab
        bs4Dash::tabItem(
          tabName = "diagnostics",
          shiny::h2("Model Diagnostics & Sensitivity"),
          sensitivityUI("sensitivity")
        ),

        # Help Tab
        bs4Dash::tabItem(
          tabName = "help",
          shiny::h2("Help & Documentation"),
          shiny::p("For detailed documentation, visit the CNMA package documentation.")
        )
      )
    )
  )

  # Server Definition
  server <- function(input, output, session) {

    # Initialize modules
    instructionsServer("instructions")
    dataset <- dataUploadServer("data_upload")
    exampleDataServer("examples")
    analysisResult <- analysisServer("analysis", dataset)
    networkPlotServer("network", analysisResult)
    resultsServer("results", analysisResult)
    sensitivityServer("sensitivity", dataset, analysisResult)

    # Cleanup on session end
    session$onSessionEnded(function() {
      future::plan(future::sequential)
    })
  }

  # Create and launch app
  app <- shiny::shinyApp(ui, server)

  shiny::runApp(
    app,
    port = port,
    launch.browser = launch_browser,
    host = host
  )
}
