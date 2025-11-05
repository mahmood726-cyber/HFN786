# =========================================================
# Comprehensive Shiny Dashboard for CNMA
# Interactive web application for network meta-analysis
# =========================================================

#' Launch Comprehensive CNMA Dashboard
#'
#' Launches interactive Shiny dashboard with complete NMA workflow including
#' data upload, analysis, visualization, validation, and manuscript generation.
#'
#' @param port Port number for Shiny app (default 3838)
#' @param launch_browser Launch in browser (default TRUE)
#' @param host Host address (default "127.0.0.1")
#' @return Shiny app object
#' @export
#' @examples
#' \dontrun{
#' # Launch comprehensive dashboard
#' launch_cnma_dashboard()
#'
#' # Launch on specific port
#' launch_cnma_dashboard(port = 8080)
#' }
launch_cnma_dashboard <- function(port = 3838,
                                  launch_browser = TRUE,
                                  host = "127.0.0.1") {

  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' required. Install with: install.packages('shiny')")
  }

  if (!requireNamespace("shinydashboard", quietly = TRUE)) {
    stop("Package 'shinydashboard' required. Install with: install.packages('shinydashboard')")
  }

  msg("Launching CNMA Comprehensive Dashboard...")
  msg("  Port: %d", port)
  msg("  Access at: http://%s:%d", host, port)

  app <- .create_comprehensive_dashboard_app()

  shiny::runApp(
    app,
    port = port,
    host = host,
    launch.browser = launch_browser
  )
}

#' @keywords internal
.create_comprehensive_dashboard_app <- function() {

  # UI
  ui <- shinydashboard::dashboardPage(
    skin = "blue",

    # Header
    shinydashboard::dashboardHeader(
      title = "CNMA: AI-Powered Network Meta-Analysis",
      titleWidth = 400
    ),

    # Sidebar
    shinydashboard::dashboardSidebar(
      width = 250,
      shinydashboard::sidebarMenu(
        id = "sidebar",
        shinydashboard::menuItem("Dashboard", tabName = "dashboard", icon = shiny::icon("dashboard")),
        shinydashboard::menuItem("Data Upload", tabName = "upload", icon = shiny::icon("upload")),
        shinydashboard::menuItem("Analysis", tabName = "analysis", icon = shiny::icon("chart-line")),
        shinydashboard::menuItem("Visualizations", tabName = "visualizations", icon = shiny::icon("images")),
        shinydashboard::menuItem("Results", tabName = "results", icon = shiny::icon("table")),
        shinydashboard::menuItem("Rules Validation", tabName = "validation", icon = shiny::icon("check-circle")),
        shinydashboard::menuItem("Manuscript", tabName = "manuscript", icon = shiny::icon("file-alt")),
        shinydashboard::menuItem("AI Assistant", tabName = "ai", icon = shiny::icon("robot")),
        shinydashboard::menuItem("Settings", tabName = "settings", icon = shiny::icon("cog")),
        shinydashboard::menuItem("Help", tabName = "help", icon = shiny::icon("question-circle"))
      )
    ),

    # Body
    shinydashboard::dashboardBody(
      # Custom CSS
      shiny::tags$head(
        shiny::tags$style(shiny::HTML("
          .box-header { font-weight: bold; }
          .info-box { margin-bottom: 15px; }
          .progress-bar { background-color: #3c8dbc; }
          .btn-primary { background-color: #3c8dbc; border-color: #367fa9; }
          .btn-success { background-color: #00a65a; border-color: #008d4c; }
        "))
      ),

      shinydashboard::tabItems(
        # Dashboard Tab
        shinydashboard::tabItem(
          tabName = "dashboard",
          shiny::fluidRow(
            shinydashboard::infoBox(
              "Studies", shiny::textOutput("n_studies"), icon = shiny::icon("book"),
              color = "blue", width = 3
            ),
            shinydashboard::infoBox(
              "Treatments", shiny::textOutput("n_treatments"), icon = shiny::icon("pills"),
              color = "green", width = 3
            ),
            shinydashboard::infoBox(
              "Comparisons", shiny::textOutput("n_comparisons"), icon = shiny::icon("exchange-alt"),
              color = "yellow", width = 3
            ),
            shinydashboard::infoBox(
              "Participants", shiny::textOutput("n_participants"), icon = shiny::icon("users"),
              color = "red", width = 3
            )
          ),

          shiny::fluidRow(
            shinydashboard::box(
              title = "Network Visualization", status = "primary", solidHeader = TRUE,
              width = 8, height = 500,
              plotly::plotlyOutput("network_plot", height = "450px")
            ),
            shinydashboard::box(
              title = "Quick Stats", status = "info", solidHeader = TRUE,
              width = 4,
              shiny::verbatimTextOutput("quick_stats")
            )
          ),

          shiny::fluidRow(
            shinydashboard::box(
              title = "Analysis Progress", status = "success", solidHeader = TRUE,
              width = 12,
              shiny::uiOutput("analysis_progress")
            )
          )
        ),

        # Data Upload Tab
        shinydashboard::tabItem(
          tabName = "upload",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Upload Data", status = "primary", solidHeader = TRUE,
              width = 12,
              shiny::fileInput("data_file", "Choose CSV File",
                             accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
              shiny::checkboxInput("header", "Header", TRUE),
              shiny::radioButtons("sep", "Separator",
                                choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),
                                selected = ","),
              shiny::actionButton("load_example", "Load Example Data", class = "btn-info"),
              shiny::hr(),
              DT::dataTableOutput("data_preview")
            )
          ),

          shiny::fluidRow(
            shinydashboard::box(
              title = "Data Validation", status = "warning", solidHeader = TRUE,
              width = 12,
              shiny::verbatimTextOutput("data_validation_output")
            )
          )
        ),

        # Analysis Tab
        shinydashboard::tabItem(
          tabName = "analysis",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Analysis Configuration", status = "primary", solidHeader = TRUE,
              width = 4,
              shiny::selectInput("sm", "Summary Measure",
                               choices = c("HR" = "HR", "OR" = "OR", "RR" = "RR",
                                         "MD" = "MD", "SMD" = "SMD"),
                               selected = "HR"),
              shiny::selectInput("model", "Model Type",
                               choices = c("Random Effects" = "random",
                                         "Fixed Effect" = "fixed",
                                         "Both" = "both"),
                               selected = "random"),
              shiny::selectInput("ref_treatment", "Reference Treatment",
                               choices = NULL),
              shiny::checkboxInput("use_ai", "Use AI Enhancement", TRUE),
              shiny::checkboxInput("validate_rules", "Run Rules Validation", TRUE),
              shiny::actionButton("run_analysis", "Run Analysis", class = "btn-success btn-lg"),
              shiny::hr(),
              shiny::downloadButton("download_results", "Download Results", class = "btn-primary")
            ),

            shinydashboard::box(
              title = "Analysis Log", status = "info", solidHeader = TRUE,
              width = 8,
              shiny::verbatimTextOutput("analysis_log", placeholder = TRUE)
            )
          ),

          shiny::fluidRow(
            shinydashboard::box(
              title = "Treatment Rankings", status = "success", solidHeader = TRUE,
              width = 6,
              DT::dataTableOutput("rankings_table")
            ),

            shinydashboard::box(
              title = "Network Statistics", status = "warning", solidHeader = TRUE,
              width = 6,
              shiny::verbatimTextOutput("network_stats")
            )
          )
        ),

        # Visualizations Tab
        shinydashboard::tabItem(
          tabName = "visualizations",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Visualization Options", status = "primary", solidHeader = TRUE,
              width = 12,
              shiny::selectInput("viz_type", "Select Visualization",
                               choices = c(
                                 "Network Plot" = "network",
                                 "Forest Plot" = "forest",
                                 "Contribution Heatmap" = "contribution",
                                 "Effect Matrix" = "effect_matrix",
                                 "Harvest Plot" = "harvest",
                                 "Ranking Heatmap" = "ranking",
                                 "Temporal Trends" = "temporal",
                                 "Risk of Bias" = "rob",
                                 "3D Network" = "network_3d",
                                 "Evidence Gaps" = "gaps"
                               ),
                               selected = "network"),
              shiny::checkboxInput("interactive_viz", "Make Interactive", TRUE),
              shiny::actionButton("generate_viz", "Generate Visualization", class = "btn-success")
            )
          ),

          shiny::fluidRow(
            shinydashboard::box(
              title = "Visualization Output", status = "info", solidHeader = TRUE,
              width = 12, height = 600,
              plotly::plotlyOutput("viz_output", height = "550px")
            )
          )
        ),

        # Results Tab
        shinydashboard::tabItem(
          tabName = "results",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Treatment Effects", status = "primary", solidHeader = TRUE,
              width = 12,
              DT::dataTableOutput("effects_table")
            )
          ),

          shiny::fluidRow(
            shinydashboard::box(
              title = "Heterogeneity & Inconsistency", status = "warning", solidHeader = TRUE,
              width = 6,
              shiny::verbatimTextOutput("heterogeneity_output")
            ),

            shinydashboard::box(
              title = "Publication Bias", status = "danger", solidHeader = TRUE,
              width = 6,
              plotly::plotlyOutput("funnel_plot", height = "300px")
            )
          )
        ),

        # Rules Validation Tab
        shinydashboard::tabItem(
          tabName = "validation",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Rules Validation Dashboard", status = "primary", solidHeader = TRUE,
              width = 12,
              shiny::actionButton("run_validation", "Run Validation", class = "btn-success"),
              shiny::hr(),
              shiny::uiOutput("validation_summary")
            )
          ),

          shiny::fluidRow(
            shinydashboard::box(
              title = "Violations", status = "danger", solidHeader = TRUE,
              width = 12,
              DT::dataTableOutput("violations_table")
            )
          )
        ),

        # Manuscript Tab
        shinydashboard::tabItem(
          tabName = "manuscript",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Manuscript Generation", status = "primary", solidHeader = TRUE,
              width = 4,
              shiny::selectInput("journal_style", "Journal Style",
                               choices = c("BMJ" = "BMJ", "Lancet" = "Lancet",
                                         "JAMA" = "JAMA", "Generic" = "generic"),
                               selected = "BMJ"),
              shiny::numericInput("methods_word_limit", "Methods Word Limit", 500, min = 100, max = 2000),
              shiny::numericInput("results_word_limit", "Results Word Limit", 800, min = 100, max = 3000),
              shiny::actionButton("generate_methods", "Generate Methods", class = "btn-info"),
              shiny::actionButton("generate_results", "Generate Results", class = "btn-info"),
              shiny::actionButton("generate_complete", "Generate Complete Manuscript", class = "btn-success"),
              shiny::hr(),
              shiny::downloadButton("download_manuscript", "Download Manuscript")
            ),

            shinydashboard::box(
              title = "Manuscript Preview", status = "info", solidHeader = TRUE,
              width = 8,
              shiny::tabsetPanel(
                shiny::tabPanel("Methods", shiny::verbatimTextOutput("methods_text")),
                shiny::tabPanel("Results", shiny::verbatimTextOutput("results_text")),
                shiny::tabPanel("Compliance", shiny::verbatimTextOutput("compliance_text"))
              )
            )
          )
        ),

        # AI Assistant Tab
        shinydashboard::tabItem(
          tabName = "ai",
          shiny::fluidRow(
            shinydashboard::box(
              title = "AI Assistant Configuration", status = "primary", solidHeader = TRUE,
              width = 12,
              shiny::selectInput("ai_model", "LLama 3 Model",
                               choices = c("llama3" = "llama3", "llama3:70b" = "llama3:70b"),
                               selected = "llama3"),
              shiny::actionButton("configure_ai", "Configure AI", class = "btn-success"),
              shiny::textInput("ai_query", "Ask AI Assistant:", placeholder = "E.g., Interpret my results..."),
              shiny::actionButton("ask_ai", "Ask", class = "btn-primary"),
              shiny::hr(),
              shiny::verbatimTextOutput("ai_response")
            )
          )
        ),

        # Settings Tab
        shinydashboard::tabItem(
          tabName = "settings",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Application Settings", status = "primary", solidHeader = TRUE,
              width = 12,
              shiny::sliderInput("figure_width", "Figure Width", min = 400, max = 1600, value = 800),
              shiny::sliderInput("figure_height", "Figure Height", min = 300, max = 1200, value = 600),
              shiny::selectInput("theme", "Dashboard Theme",
                               choices = c("Blue" = "blue", "Black" = "black",
                                         "Purple" = "purple", "Green" = "green"),
                               selected = "blue"),
              shiny::actionButton("apply_settings", "Apply Settings", class = "btn-success")
            )
          )
        ),

        # Help Tab
        shinydashboard::tabItem(
          tabName = "help",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Quick Start Guide", status = "info", solidHeader = TRUE,
              width = 12,
              shiny::HTML("
                <h3>Welcome to CNMA Dashboard!</h3>
                <ol>
                  <li><strong>Upload Data:</strong> Go to 'Data Upload' tab and load your CSV file with columns: studlab, treat1, treat2, TE, seTE</li>
                  <li><strong>Configure Analysis:</strong> In 'Analysis' tab, select your settings and reference treatment</li>
                  <li><strong>Run Analysis:</strong> Click 'Run Analysis' to perform network meta-analysis</li>
                  <li><strong>Explore Results:</strong> View visualizations, rankings, and treatment effects</li>
                  <li><strong>Validate:</strong> Check 'Rules Validation' for quality assurance</li>
                  <li><strong>Generate Manuscript:</strong> Use 'Manuscript' tab to create publication-ready text</li>
                  <li><strong>AI Assistance:</strong> Get intelligent insights from local LLama 3</li>
                </ol>
                <hr>
                <h4>Example Data Format:</h4>
                <table class='table table-striped'>
                  <tr><th>studlab</th><th>treat1</th><th>treat2</th><th>TE</th><th>seTE</th></tr>
                  <tr><td>Study1</td><td>Placebo</td><td>DrugA</td><td>-0.15</td><td>0.12</td></tr>
                  <tr><td>Study2</td><td>Placebo</td><td>DrugB</td><td>-0.28</td><td>0.15</td></tr>
                </table>
              ")
            )
          )
        )
      )
    )
  )

  # Server
  server <- function(input, output, session) {

    # Reactive values
    values <- shiny::reactiveValues(
      data = NULL,
      nma_results = NULL,
      validation_report = NULL,
      manuscript = NULL
    )

    # Load example data
    shiny::observeEvent(input$load_example, {
      values$data <- simulate_cnma_data(30, seed = 123)
      shiny::showNotification("Example data loaded successfully!", type = "message")
    })

    # File upload
    shiny::observeEvent(input$data_file, {
      req(input$data_file)
      values$data <- read.csv(input$data_file$datapath,
                             header = input$header,
                             sep = input$sep)
      shiny::showNotification("Data uploaded successfully!", type = "message")
    })

    # Data preview
    output$data_preview <- DT::renderDataTable({
      req(values$data)
      DT::datatable(values$data, options = list(pageLength = 10))
    })

    # Dashboard info boxes
    output$n_studies <- shiny::renderText({
      req(values$data)
      as.character(length(unique(values$data$studlab)))
    })

    output$n_treatments <- shiny::renderText({
      req(values$data)
      as.character(length(unique(c(values$data$treat1, values$data$treat2))))
    })

    output$n_comparisons <- shiny::renderText({
      req(values$data)
      as.character(nrow(values$data))
    })

    output$n_participants <- shiny::renderText({
      req(values$data)
      if ("n1" %in% names(values$data) && "n2" %in% names(values$data)) {
        as.character(sum(values$data$n1 + values$data$n2, na.rm = TRUE))
      } else {
        "N/A"
      }
    })

    # Update reference treatment choices
    shiny::observe({
      req(values$data)
      treatments <- unique(c(values$data$treat1, values$data$treat2))
      shiny::updateSelectInput(session, "ref_treatment", choices = treatments)
    })

    # Run analysis
    shiny::observeEvent(input$run_analysis, {
      req(values$data)

      shiny::withProgress(message = 'Running analysis...', value = 0, {
        shiny::incProgress(0.2, detail = "Validating data...")

        # Run NMA
        shiny::incProgress(0.4, detail = "Performing network meta-analysis...")
        values$nma_results <- run_cnma_analysis(
          values$data,
          ref_treatment = input$ref_treatment,
          config = setup_cnma(sm = input$sm)
        )

        # Validation
        if (input$validate_rules) {
          shiny::incProgress(0.2, detail = "Running rules validation...")
          values$validation_report <- run_rules_validation(values$data, values$nma_results)
        }

        shiny::incProgress(0.2, detail = "Complete!")
      })

      shiny::showNotification("Analysis completed successfully!", type = "message", duration = 5)
    })

    # Network plot
    output$network_plot <- plotly::renderPlotly({
      req(values$nma_results)
      plot_interactive_network(values$nma_results$results$main_nma)
    })

    # Rankings table
    output$rankings_table <- DT::renderDataTable({
      req(values$nma_results)
      rankings <- calculate_rankings(values$nma_results$results$main_nma)
      DT::datatable(rankings, options = list(pageLength = 10))
    })

    # Generate manuscript
    shiny::observeEvent(input$generate_complete, {
      req(values$nma_results)

      shiny::withProgress(message = 'Generating manuscript...', value = 0, {
        shiny::incProgress(0.5, detail = "Creating methods and results sections...")

        values$manuscript <- generate_complete_manuscript(
          values$nma_results,
          values$data,
          journal_style = input$journal_style
        )

        shiny::incProgress(0.5, detail = "Complete!")
      })

      shiny::showNotification("Manuscript generated!", type = "message")
    })

    # Methods text output
    output$methods_text <- shiny::renderText({
      req(values$manuscript)
      values$manuscript$methods$methods_text
    })

    # Results text output
    output$results_text <- shiny::renderText({
      req(values$manuscript)
      values$manuscript$results$results_text
    })
  }

  list(ui = ui, server = server)
}

#' Launch Quick Analysis Dashboard
#'
#' Simplified dashboard for quick NMA analysis.
#'
#' @param data Optional data frame to start with
#' @return Shiny app object
#' @export
launch_quick_analysis <- function(data = NULL) {

  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' required.")
  }

  msg("Launching Quick Analysis Dashboard...")

  # Simplified UI
  ui <- shiny::fluidPage(
    theme = NULL,
    shiny::titlePanel("CNMA: Quick Analysis"),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,
        shiny::fileInput("data_file", "Upload CSV"),
        shiny::actionButton("run", "Run Analysis", class = "btn-success btn-block"),
        shiny::hr(),
        shiny::downloadButton("download", "Download Results")
      ),

      shiny::mainPanel(
        width = 9,
        shiny::tabsetPanel(
          shiny::tabPanel("Network", plotly::plotlyOutput("network")),
          shiny::tabPanel("Rankings", DT::dataTableOutput("rankings")),
          shiny::tabPanel("Results", shiny::verbatimTextOutput("results"))
        )
      )
    )
  )

  server <- function(input, output, session) {
    # Simplified server logic
    values <- shiny::reactiveValues(data = data, results = NULL)

    # File upload
    shiny::observeEvent(input$data_file, {
      req(input$data_file)
      values$data <- read.csv(input$data_file$datapath)
    })

    # Run analysis
    shiny::observeEvent(input$run, {
      req(values$data)
      values$results <- run_cnma_analysis(values$data)
    })

    # Network plot
    output$network <- plotly::renderPlotly({
      req(values$results)
      plot_interactive_network(values$results$results$main_nma)
    })

    # Rankings
    output$rankings <- DT::renderDataTable({
      req(values$results)
      calculate_rankings(values$results$results$main_nma)
    })
  }

  shiny::shinyApp(ui, server)
}

#' Launch Visualization Explorer
#'
#' Interactive dashboard for exploring all visualization types.
#'
#' @param nma_results NMA results object
#' @param data Data frame
#' @return Shiny app object
#' @export
launch_visualization_explorer <- function(nma_results = NULL, data = NULL) {

  msg("Launching Visualization Explorer...")

  ui <- shiny::fluidPage(
    shiny::titlePanel("CNMA: Visualization Explorer"),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,
        shiny::selectInput("viz_type", "Visualization Type",
                         choices = c("Network", "Forest", "Contribution Heatmap",
                                   "Effect Matrix", "Harvest Plot", "3D Network",
                                   "Ranking Heatmap", "Temporal Trends",
                                   "Evidence Gaps", "Confidence Ellipses")),
        shiny::actionButton("generate", "Generate", class = "btn-primary"),
        shiny::hr(),
        shiny::downloadButton("download_viz", "Download Plot")
      ),

      shiny::mainPanel(
        width = 9,
        plotly::plotlyOutput("viz_display", height = "700px")
      )
    )
  )

  server <- function(input, output, session) {
    values <- shiny::reactiveValues(
      nma_results = nma_results,
      data = data,
      current_plot = NULL
    )

    # Generate visualization
    shiny::observeEvent(input$generate, {
      req(values$nma_results)

      values$current_plot <- switch(
        input$viz_type,
        "Network" = plot_interactive_network(values$nma_results),
        "Contribution Heatmap" = plot_contribution_heatmap(values$nma_results, TRUE),
        "Effect Matrix" = plot_effect_matrix(values$nma_results, TRUE, TRUE),
        "3D Network" = plot_network_3d(values$nma_results, values$data),
        plot_interactive_network(values$nma_results)
      )
    })

    output$viz_display <- plotly::renderPlotly({
      req(values$current_plot)
      values$current_plot
    })
  }

  shiny::shinyApp(ui, server)
}
