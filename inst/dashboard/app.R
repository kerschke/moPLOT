library(shiny)
library(shinyjs)

source("dashboardFunctions.R")

# Maximum Upload Size = 50MB
# Might need to be changed for larger files
options("shiny.maxRequestSize" = 50L * (1024L**2L))

ui <- fluidPage(
  useShinyjs(),
  
  h1(code("moPLOT"), "Landscape Explorer"),
  
  fluidRow(
    
    column(4L,
           tabsetPanel(
             tabPanel(
               "Select MOP",
               p(),
               wellPanel(
                 selectInput("benchmark_set", "Benchmark set", c("Select a benchmark set"="", benchmark_sets)),
                 selectInput("fn_name", "Function", c("Select a benchmark set first"="")),
                 uiOutput("fn_args")
               ),
               
               div(
                 h3("Generate Data"),
                 wellPanel(
                   numericInput("grid_size", "Resolution per dimension", 100, min=20, max=3000, step=1),
                   checkboxInput("compute_plot", "Compute PLOT and heatmap", TRUE),
                   checkboxInput("compute_cost_landscape", "Compute cost landscape", FALSE),
                   splitLayout(
                     actionButton("evaluate_grid", "Evaluate", style = "width: 100%"), # icon = icon('th')
                     downloadButton("download_data", "Download", style = "width: 100%")
                   )
                 ),
                 id = "evaluate_grid_panel"
               )
             ),
             tabPanel(
               "Upload Data",
               p(),
               wellPanel(
                 fileInput("upload_data", "Upload Data", accept = c(".Rds"))
               )
             ),
             type = "pills",
             id = "fn_select"
           ),
    ),
    
    column(8,
           tabsetPanel(
             tabPanel(
               "PLOT",
               plotly::plotlyOutput("plot", height = "500px"),
               value = "tab_plot"
             ),
             tabPanel(
               "Gradient Field Heatmap",
               plotly::plotlyOutput("heatmap", height = "500px"),
               value = "tab_heatmap"
             ),
             tabPanel(
               "Cost Landscape",
               plotly::plotlyOutput("cost_landscape", height = "500px"),
               value = "tab_cost_landscape"
             ),
             id = "tabset_plots"
           ),
           div(
             h3("Plot Options"),
             wellPanel(
               selectInput("space", "Space to plot", c("Decision Space" = "decision.space", "Objective Space" = "objective.space", "Decision + Objective Space" = "both")),
               div(
                 selectInput("three_d_approach", "3D approach", c("MRI Scan" = "scan", "Onion Layers" = "layers", "Nondominated" = "pareto")),
                 conditionalPanel("input.three_d_approach == 'scan'",
                                  selectInput("scan_direction", "Scan direction", c("x₁" = "x1", "x₂" = "x2", "x₃" = "x3"), selected = "x3")
                 ),
                 id = "three_d_only"
               )
             ),
             id = "plot_options"
           )
    )
  )
)

server <- function(input, output, session) {
  plot_data <- reactiveValues()
  
  # outputOptions(output, suspendWhenHidden = FALSE)
  
  reset_plots <- function() {
    hide("tabset_plots")
    hide("plot_options")
    
    hideTab("tabset_plots", "tab_plot")
    hideTab("tabset_plots", "tab_heatmap")
    hideTab("tabset_plots", "tab_cost_landscape")
    
    enable("compute_plot")
    enable("compute_cost_landscape")
    
    plot_data$design <- NULL
    plot_data$less <- NULL
    plot_data$domination_counts <- NULL
  }
  
  # Hide some parts per default
  
  hide("evaluate_grid_panel")
  hide("fn_name")
  
  reset_plots()
  
  # Reactive getters
  
  get_benchmark_set <- reactive({
    switch(
      input$benchmark_set,
      biobj_bbob = biobj_bbob_functions,
      dtlz = dtlz_functions,
      mindist = mindist_functions,
      mmf = mmf_functions,
      mop = mop_functions,
      mpm2 = mpm2_functions,
      zdt = zdt_functions,
      other = other_functions
    )
  })
  
  get_generator_fn <- reactive({
    get_benchmark_set()[[input$fn_name]]
  })
  
  get_default_args <- reactive({
    generator_fn <- get_generator_fn()
    
    if (is.null(generator_fn)) {
      list()
    } else {
      args <- formals(generator_fn)
      
      if (length(args) > 0) args else list()
    }
  })
  
  get_selected_args <- reactive({
    args <- lapply(get_args_names_ui(), function(argument_name) {
      arg <- input[[argument_name]]
      
      if (length(arg) > 0) {
        if (!is.numeric(arg)) {
          tryCatch(
            eval(parse(text = arg)),
            error = function(e) {
              tryCatch(
                parse(text = arg),
                error = function(e) NULL
              )
            })
        } else {
          arg
        }
      } else {
        NULL
      }
    })
    
    # Only return args if correctly connected to output
    if (all(sapply(args, is.null))) return(list())
    
    names(args) <- names(get_default_args())
    args <- args[!sapply(args, function(arg) (is.null(arg) || length(as.character(arg)) == 0))]
    
    args
  })
  
  get_args_names_ui <- reactive({
    args_names <- names(get_default_args())
    
    if (is.null(args_names)) list() else paste0("args_", args_names)
  })
  
  get_fn <- reactive({
    generator_fn <- get_generator_fn()
    args <- get_selected_args()
    
    req(generator_fn)
    req(args)
    
    tryCatch({
      do.call(generator_fn, args)
    }, error = function(e) {
      print(e)
      
      NULL
    })
  })
  
  # Observers that change the UI dynamically
  
  observe({
    fn <- get_fn()
    
    req(fn)
    
    hide("evaluate_grid_panel")
    
    if (smoof::getNumberOfParameters(fn) == 2) {
      updateSliderInput("grid_size", session = session, value = 200, min=50, max=3000, step=50)
      hide("three_d_only")
    } else {
      updateSliderInput("grid_size", session = session, value = 50, min=20, max=200, step=10)
      show("three_d_only")
    }
    
    show("evaluate_grid_panel")
  })
  
  observe({
    bench_set <- get_benchmark_set()
    
    updateSelectInput(session, "fn_name", choices = names(bench_set))
    
    if (is.null(bench_set) || input$benchmark_set == "biobj_bbob") {
      hide("fn_name")
    } else {
      show("fn_name")
    }
  })
  
  # Dynamically created view with function parameters
  
  output$fn_args = renderUI({
    args <- get_default_args()
    
    ui_inputs <- lapply(names(args), function(argument_name) {
      arg_char <- as.character(args[[argument_name]])
      
      if (length(arg_char) > 0) {
        if (is.symbol(args[[argument_name]])) {
          value <- arg_char
        } else if (is.numeric(args[[argument_name]])){
          value <- args[[argument_name]]
        } else {
          value <- deparse(args[[argument_name]])
        }
      } else {
        value <- 0
      }
      
      input_args <- list(
        inputId = paste0("args_", argument_name), 
        value = value,
        label = code(argument_name)
      )
      
      # General special cases
      
      if (argument_name == "dimensions") {
        input_args$label = "Dimensions"
        input_args$value = 2
        input_args$min = 2
        input_args$max = 3
      } else if (argument_name == "n.objectives") {
        input_args$label = "Objectives"
        input_args$value = 2
        input_args$min = 2
        input_args$max = 3
      }
      
      # Special cases for Bi-Obj BBOB
      
      if (input$benchmark_set == "biobj_bbob") {
        if (argument_name == "fid") {
          input_args$label = "Function ID"
          input_args$value = 1
          input_args$min = 1
          input_args$max = 55
        }
        
        if (argument_name == "iid") {
          input_args$label = "Instance ID"
          input_args$value = 1
          input_args$min = 1
        }
      }
      
      input_type <- ifelse(is.numeric(input_args$value), numericInput, textInput)
      
      do.call(input_type, input_args)
    })
    
    ui_inputs$cellArgs <- list(style = "width: 49%")
    
    do.call(flowLayout, ui_inputs)
  })
  
  # Generating data
  
  observeEvent(c(input$grid_size, get_fn()), {
    # reset stored plot data when resolution or function changes
    reset_plots()
  })
  
  observeEvent(input$evaluate_grid, {
    fn <- get_fn()
    req(fn)
    
    disable(selector = "button")
    
    if (is.null(plot_data$design)) {
      # generate design and evaluate objective space
      
      design <- moPLOT::generateDesign(fn, points.per.dimension = input$grid_size)
      design$obj.space <- calculateObjectiveValues(design$dec.space, fn, parallelize = T)
      
      plot_data$design <<- design
    }
    
    if (input$compute_plot && is.null(plot_data$less)) {
      design <- plot_data$design
      
      gradients <- computeGradientFieldGrid(design) #, prec.angle = 1
      
      divergence <- computeDivergenceGrid(gradients$multi.objective, design$dims, design$step.sizes)
      
      # Calculate locally efficient points
      plot_data$less <<- localEfficientSetSkeleton(design, gradients, divergence, integration="fast")
      
      if (smoof::getNumberOfParameters(fn) == 2) {
        # PLOT only in 2D for now
        showTab("tabset_plots", "tab_plot")
      }
      
      showTab("tabset_plots", "tab_heatmap")
      
      disable("compute_plot")
    }
    
    if (input$compute_cost_landscape && is.null(plot_data$domination_counts)) {
      nds <- ecr::doNondominatedSorting(t(plot_data$design$obj.space))
      plot_data$domination_counts <<- cbind(height = nds$dom.counter + 1)
      
      showTab("tabset_plots", "tab_cost_landscape")
      
      disable("compute_cost_landscape")
    }
    
    show("tabset_plots")
    show("plot_options")
    
    enable(selector = "button")
  })
  
  # Plotting-related functions
  
  get_plot = function(plot_type, space, three_d_approach) {
    fn <- get_fn()
    
    req(fn)
    req(plot_data$design)
    
    grid <- plot_data$design
    
    grid$height <- switch (
      plot_type,
      heatmap = {
        req(plot_data$less)
        plot_data$less$height
      },
      cost_landscape = {
        req(plot_data$domination_counts)
        plot_data$domination_counts
      },
      PLOT = {
        req(plot_data$less)
        plot_data$less$height
      }
    )
    
    less <- plot_data$less
    
    d = smoof::getNumberOfParameters(fn)
    
    if (d == 2) {
      switch (plot_type,
              heatmap = plotly2DHeatmap(grid, fn, mode = space),
              cost_landscape = plotly2DHeatmap(grid, fn, mode = space),
              PLOT = plotly2DPLOT(grid$dec.space, grid$obj.space, less$sinks, less$height, fn, mode = space),
              NULL # if plot_type is invalid
      )
    } else if (d == 3) {
      switch (three_d_approach,
              pareto = plotly3DPareto(grid, fn, mode = space),
              layers = plotly3DLayers(grid, fn, mode = space),
              scan = plotly3DScan(grid, fn, mode = space, frame = input$scan_direction),
              NULL # if plot_type is invalid
      )
    } else {
      NULL
    }
  }
  
  output$plot = plotly::renderPlotly(
    reactive({
      disable("plot")
      print("Updating PLOT")
      
      p <- get_plot("PLOT", input$space, input$three_d_approach)
      enable("plot")
      p
    }, quoted = TRUE)()
  )
  
  output$heatmap = plotly::renderPlotly({
    reactive({
      disable("heatmap")
      print("Updating Heatmap")
      
      p <- get_plot("heatmap", input$space, input$three_d_approach)
      enable("heatmap")
      p
    }, quoted = TRUE)()
  })
  
  output$cost_landscape = plotly::renderPlotly({
    reactive({
      disable("cost_landscape")
      print("Updating Cost Landscape")
      
      p <- get_plot("cost_landscape", input$space, input$three_d_approach)
      enable("cost_landscape")
      p
    }, quoted = TRUE)()
  })
  
  # Up- and Download
  
  output$download_data = downloadHandler(
    filename = function() {
      fn = get_fn()
      req(fn)
      
      paste0(smoof::getName(fn), "-data-", input$grid_size, '.Rds')
    },
    
    content = function(file) {
      input_list <- reactiveValuesToList(input)
      plot_data_list <- reactiveValuesToList(plot_data)
      
      data = list()
      
      data$input <- input_list
      data$plot_data <- plot_data_list
      
      saveRDS(data, file)
    }
  )
  
  observeEvent(input$upload_data, {
    if (!endsWith(input$upload_data$datapath, ".Rds")) {
      print(input$upload_data)
      return()
    }
    
    data <- readRDS(input$upload_data$datapath)
    
    reset_plots()
    
    updateTabsetPanel(session, "fn_select", selected = "Select MOP")
    updateSelectInput(session, "benchmark_set", selected = data$input$benchmark_set)

    updateNumericInput(session, "grid_size", value = data$input$grid_size)
    updateCheckboxInput(session, "compute_plot", value = data$input$compute_plot)
    updateCheckboxInput(session, "compute_cost_landscape", value = data$input$compute_cost_landscape)
    
    delay(500, {
      updateSelectInput(session, "fn_name", selected = data$input$fn_name)
    })
    
    delay(1000, {
      
      # Set values of dynamic UI
      
      args <- get_args_names_ui()
      
      lapply(args, function(arg_name) {
        session$sendInputMessage(arg_name, list(value = data$input[[arg_name]]))
      })
    })
    
    delay(1500, {
      
      # Set data generation enabled / disabled accordingly
      
      if (data$input$compute_plot) {
        disable("compute_plot")
      } else {
        enable("compute_plot")
      }
      
      if (data$input$compute_cost_landscape) {
        disable("compute_cost_landscape")
      } else {
        enable("compute_cost_landscape")
      }
      
      # Update plot_data object
      
      plot_data$design <<- data$plot_data$design
      plot_data$less <<- data$plot_data$less
      plot_data$domination_counts <<- data$plot_data$domination_counts
      
      if (!is.null(plot_data$design)) {
        show("tabset_plots")
        show("plot_options")
      } else {
        hide("tabset_plots")
        hide("plot_options")
      }
      
      if (!is.null(plot_data$less)) {
        showTab("tabset_plots", "tab_heatmap")
        
        if (smoof::getNumberOfParameters(get_fn()) == 2) {
          showTab("tabset_plots", "tab_plot")
        } else {
          hideTab("tabset_plots", "tab_plot")
        }
      } else {
        hideTab("tabset_plots", "tab_plot")
        hideTab("tabset_plots", "tab_heatmap")
      }
      
      if (!is.null(plot_data$domination_counts)) {
        showTab("tabset_plots", "tab_cost_landscape")
      } else {
        hideTab("tabset_plots", "tab_cost_landscape")
      }
      
    })
  })
}

shinyApp(ui, server)
