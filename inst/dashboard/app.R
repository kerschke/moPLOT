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
      h3("Select MOP"),
      
      wellPanel(
        selectInput("function_family", "Function family", c("Select a function family"="", function_families)),
        selectInput("fn_name", "Function", c("Select a function family first"="")),
        
        uiOutput("fn_args")
      ),
      
      div(
        h3("Generate Data"),
        wellPanel(
          numericInput("grid_size", "Resolution per dimension", 100, min=20, max=3000, step=1),
          checkboxInput("compute_plot", "Compute PLOT and heatmap", TRUE),
          checkboxInput("compute_cost_landscape", "Compute cost landscape", FALSE),
          actionButton("evaluate_grid", "Evaluate"), # icon = icon('th')
          # downloadButton('download_grid', "Download Grid"),
        ),
        id = "evaluate_grid_panel"
      )
      
      # tabsetPanel(
      #   tabPanel("Upload Grid",
      #            fileInput('upload_grid', "Upload Grid", accept = c(".Rds"))
      #   ),
      #   id = "grid.tabs",
      #   type = "pills"
      # ),
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
          selectInput("space", "Select space to plot", c("Decision space" = "decision.space", "Objective space" = "objective.space", "Decision + objective space" = "both")),
          div(
            selectInput("three_d_approach", "3D Approach", c("MRI Scan" = "scan", "Onion Layers" = "layers", "Nondominated" = "pareto")),
            conditionalPanel("input.three_d_approach == 'scan'",
              selectInput("scan_direction", "Scan Direction", c("x₁" = "x1", "x₂" = "x2", "x₃" = "x3"), selected = "x3")
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
  plot_data <- list()
  
  outputOptions(output, suspendWhenHidden = FALSE)
  
  reset_plots <- function() {
    plot_data <<- list()
    
    hideTab("tabset_plots", "tab_plot")
    hideTab("tabset_plots", "tab_heatmap")
    hideTab("tabset_plots", "tab_cost_landscape")
    
    hide("tabset_plots")
    hide("plot_options")
    
    enable("compute_plot")
    enable("compute_cost_landscape")
  }
  
  # Hide some parts per default
  
  hide("evaluate_grid_panel")
  hide("fn_name")
  
  reset_plots()
  
  # Reactive getters
  
  get_fn_family <- reactive({
    # TODO MinDist
    switch(
      input$function_family,
      biobj_bbob = biobj_bbob_functions,
      dtlz = dtlz_functions,
      mmf = mmf_functions,
      mop = mop_functions,
      zdt = zdt_functions,
      other = other_functions
    )
  })
  
  get_generator_fn <- reactive({
    get_fn_family()[[input$fn_name]]
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
      input[[argument_name]]
    })
    
    # Only return args if correctly connected to output
    if (all(sapply(args, is.null))) return(list())
    
    names(args) <- names(get_default_args())
    args <- args[!sapply(args, function(arg) (is.null(arg) || arg == ""))]

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
      updateSelectInput(session = session, inputId = "plot_type", choices = list("PLOT" = "PLOT", "Heatmap" = "heatmap", "Cost Landscape" = "cost_landscape"))
      hide("three_d_only")
    } else {
      updateSliderInput("grid_size", session = session, value = 50, min=20, max=200, step=10)
      show("three_d_only")
    }
    
    show("evaluate_grid_panel")
  })
  
  observe({
    fn_family <- get_fn_family()
    
    updateSelectInput(session, "fn_name", choices = names(fn_family))

    if (is.null(fn_family) || input$function_family == "biobj_bbob") {
      hide("fn_name")
    } else {
      show("fn_name")
    }
  })
  
  # Dynamically created view with function parameters
  
  output$fn_args = renderUI({
    args <- get_default_args()
    
    ui_inputs <- lapply(names(args), function(argument_name) {
      input_args <- list(
        inputId = paste0("args_", argument_name), 
        value = if (is.symbol(args[[argument_name]])) NULL else tryCatch(eval(args[[argument_name]]), error = function(e) NULL),
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
      
      if (input$function_family == "biobj_bbob") {
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
      
      design <- generateDesign(fn, points.per.dimension = input$grid_size)
      design$obj.space <- calculateObjectiveValues(design$dec.space, fn, parallelize = T)
      
      plot_data$design <<- design
    }
    
    if (input$compute_plot && is.null(plot_data$less)) {
      design <- plot_data$design
      
      gradients <- computeGradientFieldGrid(design, prec.angle = 1)
      
      divergence <- computeDivergenceGrid(gradients$multi.objective, design$dims, design$step.sizes)
      
      # Calculate locally efficient points
      plot_data$less <<- localEfficientSetSkeleton(design, gradients, divergence, integration="fast")
      
      showTab("tabset_plots", "tab_plot")
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
  
  get_plot = function(plot_data, plot_type, space, three_d_approach) {
    fn <- get_fn()
    
    req(fn)

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
      get_plot(plot_data, "PLOT", input$space, input$three_d_approach)
    }, quoted = TRUE)()
  )
  
  output$heatmap = plotly::renderPlotly({
    reactive({
      get_plot(plot_data, "heatmap", input$space, input$three_d_approach)
    }, quoted = TRUE)()
  })
  
  output$cost_landscape = plotly::renderPlotly({
    reactive({
      get_plot(plot_data, "cost_landscape", input$space, input$three_d_approach)
    }, quoted = TRUE)()
  })
  
  # output$download_grid = downloadHandler(
  #   filename = function() {
  #     fn = test.functions[[as.numeric(input$fn.id)]]
  #     paste0(smoof::getName(fn), "-grid-", input$grid_size, '.Rds')
  #   },
  #   
  #   content = function(file) {
  #     saveRDS(grid, file)
  #   }
  # )
  
  # observeEvent(input$upload_grid, {
  #   if (!endsWith(input$upload_grid$datapath, ".Rds")) {
  #     print(input$upload_grid)
  #     return()
  #   }
  #   
  #   disable('evaluate_grid')
  #   disable('update_plot')
  #   
  #   grid <<- readRDS(input$upload_grid$datapath)
  #   
  #   show('plot_type')
  #   show('space')
  #   show('update_plot')
  #   
  #   enable('evaluate_grid')
  #   enable('update_plot')
  #   enable('download_grid')
  # })
}

shinyApp(ui, server)
