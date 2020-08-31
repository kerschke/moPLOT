library(shiny)
library(shinyjs)

source("dashboardFunctions.R")

# Maximum Upload Size = 50MB
# Might need to be changed for larger files
options("shiny.maxRequestSize" = 50L * (1024L**2L))

ui <- fluidPage(
  useShinyjs(),
  
  h1(markdown("`moPLOT` Landscape Explorer")),
  
  fluidRow(
    
    column(4L,
      h3("Choose MOP"),
      wellPanel(
        selectInput("function_family", "Function Family", c("Select a function family"="", function_families)),
        selectInput("fn_name", "Function", c("Select a function family first"="")),
        
        uiOutput("fn_args")
      ),
      
      wellPanel(
        numericInput("grid_size", "Resolution per dimension", 100, min=20, max=3000, step=1),
        actionButton("evaluate_grid", "Evaluate Grid"),
        # downloadButton('download_grid', "Download Grid"),
        id = "evaluate_grid_panel"
      ),
      
      wellPanel(
        selectInput("plot_type", "Type of plot", c("Select a function first" = "")),
        selectInput("space", "Select space to plot", c("Decision Space"="decision.space", "Objective Space"="objective.space", "Decision + Objective Space"="both")),
        actionButton("update_plot", "Update Plot"),
        id = "update_plot_panel"
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
      plotly::plotlyOutput(outputId = "plot", height = "100%")
    )
  )
)

server <- function(input, output, session) {
  plot_data <- list()
  
  # Hide some parts per default
  
  hide("evaluate_grid_panel")
  hide("update_plot_panel")
  hide("fn_name")
  
  # Reactive getters
  
  get_fn_family <- reactive({
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
    if (any(sapply(args, is.null))) return(NULL)
    
    names(args) <- names(get_default_args())

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
    hide("update_plot_panel")
    
    if (smoof::getNumberOfParameters(fn) == 2) {
      updateSliderInput("grid_size", session = session, value = 200, min=50, max=3000, step=50)
      updateSelectInput(session = session, inputId = "plot_type", choices = list("PLOT" = "PLOT", "Heatmap" = "heatmap", "Cost Landscape (TODO)" = "cost_landscape"))
    } else {
      updateSliderInput("grid_size", session = session, value = 50, min=20, max=200, step=10)
      updateSelectInput(session = session, inputId = "plot_type", choices = list("Onion Layers" = "layers", "MRI Scan" = "scan", "Nondominated" = "pareto"))
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
    print(input$fn_name)
    
    args <- get_default_args()
    
    ui_inputs <- lapply(names(args), function(argument_name) {
      input_args <- list(
        inputId = paste0("args_", argument_name), 
        value = if (is.symbol(args[[argument_name]])) 0 else eval(args[[argument_name]]),
        label = argument_name
      )
      
      # General special cases
      
      if (input_args$label == "dimensions") {
        input_args$label = "Dimensions"
        input_args$value = 2
        input_args$min = 2
        input_args$max = 3
      } else if (input_args$label == "n.objectives") {
        input_args$label = "Objectives"
        input_args$value = 2
        input_args$min = 2
        input_args$max = 3
      }
      
      # Special cases for Bi-Obj BBOB
      
      if (input$function_family == "biobj_bbob") {
        if (input_args$label == "fid") {
          input_args$label = "Function ID"
          input_args$value = 1
          input_args$min = 1
          input_args$max = 55
        }
        
        if (input_args$label == "iid") {
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
  
  # Plotting-related functions
  
  evaluate_grid = function() {
    fn <- get_fn()
    
    req(fn)
    
    print(fn)
    
    # reset stored plot data
    
    plot_data <<- list()
    
    # generate design and evaluate objective space
    
    design <- generateDesign(fn, points.per.dimension = input$grid_size)
    design$obj.space <- calculateObjectiveValues(design$dec.space, fn, parallelize = T)
    
    plot_data$design <<- design
  }
  
  get.plot = function() {
    fn <- get_fn()
    
    req(fn)
    
    if (is.null(plot_data$less)) {
      design <- plot_data$design
      
      gradients <- computeGradientFieldGrid(design)
      
      divergence <- computeDivergenceGrid(gradients$multi.objective, design$dims, design$step.sizes)
      
      # Calculate locally efficient points
      plot_data$less <<- localEfficientSetSkeleton(design, gradients, divergence, integration="fast")
    }
    
    grid <- plot_data$design
    grid$height <- plot_data$less$height
    
    less <- plot_data$less

    d = smoof::getNumberOfParameters(fn)

    if (d == 2) {
      switch (input$plot_type,
              heatmap = plotly2DHeatmap(grid, fn, mode = input$space),
              # cost_landscape = plotly3DLayers(grid, fn, mode = input$space),
              PLOT = plotly::ggplotly(ggplotPLOT(grid$dec.space, grid$obj.space, less$sinks, less$height)),
              NULL # if plot_type is invalid
      )
    } else if (d == 3) {
      switch (input$plot_type,
              pareto = plotly3DPareto(grid, fn, mode = input$space),
              layers = plotly3DLayers(grid, fn, mode = input$space),
              scan = plotly3DScan(grid, fn, mode = input$space),
              NULL # if plot_type is invalid
      )
    } else {
      NULL
    }
  }
  
  observeEvent(input$evaluate_grid, {
    disable('evaluate_grid')
    disable('update_plot')
    
    evaluate_grid()
    
    show("update_plot_panel")
    
    enable('update_plot')
    enable('evaluate_grid')
  })
  
  output$plot = plotly::renderPlotly(
    eventReactive(input$update_plot, {
      disable('evaluate_grid')

      options(warn = -1)
      p = get.plot()
      
      enable('evaluate_grid')

      p
    }, event.quoted = T)()
  )
  
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
