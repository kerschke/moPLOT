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
        
        # Bi-Objective BBOB
        
        conditionalPanel(
          condition = "input.function_family == 'biobj_bbob'",
          splitLayout(
            numericInput("biobj_bbob_fid", "Function ID", 1L, min = 1L, max = 55L, step = 1L),
            numericInput("biobj_bbob_iid", "Instance ID", 1L, min = 1L),
            numericInput("biobj_bbob_dim", "Dimensions", 2L, min = 2L, max = 3L, step = 1L)
          )
        ),
        
        # DTLZ
        
        conditionalPanel(
          condition = "input.function_family == 'dtlz'",
          selectInput("fn_name_dtlz", "Function", names(dtlz_functions)),
          splitLayout(
            numericInput("dtlz_dim", "Dimensions", 2L, min = 2L, max = 3L, step = 1L),
            numericInput("dtlz_obj", "Objectives", 2L, min = 2L, max = 3L, step = 1L)
          ),
          
          conditionalPanel(
            condition = "input.fn_name_dtlz == 'DTLZ4'",
            numericInput("dtlz_alpha", "Alpha", 100L)
          )
        ),
        
        # MMF
        
        conditionalPanel(
          condition = "input.function_family == 'mmf'",
          selectInput("fn_name_mmf", "Function", names(mmf_functions)),
          
          conditionalPanel(
            condition = "input.fn_name_mmf == 'MMF14' ||
                         input.fn_name_mmf == 'MMF14a' ||
                         input.fn_name_mmf == 'MMF15' ||
                         input.fn_name_mmf == 'MMF15a'",
            splitLayout(
              numericInput("mmf_dim", "Dimensions", 2L, min = 2L, max = 3L, step = 1L),
              numericInput("mmf_obj", "Objectives", 2L, min = 2L, max = 3L, step = 1L)
            )
          ),
          
          conditionalPanel(
            condition = "input.fn_name_mmf == 'MMF1e'",
            numericInput("mmf_a", "a", exp(1L))
          ),
          
          conditionalPanel(
            condition = "input.fn_name_mmf == 'MMF1z'",
            numericInput("mmf_k", "k", 3L)
          ),
          
          conditionalPanel(
            condition = "input.fn_name_mmf == 'MMF9' ||
                         input.fn_name_mmf == 'MMF11' ||
                         input.fn_name_mmf == 'MMF12' ||
                         input.fn_name_mmf == 'MMF13' ||
                         input.fn_name_mmf == 'MMF14' ||
                         input.fn_name_mmf == 'MMF14a' ||
                         input.fn_name_mmf == 'MMF15' ||
                         input.fn_name_mmf == 'MMF15a'",
            numericInput("mmf_np", "np", 2L, min = 1L, step = 1L)
          ),
          
          conditionalPanel(
            condition = "input.fn_name_mmf == 'MMF12'",
            numericInput("mmf_q", "q", 4L, min = 1L, step = 1L)
          ),
          
          conditionalPanel(
            condition = "input.fn_name_mmf == 'SYMPART-rotated'",
            numericInput("sympart_w", "w", pi / 4, step = pi / 8),
          ),
          
          conditionalPanel(
            condition = "input.fn_name_mmf == 'SYMPART-simple' ||
                         input.fn_name_mmf == 'SYMPART-rotated'",
            splitLayout(
              numericInput("sympart_a", "a", 1L),
              numericInput("sympart_b", "b", 10L),
              numericInput("sympart_c", "c", 8L)
            )
          )
        ),
        
        # ZDT
        
        conditionalPanel(
          condition = "input.function_family == 'zdt'",
          selectInput("fn_name_zdt", "Function", names(zdt_functions)),
          numericInput("zdt_dim", "Dimensions", 2L, min = 2L, max = 3L, step = 1L)
        ),
        
        # Other
        
        conditionalPanel(
          condition = "input.function_family == 'other'",
          selectInput("fn_name_other", "Function", names(other_functions)),
        ),
      ),
      
      wellPanel(
        numericInput("grid.size", "Resolution per dimension", 100, min=20, max=3000, step=1),
        actionButton("evaluate.grid", "Evaluate Grid"),
        # downloadButton('download.grid', "Download Grid"),
        id = "evaluate_grid_panel"
      ),
      
      wellPanel(
        selectInput("plot.type", "Type of plot", c("Select a function first" = "")),
        selectInput("space", "Select space to plot", c("Decision Space"="decision.space", "Objective Space"="objective.space", "Decision + Objective Space"="both")),
        actionButton("update.plot", "Update Plot"),
        id = "update_plot_panel"
      )
      
      # tabsetPanel(
      #   tabPanel("Upload Grid",
      #            fileInput('upload.grid', "Upload Grid", accept = c(".Rds"))
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
  fn <- NULL
  
  hide("evaluate_grid_panel")
  hide("update_plot_panel")

  observe({
    generator_fn <- NULL
    
    hide("evaluate_grid_panel")
    hide("update_plot_panel")

    shiny::validate(
      shiny::need(try({
        if (input$function_family == "biobj_bbob") {
          # Bi-objective BBOB
          
          generator_fn <- smoof::makeBiObjBBOBFunction
          
          args <- list(
            dimensions = input$biobj_bbob_dim,
            fid = input$biobj_bbob_fid,
            iid = input$biobj_bbob_iid
          )
        } else if (input$function_family == "dtlz") {
          # DTLZ
          
          generator_fn <- dtlz_functions[[input$fn_name_dtlz]]
          
          args <- list(
            dimensions = input$dtlz_dim,
            n.objectives = input$dtlz_obj
          )
          
          if (input$fn_name_dtlz == "DTLZ4") {
            args$alpha <- input$dtlz_alpha
          }
        } else if (input$function_family == "zdt") {
          # ZDT
          
          generator_fn <- zdt_functions[[input$fn_name_zdt]]
          
          args <- list(
            dimensions = input$zdt_dim
          )
        } else if (input$function_family == "mmf") {
          # MMF
          
          generator_fn <- mmf_functions[[input$fn_name_mmf]]
          
          args <- list()
          
          if (input$fn_name_mmf == "MMF1e") {
            args$a <- input$mmf_a
          }
          
          if (input$fn_name_mmf == "MMF1z") {
            args$k <- input$mmf_k
          }
          
          if (input$fn_name_mmf %in% c("MMF14", "MMF14a", "MMF15", "MMF15a")) {
            args$dimensions <- input$mmf_dim
            args$n.objectives <- input$mmf_obj
          }
          
          if (input$fn_name_mmf %in% c("MMF9", "MMF11", "MMF12", "MMF13", "MMF14", "MMF14a", "MMF15", "MMF15a")) {
            args$np <- input$mmf_np
          }
          
          if (input$fn_name_mmf == "MMF12") {
            args$q <- input$mmf_q
          }
          
          if (input$fn_name_mmf == "SYMPART-rotated") {
            args$w <- input$sympart_w
          }
          
          if (input$fn_name_mmf %in% c("SYMPART-simple", "SYMPART-rotated")) {
            args$a <- input$sympart_a
            args$b <- input$sympart_b
            args$c <- input$sympart_c
          }
        } else if (input$function_family == "other") {
          generator_fn <- other_functions[[input$fn_name_other]]
          
          args <- list()
        }
        
        if (!is.null(generator_fn)) {
          fn <<- do.call(generator_fn, args)
          
          if (smoof::getNumberOfParameters(fn) == 2) {
            updateSliderInput("grid.size", session = session, value = 100, min=50, max=3000, step=50)
            updateSelectInput(session = session, inputId = "plot.type", choices = list("PLOT" = "PLOT", "Heatmap" = "heatmap", "Cost Landscape (TODO)" = "cost_landscape"))
          } else {
            updateSliderInput("grid.size", session = session, value = 50, min=20, max=200, step=10)
            updateSelectInput(session = session, inputId = "plot.type", choices = list("Onion Layers" = "layers", "MRI Scan" = "scan", "Nondominated" = "pareto"))
          }
          
          show("evaluate_grid_panel")
        }
      }), "Please define a valid function.")
    )
  })
  
  evaluate.grid = function() {
    if (is.null(fn)) {
      return()
    }
    
    # reset stored plot data
    
    plot_data <<- list()
    
    # generate design and evaluate objective space
    
    design <- generateDesign(fn, points.per.dimension = input$grid.size)
    design$obj.space <- calculateObjectiveValues(design$dec.space, fn, parallelize = T)
    
    plot_data$design <<- design
  }
  
  get.plot = function() {
    if (is.null(fn)) {
      return(NULL)
    }
    
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
      switch (input$plot.type,
              heatmap = plotly2DHeatmap(grid, fn, mode = input$space),
              # cost_landscape = plotly3DLayers(grid, fn, mode = input$space),
              PLOT = plotly::ggplotly(ggplotPLOT(grid$dec.space, grid$obj.space, less$sinks, less$height)),
              NULL # if plot.type is invalid
      )
    } else if (d == 3) {
      switch (input$plot.type,
              pareto = plotly3DPareto(grid, fn, mode = input$space),
              layers = plotly3DLayers(grid, fn, mode = input$space),
              scan = plotly3DScan(grid, fn, mode = input$space),
              NULL # if plot.type is invalid
      )
    } else {
      NULL
    }
  }
  
  observeEvent(input$evaluate.grid, {
    disable('evaluate.grid')
    disable('update.plot')
    
    evaluate.grid()
    
    show("update_plot_panel")
    
    enable('update.plot')
    enable('evaluate.grid')
  })
  
  output$plot = plotly::renderPlotly(
    eventReactive(input$update.plot, {
      disable('evaluate.grid')

      options(warn = -1)
      p = get.plot()
      
      enable('evaluate.grid')

      p
    }, event.quoted = T)()
  )
  
  # output$download.grid = downloadHandler(
  #   filename = function() {
  #     fn = test.functions[[as.numeric(input$fn.id)]]
  #     paste0(smoof::getName(fn), "-grid-", input$grid.size, '.Rds')
  #   },
  #   
  #   content = function(file) {
  #     saveRDS(grid, file)
  #   }
  # )
  
  # observeEvent(input$upload.grid, {
  #   if (!endsWith(input$upload.grid$datapath, ".Rds")) {
  #     print(input$upload.grid)
  #     return()
  #   }
  #   
  #   disable('evaluate.grid')
  #   disable('update.plot')
  #   
  #   grid <<- readRDS(input$upload.grid$datapath)
  #   
  #   show('plot.type')
  #   show('space')
  #   show('update.plot')
  #   
  #   enable('evaluate.grid')
  #   enable('update.plot')
  #   enable('download.grid')
  # })
}

shinyApp(ui, server)
