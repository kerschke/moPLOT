library(shiny)
library(shinyjs)

source("dashboardFunctions.R")

# Maximum Upload Size = 50MB
# Might need to be changed for larger files
options("shiny.maxRequestSize" = 50L * (1024L**2L))

ui <- fluidPage(
  useShinyjs(),
  
  h1(markdown("`moPLOT` Dashboard <Working Title>")),
  
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
      ),
      
      wellPanel(
        numericInput("grid.size", "Resolution per dimension", 100, min=20, max=3000, step=1),
        actionButton("update.grid", "Update Plots"),
        # downloadButton('download.grid', "Download Grid"),
        id = "update_plots"
      )
      
      # tabsetPanel(
      #   tabPanel("Calculate Grid",
      #            numericInput("grid.size", "Resolution per dimension", 100, min=20, max=3000, step=1),
      #            actionButton("update.grid", "Calculate Grid"),
      #            downloadButton('download.grid', "Download Grid"),
      #   ),
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
  grid <- NULL
  fn <- NULL
  hide("update_plots")

  observe({
    generator_fn <- NULL
    
    hide("update_plots")

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
        }
        
        if (!is.null(generator_fn)) {
          fn <<- do.call(generator_fn, args)
          show("update_plots")
        }
      }), "Please define a valid function.")
    )
    
    print(fn)
  })
  
  update.grid = function() {
    if (is.null(fn)) {
      return()
    }
    
    grid <<- moPLOT::generateDesign(fn, points.per.dimension = input$grid.size)
    grid$obj.space <<- moPLOT::calculateObjectiveValues(grid$dec.space, fn, parallelize = T)
    
    gradients <- moPLOT::computeGradientFieldGrid(grid)
    
    divergence <- computeDivergenceGrid(gradients$multi.objective, grid$dims, grid$step.sizes)
    
    # Calculate locally efficient points
    less <<- localEfficientSetSkeleton(grid, gradients, divergence, integration="fast")
    grid$height <<- less$height
  }
  
  get.plot = function() {
    if (is.null(fn)) {
      return(NULL)
    }
    
    d = smoof::getNumberOfParameters(fn)
    n = smoof::getNumberOfObjectives(fn)
    
    if (d == 2) {
      moPLOT::plotly2DHeatmap(grid, fn, mode="decision.space")
      # ggplotPLOT(grid$dec.space, grid$obj.space, less$sinks, less$height) %>% ggplotly()
    } else {
      # switch (input$plot.type,
      #         pareto = moPLOT::plotly3DPareto(grid, fn, mode=input$space),
      #         layers = moPLOT::plotly3DLayers(grid, fn, mode=input$space),
      #         scan = moPLOT::plotly3DScan(grid, fn, mode=input$space),
      #         NULL # if plot.type is invalid
      # )
    }
  }
  
  output$plot = plotly::renderPlotly(
    eventReactive(input$update.grid, {
      disable('update.grid')

      update.grid()
      
      options(warn = -1)
      p = get.plot()
      
      enable('update.grid')

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
  
  observeEvent(input$upload.grid, {
    if (!endsWith(input$upload.grid$datapath, ".Rds")) {
      print(input$upload.grid)
      return()
    }
    
    disable('update.grid')
    disable('update.plot')
    
    grid <<- readRDS(input$upload.grid$datapath)
    
    show('plot.type')
    show('space')
    show('update.plot')
    
    enable('update.grid')
    enable('update.plot')
    enable('download.grid')
  })
}

shinyApp(ui, server)
