library(shiny)
library(shinyjs)

source("benchmarkFunctions.R")

# Maximum Upload Size = 50MB
# Might need to be changed for larger files
options("shiny.maxRequestSize" = 50 * (1024**2))

ui <- fluidPage(
  useShinyjs(),
  
  h1(markdown("`moPLOT` Dashboard <Working Title>")),
  
  fluidRow(
    
    column(4,
      wellPanel(
        selectInput("fn.id", "Multi-objective function", c("Select a function"="",test.function.ids)),
        tabsetPanel(
          tabPanel("Calculate Grid",
                   numericInput("grid.size", "Resolution per dimension", 100, min=20, max=3000, step=1),
                   actionButton("update.grid", "Calculate Grid"),
                   downloadButton('download.grid', "Download Grid"),
          ),
          tabPanel("Upload Grid",
                   fileInput('upload.grid', "Upload Grid", accept = c(".Rds"))
          ),
          id = "grid.tabs",
          type = "pills"
        ),
      ),
      
      conditionalPanel(condition = "input['update.grid']",
        wellPanel(
          selectInput("plot.type", "Type of plot", c("Select a function first" = "")),
          selectInput("space", "Select space to plot", c("Decision Space"="decision.space", "Objective Space"="objective.space", "Decision + Objective Space"="both")),
          actionButton("update.plot", "Update Plot", class="btn-primary")
        )
      )
    ),
    
    column(8,
      plotly::plotlyOutput(outputId = "plot", height = "100%")
    )
  )
)

server <- function(input, output, session) {
  # initially, hide all these elements:
  hide('grid.size')
  hide('grid.tabs')
  hide('update.grid')
  hide('download.grid')
  hide('upload.grid')
  hide('plot.type')
  hide('space')
  hide('update.plot')
  
  grid = NULL
  
  observe({
    fn = test.functions[[as.numeric(input$fn.id)]]
    
    if (is.null(fn)) {
      return()
    }
    
    hide('plot.type')
    hide('update.plot')
    hide('space')
    
    disable('download.grid')
    
    grid <<- NULL
    
    show('grid.tabs')
    show('grid.size')
    show('update.grid')
    show('download.grid')
    show('upload.grid')
    
    if (smoof::getNumberOfParameters(fn) == 2) {
      updateSliderInput("grid.size", session = session, value = 100, min=20, max=3000, step=1)
      updateSelectInput(session = session, inputId = "plot.type", choices = list("Heatmap" = "heatmap"))
    } else {
      updateSliderInput("grid.size", session = session, value = 50, min=20, max=200, step=1)
      updateSelectInput(session = session, inputId = "plot.type", choices = list("Onion Layers" = "layers", "MRI Scan" = "scan", "Nondominated" = "pareto"))
    }
  })
  
  update.grid = function() {
    fn = test.functions[[as.numeric(input$fn.id)]]
    if (is.null(fn)) {
      return()
    }
    
    grid <<- moPLOT::generateDesign(fn, points.per.dimension = input$grid.size)
    grid$obj.space <<- moPLOT::calculateObjectiveValues(grid$dec.space, fn, parallelize = T)
    
    gradients <<- moPLOT::computeGradientFieldGrid(grid)
    grid$height <<- moPLOT::computeCumulatedPathLengths(grid$dec.space, gradients$multi.objective, fix.diagonals = T)$height
  }
  
  observeEvent(input$update.grid, {
    disable('update.grid')
    disable('update.plot')
    
    update.grid()
    
    show('plot.type')
    show('space')
    show('update.plot')
    
    enable('update.grid')
    enable('update.plot')
    enable('download.grid')
  })
  
  get.plot = function() {
    fn = test.functions[[as.numeric(input$fn.id)]]
    if (is.null(fn)) {
      return(NULL)
    }
    
    d = smoof::getNumberOfParameters(fn)
    n = smoof::getNumberOfObjectives(fn)
    
    if (d == 2) {
      moPLOT::plotly2DHeatmap(grid, fn, mode=input$space)
    } else {
      switch (input$plot.type,
              pareto = moPLOT::plotly3DPareto(grid, fn, mode=input$space),
              layers = moPLOT::plotly3DLayers(grid, fn, mode=input$space),
              scan = moPLOT::plotly3DScan(grid, fn, mode=input$space),
              NULL # if plot.type is invalid
      )
    }
  }
  
  output$plot = plotly::renderPlotly(
    eventReactive(input$update.plot, {
      disable('update.grid')
      disable('update.plot')
      
      options(warn = -1)
      p = get.plot()
      
      enable('update.grid')
      enable('update.plot')
      
      p
    }, event.quoted = T)()
  )
  
  output$download.grid = downloadHandler(
    filename = function() {
      fn = test.functions[[as.numeric(input$fn.id)]]
      paste0(smoof::getName(fn), "-grid-", input$grid.size, '.Rds')
    },
    
    content = function(file) {
      saveRDS(grid, file)
    }
  )
  
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
