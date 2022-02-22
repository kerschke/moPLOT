library(shiny)
library(shinyjs)
library(fields)
library(smoof)
library(moPLOT)

# ==== Setup Function Sets ====

benchmark_sets <- c(
  "(Extended) Bi-objective BBOB" = "biobj_bbob",
  "DTLZ Functions" = "dtlz",
  "MinDist Functions" = "mindist",
  "MMF Functions" = "mmf",
  "MOP Functions" = "mop",
  "MPM2 Generator" = "mpm2",
  "ZDT Functions" = "zdt",
  "Other" = "other"
)

# (Extended) Bi-objective BBOB ====

biobj_bbob_functions <- list(
  # "Bi-objective BBOB" = smoof::makeBiObjBBOBFunction
  "(Extended) Bi-objective BBOB" = makeExtendedBiObjBBOBFunction
)

# MinDist Functions ====

mindist_functions = list(
  "Bi-objective MinDist" = makeBiObjMinDistFunction,
  "Tri-objective MinDist" = makeMinDistFunction
)

# MPM2 Functions ====

mpm2_functions = list(
  "Bi-objective MPM2" = makeBiObjMPM2Function,
  "Tri-objective MPM2" = makeTriObjMPM2Function
)

# Benchmark Sets Extracted From smoof ====

dtlz_functions <- list(
  "DTLZ1" = smoof::makeDTLZ1Function,
  "DTLZ2" = smoof::makeDTLZ2Function,
  "DTLZ3" = smoof::makeDTLZ3Function,
  "DTLZ4" = smoof::makeDTLZ4Function,
  "DTLZ5" = smoof::makeDTLZ5Function,
  "DTLZ6" = smoof::makeDTLZ6Function,
  "DTLZ7" = smoof::makeDTLZ7Function
)

mmf_functions <- list(
  "MMF1" = smoof::makeMMF1Function,
  "MMF1e" = smoof::makeMMF1eFunction,
  "MMF1z" = smoof::makeMMF1zFunction,
  "MMF2" = smoof::makeMMF2Function,
  "MMF3" = smoof::makeMMF3Function,
  "MMF4" = smoof::makeMMF4Function,
  "MMF5" = smoof::makeMMF5Function,
  "MMF6" = smoof::makeMMF6Function,
  "MMF7" = smoof::makeMMF7Function,
  "MMF8" = smoof::makeMMF8Function,
  "MMF9" = smoof::makeMMF9Function,
  "MMF10" = smoof::makeMMF10Function,
  "MMF11" = smoof::makeMMF11Function,
  "MMF12" = smoof::makeMMF12Function,
  "MMF13" = smoof::makeMMF13Function,
  "MMF14" = smoof::makeMMF14Function,
  "MMF14a" = smoof::makeMMF14aFunction,
  "MMF15" = smoof::makeMMF15Function,
  "MMF15a" = smoof::makeMMF15aFunction,
  "SYMPART-simple" = smoof::makeSYMPARTsimpleFunction,
  "SYMPART-rotated" = smoof::makeSYMPARTrotatedFunction,
  "Omni test" = smoof::makeOmniTestFunction
)

zdt_functions = list(
  "ZDT1" = smoof::makeZDT1Function,
  "ZDT2" = smoof::makeZDT2Function,
  "ZDT3" = smoof::makeZDT3Function,
  "ZDT4" = smoof::makeZDT4Function,
  "ZDT6" = smoof::makeZDT6Function
)

mop_functions = list(
  "MOP1" = smoof::makeMOP1Function,
  "MOP2" = smoof::makeMOP2Function,
  "MOP3" = smoof::makeMOP3Function,
  "MOP4" = smoof::makeMOP4Function,
  "MOP5" = smoof::makeMOP5Function,
  "MOP6" = smoof::makeMOP6Function,
  "MOP7" = smoof::makeMOP7Function
)

# Summary of "other" functions ====

other_functions = list(
  "Aspar" = makeAsparFunction,
  "BiRosenbrock" = makeBiRosenbrockFunction,
  "BiSphere" = smoof::makeBiSphereFunction, # does not work as expected and is kinda boring
  "BK1" = smoof::makeBK1Function,
  "Dent" = smoof::makeDentFunction,
  "ED1" = smoof::makeED1Function,
  "ED2" = smoof::makeED2Function,
  "Kursawe" = smoof::makeKursaweFunction,
  "SGK" = makeSGKFunction,
  "Viennet" = smoof::makeViennetFunction
)

# Setup reticulate for MPM2 function interfaces

if (!reticulate::virtualenv_exists("moPLOT")) {
  reticulate::virtualenv_create("moPLOT")
}

reticulate::use_virtualenv("moPLOT")

if (!reticulate::py_numpy_available()) {
  reticulate::py_install("numpy")
}

# Maximum Upload Size = 50MB
# Might need to be changed for larger files
options("shiny.maxRequestSize" = 50L * (1024L**2L))

ui <- fluidPage(
  useShinyjs(),
  
  h1(code("moPLOT"), "Landscape Explorer"),
  
  fluidRow(
    
    column(
      4L,
      tabsetPanel(
       type = "pills",
       id = "fn_select",
       
       tabPanel(
         "Select MOP",
         value = "smoof_mop",
         p(),
         wellPanel(
           selectInput(
             "benchmark_set",
             "Benchmark set",
             c("Select a benchmark set"="", benchmark_sets)
           ),
           selectInput(
             "fn_name",
             "Function",
             c("Select a benchmark set first"="")
           ),
           uiOutput("fn_args"),
           div(
             id = "evaluate_design_panel",
             numericInput(
               "grid_size",
               "Resolution per dimension",
               100,
               min = 20,
               max = 3000,
               step = 1
             ),
             splitLayout(
               actionButton("evaluate_design", "Evaluate", style = "width: 100%; bottom: 0"),
               downloadButton("download_data", "Download", style = "width: 100%")
             )
           )
         ),
       ),
       
       tabPanel(
         "Upload Data",
         value = "upload_data",
         p(),
         wellPanel(
           fileInput("upload_data", "Upload Data", accept = c(".csv")),
           helpText("Please provide a CSV with a column for each decision space dimension (x1, x2 [and x3])",
                    "and the corresponding objective values (y1, y2 [and y3]), with decision space variables",
                    "covering a grid of evenly spaced values per dimension.")
         )
       )
      ),
      
      div(
       h3("Visualization Options"),
       wellPanel(
         checkboxInput("compute_plot", "Enable PLOT and heatmap", TRUE),
         checkboxInput("compute_local_dominance", "Enable local dominance", TRUE),
         checkboxInput("compute_cost_landscape", "Enable cost landscape", FALSE),
       ),
       id = "viz_options_panel"
      )
    ),
    
    column(
      8,
      tabsetPanel(
       id = "tabset_plots",
       selected = "tab_plot",
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
         "Set Transitions",
         plotOutput("set_transitions", height = "500px"),
         value = "tab_set_transitions"
       ),
       tabPanel(
         "Contours",
         plotly::plotlyOutput("contours", height = "500px"),
         value = "tab_contours"
       ),
       tabPanel(
         "Local Dominance",
         plotly::plotlyOutput("local_dominance", height = "500px"),
         value = "tab_local_dominance"
       ),
       tabPanel(
         "Cost Landscape",
         plotly::plotlyOutput("cost_landscape", height = "500px"),
         value = "tab_cost_landscape"
       )
      ),
      div(
       h3("Plot Options"),
       wellPanel(
         selectInput(
           "space",
           "Space to plot",
           c(
             "Decision Space" = "decision.space",
             "Objective Space" = "objective.space",
             "Decision + Objective Space" = "both"
           )
         ),
         div(
           selectInput(
             "three_d_approach",
             "3D approach",
             c(
               "MRI Scan" = "scan",
               "Onion Layers" = "layers",
               "Nondominated" = "pareto"
             )
           ),
           conditionalPanel("input.three_d_approach == 'scan'",
                            selectInput(
                              "scan_direction",
                              "Scan direction",
                              c("x₁" = "x1", "x₂" = "x2", "x₃" = "x3"),
                              selected = "x3"
                            )
           ),
           id = "three_d_only"
         ),
         selectInput("show_nondominated",
                     "Contour mode",
                     c("Show nondominated points" = "TRUE",
                       "Contours only" = "FALSE"),
                     selected = "TRUE")
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

    enable("compute_plot")
    enable("compute_cost_landscape")
    enable(selector = "button")
    
    plot_data$design <- NULL
    plot_data$less <- NULL
    plot_data$domination_counts <- NULL
    plot_data$dummy_fn_upload <- NULL
    plot_data$ld_height <- NULL
  }
  
  # Hide some parts per default
  
  hide("evaluate_design_panel")
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
    if (input$fn_select == "upload_data") {
      # We just need to return a dummy function that has the
      # correct lower, upper bounds and a suitable name
      plot_data$dummy_fn_upload
    } else {
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
    }
  })
  
  # Compute PLOT Data ====
  
  compute_plot_data <- reactive({
    design <- plot_data$design
    
    gradients <- computeGradientFieldGrid(design, prec.norm = 0)
    
    divergence <- computeDivergenceGrid(gradients$multi.objective, design$dims, design$step.sizes)
    
    # Calculate locally efficient points
    localEfficientSetSkeleton(design, gradients, divergence, integration="fast", with.basins = TRUE)
  })
  
  # Compute Local Dominance Data ====
  
  compute_ld_height <- reactive({
    design <- plot_data$design
    
    ld_data <- moPLOT:::computeLocalDominance(design$obj.space, design$dims)
    
    basins <- sapply(ld_data$basins, function(v) {
      if (length(v) == 1) v
      else NA
    })
    
    chob <- moPLOT:::changeOfBasin(basins, design$dims, ld_data$locally_efficient_ids)
    basins[setdiff(chob$ridges, ld_data$locally_efficient_ids)] <- NA
    
    display_height <- rep(0.5, length(ld_data$basins))
    display_height[ld_data$locally_efficient_ids] <- 0
    display_height[is.na(basins)] <- 1
    
    display_height <- matrix(display_height, ncol = 1)
    colnames(display_height) <- c("height")
    
    display_height
  })
  
  # Compute Dominance Counts ====
  
  compute_dominance_counts <- reactive({
    nds <- ecr::doNondominatedSorting(t(plot_data$design$obj.space))
    dom_counts <- cbind(height = nds$dom.counter + 1)
    dom_counts
  })
  
  # Observers that change the UI dynamically ====
  
  observe({
    fn <- get_fn()
    
    req(fn)
    
    hide("evaluate_design_panel")
    
    if (smoof::getNumberOfParameters(fn) == 2) {
      updateSliderInput("grid_size", session = session, value = 500L, min = 50L, max = 600L, step = 50L)
      hide("three_d_only")
    } else {
      updateSliderInput("grid_size", session = session, value = 50L, min = 20L, max = 100L, step = 10L)
      show("three_d_only")
    }
    
    show("evaluate_design_panel")
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
  
  # Which tabs to show? ====
  
  observe({
    fn <- get_fn()
    req(fn)
    
    if (smoof::getNumberOfParameters(fn) == 2) {
      showTab("tabset_plots", "tab_contours")
    } else {
      hideTab("tabset_plots", "tab_contours")
      updateTabsetPanel(inputId = "tabset_plots", selected = NULL)
    }
    
    if (input$compute_plot) {
      showTab("tabset_plots", "tab_plot")
      showTab("tabset_plots", "tab_heatmap")
      showTab("tabset_plots", "tab_set_transitions")
    } else {
      hideTab("tabset_plots", "tab_plot")
      hideTab("tabset_plots", "tab_heatmap")
      hideTab("tabset_plots", "tab_set_transitions")
    }
    
    if (input$compute_cost_landscape) {
      showTab("tabset_plots", "tab_cost_landscape")
    } else {
      hideTab("tabset_plots", "tab_cost_landscape")
    }
    
    if (input$compute_local_dominance) {
      showTab("tabset_plots", "tab_local_dominance")
    } else {
      hideTab("tabset_plots", "tab_local_dominance")
    }
  })
  
  # Dynamically created view with function parameters ====
  
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
          input_args$value = 1L
          input_args$min = 1L
          # input_args$max = 55L
          input_args$max = 92L # Extended Bi-Objective BBOB
        }
        
        if (argument_name == "iid") {
          input_args$label = "Instance ID"
          input_args$value = 1L
          input_args$min = 1L
          input_args$max = 15L # Only the verified instances
        }
      }
      
      input_type <- ifelse(is.numeric(input_args$value), numericInput, textInput)
      
      do.call(input_type, input_args)
    })
    
    ui_inputs$cellArgs <- list(style = "width: 49%")
    
    do.call(flowLayout, ui_inputs)
  })
  
  # Generating data ====
  
  observeEvent(c(input$grid_size, get_fn()), {
    if (input$fn_select == "smoof_mop") {
      # reset stored plot data when resolution or function changes
      reset_plots()
    }
  })
  
  observeEvent(input$evaluate_design, {
    fn <- get_fn()
    req(fn)
    
    disable(selector = "button")
    
    if (is.null(plot_data$design)) {
      # generate design and evaluate objective space
      
      design <- moPLOT::generateDesign(fn, points.per.dimension = input$grid_size)
      design$obj.space <- calculateObjectiveValues(design$dec.space, fn, parallelize = FALSE)
      
      plot_data$design <<- design
    }
    
    enable(selector = "button")
    
    show("tabset_plots")
    show("plot_options")
    
    disable("evaluate_design")
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
      },
      local_dominance = {
        req(plot_data$ld_height)
        plot_data$ld_height
      },
      NULL # if plot_type is invalid
    )
    
    less <- plot_data$less
    
    if (plot_type == "PLOT") {
      sinks <- less$sinks
    } else if (plot_type == "local_dominance") {
      ld_data <- plot_data$ld_data
      sinks <- ld_data$locally_efficient_ids
    } else {
      sinks <- NULL
    }
    
    d = smoof::getNumberOfParameters(fn)
    
    if (d == 2) {
      switch (plot_type,
              heatmap = plotly2DHeatmap(
                grid,
                fn,
                mode = space
              ),
              cost_landscape = plotly2DHeatmap(
                grid,
                fn,
                mode = space
              ),
              PLOT = plotly2DPLOT(
                grid$dec.space,
                grid$obj.space,
                less$sinks,
                less$height,
                fn,
                mode = space
              ),
              local_dominance = plotly2DHeatmap(
                grid,
                fn,
                mode = space,
                colorscale = moPLOT:::plotlyColorscale(moPLOT:::gray.colorscale),
                impute.zero = FALSE,
                log.scale = FALSE
              ),
              NULL # if plot_type is invalid
      )
    } else if (d == 3) {
      if (plot_type == "local_dominance") {
        colorscale.sinks = moPLOT:::plotlyColorscale(c("#000000", "#000000"))
        colorscale.heatmap = moPLOT:::plotlyColorscale(moPLOT:::gray.colorscale)
      } else {
        colorscale.sinks = moPLOT:::plotlyColorscale()
        colorscale.heatmap = moPLOT:::plotlyColorscale(moPLOT:::gray.colorscale)
      }
      
      switch (three_d_approach,
              pareto = plotly3DPareto(grid, fn, mode = space),
              layers = plotly3DLayers(grid, fn, sinks, mode = space,
                                      colorscale.sinks = colorscale.sinks,
                                      colorscale.heatmap = colorscale.heatmap),
              scan = plotly3DScan(grid, fn, sinks, mode = space, frame = input$scan_direction,
                                  colorscale.sinks = colorscale.sinks,
                                  colorscale.heatmap = colorscale.heatmap),
              NULL # if plot_type is invalid
      )
    } else {
      NULL
    }
  }
  
  output$plot = plotly::renderPlotly({
    reactive({
      req(plot_data$design)
      
      disable("plot")
      
      if (is.null(plot_data$less) && input$compute_plot) {
        plot_data$less <<- compute_plot_data()
      }
      
      print("Updating PLOT")
      
      p <- get_plot("PLOT", input$space, input$three_d_approach)
      enable("plot")
      return(p)
    }, quoted = TRUE)()
  })
  
  output$heatmap = plotly::renderPlotly({
    reactive({
      req(plot_data$design)
      
      disable("heatmap")
      
      if (is.null(plot_data$less) && input$compute_plot) {
        plot_data$less <<- compute_plot_data()
      }
      
      print("Updating Heatmap")
      
      p <- get_plot("heatmap", input$space, input$three_d_approach)
      enable("heatmap")
      return(p)
    }, quoted = TRUE)()
  })
  
  output$contours = plotly::renderPlotly({
    reactive({
      req(plot_data$design)
      
      disable("contours")
      print("Updating Contours")
      
      p <- plotly2DContours(plot_data$design,
                            show.nondominated = input$show_nondominated == "TRUE")
      enable("contours")
      return(p)
    }, quoted = TRUE)()
  })
  
  output$local_dominance = plotly::renderPlotly({
    reactive({
      req(plot_data$design)
      
      disable("local_dominance")
      
      if (is.null(plot_data$ld_height) && input$compute_local_dominance) {
        plot_data$ld_height <<- compute_ld_height()
      }
      
      print("Updating Local Dominance")
      
      p <- get_plot("local_dominance", input$space, input$three_d_approach)
      enable("local_dominance")
      return(p)
    }, quoted = TRUE)()
  })
  
  output$cost_landscape = plotly::renderPlotly({
    reactive({
      req(plot_data$design)
      
      disable("cost_landscape")
      
      if (is.null(plot_data$domination_counts) && input$compute_cost_landscape) {
        plot_data$domination_counts <<- compute_dominance_counts()
      }
      
      print("Updating Cost Landscape")
      
      p <- get_plot("cost_landscape", input$space, input$three_d_approach)
      enable("cost_landscape")
      return(p)
    }, quoted = TRUE)()
  })
  
  output$set_transitions = renderPlot({
    reactive({
      req(plot_data$design)
      
      if (is.null(plot_data$less) && input$compute_plot) {
        plot_data$less <<- compute_plot_data()
      }
      
      less <- plot_data$less
      basins <- plot_data$less$basins
      design <- plot_data$design
      
      n_descent_target <- table(basins[basins != -1])
      transition_sinks <- intersect(which(less$set_transitions != -1), less$sinks)
      
      set_transitions <- t(sapply(transition_sinks, function(s) {
        c(from = basins[s], to = less$set_transitions[s])
      }))
      set_transitions <- set_transitions[!duplicated(set_transitions),,drop=FALSE]
      
      sets <- lapply(sort(unique(basins[basins != -1])), function(b) {
        set_ids <- intersect(which(basins == b), less$sinks)
        list(
          dec_space = design$dec.space[set_ids,],
          obj_space = design$obj.space[set_ids,]
        )
      })
      
      # === Set Visualization ===
      
      compute_reach_proportions <- function(graph, weights) {
        sapply(igraph::V(graph), function(v) {
          reachable_v <- igraph::subcomponent(graph, v, "in")
          sum(weights[reachable_v]) / sum(weights)
        })
      }
      
      compute_nondominated_sets <- function(sets) {
        obj_space <- Reduce(rbind, lapply(sets, function(set) set$obj_space))
        
        # When do we change between different sets of the trace?
        set_change <- cumsum(sapply(sets, function(set) nrow(set$obj_space)))
        
        nd <- obj_space %>% t %>% ecr::nondominated()
        
        sapply(seq_along(set_change), function(i) {
          if (i == 1) {
            lower = 1
          } else {
            lower = set_change[i - 1] + 1
          }
          
          sum(nd[lower:set_change[i]])
        })
      }
      
      library(tidygraph)
      library(ggraph)
      
      if (length(set_transitions) != 0) {
        tbl_transitions <- as_tbl_graph(set_transitions)
        prop <- compute_reach_proportions(tbl_transitions, n_descent_target)
      } else {
        tbl_transitions <- NULL
        prop <- n_descent_target / sum(n_descent_target)
      }
      
      set_nd_counts <- compute_nondominated_sets(sets)
      node_color <- ifelse(set_nd_counts > 0, "green", "red")
      
      print((1 / prop[set_nd_counts > 0]) %>% sort(decreasing = TRUE))
      
      ggraph(tbl_transitions, layout = "stress") +
        geom_edge_fan(arrow = arrow(length = unit(4, "mm")),
                      end_cap = circle(4, "mm")) +
        geom_node_point(aes(size = prop), color = node_color) +
        scale_radius(limits = c(0,1))
    }, quoted = TRUE)()
  })
  
  observeEvent(input$tabset_plots, {
    if (input$tabset_plots == "tab_contours") {
      hide("space")
      show("show_nondominated")
    } else {
      show("space")
      hide("show_nondominated")
    }
  })
  
  # Up- and Download
  
  output$download_data = downloadHandler(
    filename = function() {
      fn = get_fn()
      req(fn)
      
      paste0(smoof::getName(fn), '.csv')
    },
    
    content = function(file) {
      req(plot_data$design)
      
      df <- cbind(plot_data$design$dec.space, plot_data$design$obj.space)
      
      write.csv(df, file)
    }
  )
  
  observeEvent(input$upload_data, {
    if (!endsWith(input$upload_data$datapath, ".csv")) {
      print(input$upload_data)
      return()
    }
    
    df <- as.matrix(read.csv(input$upload_data$datapath))
    
    reset_plots()
    
    design <- list()
    
    # design$dec.space
    
    if (all(c("x1", "x2", "x3") %in% colnames(df))) {
      design$dec.space <- df[,c("x1", "x2", "x3")]
    } else if (all(c("x1", "x2") %in% colnames(df))) {
      design$dec.space <- df[,c("x1", "x2")]
    } else {
      # ERROR!
      return()
    }
    
    # design$obj.space
    
    if (all(c("y1", "y2", "y3") %in% colnames(df))) {
      design$obj.space <- df[,c("y1", "y2", "y3")]
    } else if (all(c("y1", "y2") %in% colnames(df))) {
      design$obj.space <- df[,c("y1", "y2")]
    } else {
      # ERROR!
      return()
    }
    
    # Miscellaneous
    
    design$dims <- apply(design$dec.space, 2, function(x) length(unique(x)))
    lower <- apply(design$dec.space, 2, function(x) min(x))
    upper <- apply(design$dec.space, 2, function(x) max(x))
    design$step.sizes <- (upper - lower) / (design$dims - 1)
    
    # Define a dummy function with the name and correct bounds
    
    param_set <- ParamHelpers::makeNumericParamSet(
      len = ncol(design$dec.space), lower = lower, upper = upper)
    
    plot_data$dummy_fn_upload <<- smoof::makeMultiObjectiveFunction(
      name = basename(input$upload_data$datapath),
      id = basename(input$upload_data$datapath),
      par.set = param_set,
      fn = function(x) rep(0, ncol(design$obj.space))
    )
    
    plot_data$design <<- design
  })
}

shinyApp(ui, server)
