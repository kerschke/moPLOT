#' Run the moPLOT dashboard
#'
#' `runDashboard` starts a [`shiny`](https://CRAN.R-project.org/package=shiny)
#' application, which allows the interactive exploration of the optimization landscape
#' of many integrated benchmark suites.
#'
#'@export
runDashboard = function() {
  # the app requires the smoof package
  if (!requireNamespace("smoof", quietly = TRUE)) {
    stop("The package `smoof` is required for the dashboard. Please install it.", call. = FALSE)
  }
  
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The package `shiny` is required for the dashboard. Please install it.", call. = FALSE)
  }
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("The package `reticulate` is required for the dashboard. Please install it.", call. = FALSE)
  }
  
  appDir = system.file("dashboard", package = "moPLOT")
  if (appDir == "") {
    stop("Could not find dashboard directory. Try re-installing `moPLOT`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
