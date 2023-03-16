#' Visualize Multi-Objective Gradient Landscape in the Objective Space.
#'
#' @description
#'   Use the coloring of the cumulated path lengths (see
#'   [computeCumulatedPathLengths()]) and visualize the
#'   corresponding points in the objective space.
#'
#' @param df [[data.frame()]]\cr
#'   Data frame containing the positions of the individuals (per row) and their
#'   corresponding cumulated path length (denoted: `"height"`) as returned by
#'   [computeCumulatedPathLengths()].
#' @param var1 [`[character](1L)`]\cr
#'   Name of the variable indicating the first objective (default: `"y1"`).
#' @param var2 [`[character](1L)`]\cr
#'   Name of the variable indicating the second objective (default: `"y2"`).
#' @param log.scale [`[logical](1L)`]\cr
#'   Should the resulting heights be displayed on a log-scale? The default is `TRUE`.
#' @param impute.zero [`[logical](1L)`]\cr
#'   Should height values, which are exactly zero be imputed by a value half the magnitude
#'   of the smallest non-zero height? Otherwise ggplot will automatically color the
#'   corresponding tiles by a color representing `NA` values (usually grey).
#'   Note that this parameter is only relevant in case of `log.scale = TRUE`.
#'   The default is `TRUE`.
#' @param minimalistic.image [`[logical](1L)`]\cr
#'   Should all information surrounding the image (axes, legends, background, etc.) be discarded?
#'   The default is `FALSE`.
#' @param color.palette [[character()]]\cr
#'   Vector of colors used for visualizing the different heights of the landscape. By default,
#'   this function tries to use the color palettes from `fields::tim.color` or
#'   `viridis`. However, if neither of these packages is installed, it will use
#'   `terrain.colors`.
#' @param legend.position [`[character](1L)`]\cr
#'   On which side of the plot, should the legend be located? If this information is not provided
#'   and `minimalisitic.image = FALSE`, the legend will be placed on the right side.
#' @param ... [any]\cr
#'   Further arguments to be passed to the `geom_tile` function of `ggplot`.
#' @return [`ggplot`]\cr
#'   A `ggplot` object displaying the multi-objective gradient landscape.
#' @examples
#' 
#' @export
ggplotObjectiveSpace = function(df, var1 = "y1", var2 = "y2", log.scale = TRUE,
  impute.zero = TRUE, minimalistic.image = FALSE,
  color.palette, legend.position, ...) {

  assertDataFrame(df, any.missing = FALSE, min.cols = 3L)
  assertCharacter(x = var1, any.missing = FALSE, len = 1L, null.ok = FALSE)
  assertCharacter(x = var2, any.missing = FALSE, len = 1L, null.ok = FALSE)
  assertSubset(c(var1, var2, "height"), colnames(df))
  assertLogical(log.scale, any.missing = FALSE, len = 1L, null.ok = FALSE)
  assertLogical(impute.zero, any.missing = FALSE, len = 1L, null.ok = FALSE)
  assertLogical(minimalistic.image, any.missing = FALSE, len = 1L, null.ok = FALSE)

  if (log.scale & impute.zero) {
    ## impute heights of zero for log-scale visualizations
    z = df$height
    mz = min(z[z != 0])
    z[z == 0] = mz / 2
    df$height = z
  }

  if (missing(color.palette)) {
    color.palette = viridisLite::turbo(500L, alpha = 1, begin = 0.05, end = 0.95, direction = 1)
  }

  ## order df such that the efficient sets are located on top of the inferior points
  df = df[order(df$height, decreasing = TRUE), , drop = FALSE]
  
  ## create the surface (colored points) in the objective space
  g = ggplot() +
    geom_point(
      data = df,
      mapping = aes_string(x = var1, y = var2, color = "height"), ...) +
    theme_minimal()

  if (log.scale) {
    ## if results are shown on log-scale, provide a 'pretty' height scale
    pretty_breaks_log = function(){
      function(x) {
        axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE)
      }
    }
    g = g + scale_color_gradientn(colors = color.palette, trans = "log", breaks = pretty_breaks_log())
  } else {
    g = g + scale_color_gradientn(colors = color.palette)
  }

  ## Modify axes labels (place indices as subscripts)
  g = g +
    xlab(variable.as.expression(var1)) +
    ylab(variable.as.expression(var2))
  
  if (minimalistic.image) {
    ## in case of minimalistic images, remove the color legend, axis labels and ticks
    g = g + minimalistic.theme()
  } else if (!missing(legend.position)) {
    ## position legend
    assertCharacter(legend.position, len = 1L, null.ok = FALSE)
    assertChoice(legend.position, c("none", "left", "right", "bottom", "top"), null.ok = FALSE)
    if (legend.position != "none") {
      g = g + labs(fill = "Height")
    }
    g = g + theme(legend.position = legend.position)
  }

  return(g)
}
