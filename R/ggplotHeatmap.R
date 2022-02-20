#' Visualize Heatmap of Multi-Objective Gradients.
#'
#' @description
#'   Visualization of the multi-objective gradient landscape by means of
#'   a heatmap and on the basis of the `ggplot2`-package.
#'
#' @param df [[data.frame()]]\cr
#'   Data frame as returned by [computeCumulatedPathLengths()].
#' @param var1 [`[character](1L)`]\cr
#'   Name of the variable indicating the first dimension (default: `"x1"`).
#' @param var2 [`[character](1L)`]\cr
#'   Name of the variable indicating the second dimension (default: `"x2"`).
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
ggplotHeatmap = function(df, var1 = "x1", var2 = "x2", log.scale = TRUE, impute.zero = TRUE,
  minimalistic.image = FALSE, color.palette, legend.position, ...) {

  assertDataFrame(df, any.missing = FALSE, min.cols = 3L)
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
    ## if no information on the color palette is provided, this
    ## function tries to sequentially tries to use tim.colors,
    ## viridis or at last the terrain.colors
    inst.pkgs = rownames(installed.packages())
    if ("fields" %in% inst.pkgs) {
      color.palette = fields::tim.colors(500L)
    } else if ("viridisLite" %in% inst.pkgs) {
      color.palette = viridisLite::viridis(500L,
        alpha = 1, begin = 0, end = 1, direction = 1, option = "D")
    } else {
      color.palette = terrain.colors(500L)
    }
  }
  
  ## create the heatmap (colored tiles)
  g = ggplot() +
    geom_raster(data = df, mapping = aes_string(x = var1, y = var2, fill = "height"), ...) +
    theme_minimal()

  if (log.scale) {
    ## if results are shown on log-scale, provide a 'pretty' height scale
    pretty_breaks_log = function(){
      function(x) {
        axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE)
      }
    }
    g = g + scale_fill_gradientn(colors = color.palette, trans = "log", breaks = pretty_breaks_log())
  } else {
    g = g + scale_fill_gradientn(colors = color.palette)
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
