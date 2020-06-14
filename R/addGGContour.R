#' Add Contour Lines to Gradient Field Heatmap.
#'
#' @description
#'   Add contour lines of the separate objectives to the heatmap of the
#'   cumulated path lengths.
#'
#' @template arg_ggplot
#' @template arg_lower
#' @template arg_upper
#' @template arg_fni
#' @param fn3 [\code{\link{function}}]\cr
#'   The third objective (if existing) used for computing the multi-objective gradient.
#' @param log.scale [\code{\link{logical}(1L)}]\cr
#'   Should the resulting heights be displayed on a log-scale? The default is \code{TRUE}.
#' @param col1 [\code{\link{character}(1L)}]\cr
#'   Color used for the contour lines of the first objective (default: \code{"goldenrod1"}).
#' @param col2 [\code{\link{character}(1L)}]\cr
#'   Color used for the contour lines of the second objective (default: \code{"white"}).
#' @param col3 [\code{\link{character}(1L)}]\cr
#'   Color used for the contour lines of the third objective (default: \code{"cyan3"}).
#' @param n.points [\code{\link{integer}(1L)}]\cr
#'   Number of points used for computing the contour lines. The default is \code{30L}.
#' @param ... [any]\cr
#'   Further arguments to be passed to the \code{geom_tile} function of \code{ggplot}.
#' @return [\code{ggplot}]\cr
#'   A \code{ggplot} object displaying the multi-objective gradient landscape.
#' @examples
#' # Define two single-objective test problems and a grid of points:
#' fn1 = function(x) sum((x - c(0.2, 1))^2)
#' fn2 = function(x) sum((x - c(0.5, 0.5))^2)
#' points = as.matrix(expand.grid(x1 = seq(0, 0.7, 0.005), x2 = seq(0, 1.25, 0.005)))
#' 
#' # Compute the corresponding gradients and the cumulated path lengths:
#' gradients = computeGradientField(points, fn1, fn2)
#' x = computeCumulatedPathLengths(points, gradients)
#'
#' # Visualize the resulting multi-objective "landscape":
#' g = ggplotHeatmap(x)
#' g
#'
#' # Add dashed contour lines to the plot:
#' addGGContour(g = g, lower = c(0, 0), upper = c(0.7, 1.25),
#'   fn1 = fn1, fn2 = fn2, linetype = "dashed")
#' @export
addGGContour = function(g, lower, upper, fn1, fn2, fn3, log.scale = TRUE,
  col1 = "goldenrod1", col2 = "white", col3 = "cyan3", n.points = 30L, ...) {

  assertClass(x = g, classes = "ggplot")
  assertNumeric(x = lower, len = 2L, finite = TRUE, any.missing = FALSE)
  assertNumeric(x = upper, len = 2L, finite = TRUE, any.missing = FALSE)
  assertTRUE(all(lower < upper))
  assertCharacter(x = col1, len = 1L, any.missing = FALSE, null.ok = FALSE)
  assertCharacter(x = col2, len = 1L, any.missing = FALSE, null.ok = FALSE)
  assertCharacter(x = col3, len = 1L, any.missing = FALSE, null.ok = FALSE)

  ## create a help grid, based on which the contour lines are constructed
  x.grid = expand.grid(
    x1 = seq(lower[1L], upper[1L], length.out = n.points),
    x2 = seq(lower[2L], upper[2L], length.out = n.points))
  attr(x.grid, "out.attrs") = NULL

  ## add contour lines for first objective
  if (!missing(fn1)) {
    assertFunction(fn1)
    x.grid$y1 = apply(as.matrix(x.grid[, 1:2]), 1L, fn1)
    if (log.scale) {
      x.grid$y1 = log10(x.grid$y1)
    }
    g = g +
      geom_contour(data = x.grid, mapping = aes_string(x = "x1", y = "x2", z = "y1"), colour = col1, ...)
  }

  ## add contour lines for second objective
  if (!missing(fn2)) {
    assertFunction(fn2)
    x.grid$y2 = apply(as.matrix(x.grid[, 1:2]), 1L, fn2)
    if (log.scale) {
      x.grid$y2 = log10(x.grid$y2)
    }
    g = g +
      geom_contour(data = x.grid, mapping = aes_string(x = "x1", y = "x2", z = "y2"), colour = col2, ...)
  }

  ## add contour lines for third objective
  if (!missing(fn3)) {
    assertFunction(fn3)
    x.grid$y3 = apply(as.matrix(x.grid[, 1:2]), 1L, fn3)
    if (log.scale) {
      x.grid$y3 = log10(x.grid$y3)
    }
    g = g +
      geom_contour(data = x.grid, mapping = aes_string(x = "x1", y = "x2", z = "y3"), colour = col3, ...)
  }

  return(g)
}
