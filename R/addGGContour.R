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
#' @param fn3 [[function()]]\cr
#'   The third objective (if existing) used for computing the multi-objective gradient.
#' @param log.scale [`[logical](1L)`]\cr
#'   Should the resulting heights be displayed on a log-scale? The default is `TRUE`.
#' @param col1 [`[character](1L)`]\cr
#'   Color used for the contour lines of the first objective (default: `"goldenrod1"`).
#' @param col2 [`[character](1L)`]\cr
#'   Color used for the contour lines of the second objective (default: `"white"`).
#' @param col3 [`[character](1L)`]\cr
#'   Color used for the contour lines of the third objective (default: `"cyan3"`).
#' @param n.points [`[integer](1L)`]\cr
#'   Number of points used for computing the contour lines. The default is `30L`.
#' @param ... [any]\cr
#'   Further arguments to be passed to the `geom_tile` function of `ggplot`.
#' @return [`ggplot`]\cr
#'   A `ggplot` object displaying the multi-objective gradient landscape.
#' @examples
#' 
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
