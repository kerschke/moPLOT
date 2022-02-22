#' Add Gradient Vector Arrows to Gradient Field Heatmap.
#'
#' @description
#'   Add vectors indicating the gradient field to the heatmap of the
#'   cumulated path lengths.
#'
#' @template arg_ggplot
#' @template arg_centers
#' @template arg_gradients
#' @param fac [`[numeric]`]\cr
#'   Factor used for scaling the lengths of the gradient vectors. The default is `0.025`.
#' @param arrow_len [`[unit]`]\cr
#'   Length of the arrow tips (default: `unit(0.075 / 2, "inches")`).
#' @param nColumns [`[integer]`]\cr
#'   How many columns of arrows should be drawn? The default is `10`.
#' @param nRows [`[integer]`]\cr
#'   How many rows of arrows should be drawn? The default is `10`.
#' @param ... [any]\cr
#'   Further arguments to be passed to the `geom_tile` function of `ggplot`.
#' @return [`ggplot`]\cr
#'   A `ggplot` object displaying the multi-objective gradient landscape.
#' @examples
#' 
#' @export
addGGArrows = function(g, centers, gradients, fac = 0.025,
  arrow_len = unit(0.075 / 2, "inches"), nColumns = 10L, nRows = 10L, ...) {

  assertClass(x = g, classes = "ggplot")
  assertIntegerish(x = nColumns, len = 1L, any.missing = FALSE, null.ok = FALSE)
  assertIntegerish(x = nRows, len = 1L, any.missing = FALSE, null.ok = FALSE)
  assertMatrix(x = centers, mode = "numeric", min.cols = 2L, null.ok = FALSE, any.missing = FALSE)
  assertMatrix(x = gradients, mode = "numeric", min.cols = 2L, null.ok = FALSE, any.missing = FALSE)
  assertTRUE(all(dim(centers) == dim(gradients)))

  ## ensure that there are no more arrows per dimension than cells
  x1 = sort(unique(centers[,1]))
  nc = length(x1)
  x2 = sort(unique(centers[,2]))
  nr = length(x2)
  nColumns = min(nColumns, nc)
  nRows = min(nRows, nr)

  ## center arrows
  if (nColumns < nc) {
    i1 = seq(1, nc, nc / nColumns)
    i1 = i1 + floor((nc - max(i1)) / 2)
    mx1 = x1[i1]
  }
  if (nRows < nr) {
    i2 = seq(1, nr, nr / nRows)
    i2 = i2 + floor((nr - max(i2)) / 2)
    mx2 = x2[i2]
  }

  ## reduce the data accordingly
  index = which((centers[,1] %in% mx1) & (centers[,2] %in% mx2))
  centers = centers[index, , drop = FALSE]
  gradients = gradients[index, , drop = FALSE]

  ## produce a data frame, which contains the arrow paths
  arrow.df = t(vapply(seq_row(centers), function(i) {
    st = centers[i,]
    ziel = st + fac * gradients[i,]
    return(c(st, ziel))
  }, double(4L)))
  arrow.df = as.data.frame(rbind(arrow.df[, c(1, 2)], arrow.df[, c(3, 4)]))
  colnames(arrow.df) = c("x1", "x2")
  arrow.df$arrow = rep(seq_len(nrow(arrow.df) / 2), 2)

  ## add the arrow paths to the existing ggplot
  g = g +
    geom_path(
      data = arrow.df,
      mapping = aes(x = x1, y = x2, group = arrow),
      arrow = arrow(length = arrow_len), ...
    )

  return(g)
}


## add gradient vector arrows to a base R plot
addArrowsToPlot = function(points, gradients, fac = 0.01 / 2, length = 0.01, ...) {
  for (i in seq_row(points)) {
    st = points[i,]
    ziel = st + fac * gradients[i,]
    arrows(st[1], st[2], ziel[1], ziel[2], col = "grey", length = length, ...)
  }
}
