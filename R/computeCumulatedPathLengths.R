#' Compute the Cumulated Path Lengths.
#'
#' @description
#'   Given a grid of points in the search space, along with their corresponding
#'   multi-objective gradients, this function will compute (for each point of the grid)
#'   the length of the cumulated path from a point towards its attracting local
#'   efficient point.
#'
#' @note 
#'   ATTENTION: Only turn off the sanity checks (\code{check.data = FALSE}),
#'   if you can ensure that all input parameters are provided in the correct format.
#'
#' @template arg_centers
#' @template arg_gradients
#' @template arg_precnorm
#' @template arg_precveclen
#' @template arg_checkdata
#' @return [\code{\link{data.frame}}]\cr
#'   Returns a \code{data.frame}, which appends the cumulated path lengths to the points
#'   provided by \code{centers}.
#' @examples
#' # Define two single-objective test problems:
#' fn1 = function(x) sum((x - c(0.2, 1))^2)
#' fn2 = function(x) sum((x - c(0.5, 0.5))^2)
#'
#' # Define a grid of points:
#' points = as.matrix(expand.grid(x1 = seq(0, 0.7, 0.005), x2 = seq(0, 1.25, 0.005)))
#' 
#' # Compute the corresponding gradients:
#' gradients = computeGradientField(points, fn1, fn2)
#'
#' # Now, compute the cumulated path lengths:
#' x = computeCumulatedPathLengths(points, gradients)
#'
#' # Finally, we can visualize the resulting multi-objective "landscape":
#' ggplotHeatmap(x, hide.legend = TRUE)
#' @export
computeCumulatedPathLengths = function(centers, gradients, local.efficient.ids = numeric(0), prec.vector.length = 1e-3, prec.norm = 1e-6, cumulate.gradient.length = TRUE, fix.diagonals = FALSE, check.data = TRUE) {

  if (check.data) {
    if (is.data.frame(centers)) {
      centers = as.matrix(centers)
    }
    if (is.data.frame(gradients)) {
      gradients = as.matrix(gradients)
    }
    # FIXME: so far, only supporting 2D-problems
    assertMatrix(x = centers, mode = "numeric", min.cols = 2L, null.ok = FALSE)
    assertMatrix(x = gradients, mode = "numeric", min.cols = 2L, null.ok = FALSE)
    assertTRUE(all(dim(centers) == dim(gradients)))
    # index = order(centers[,2], centers[,1])
    index = do.call(order, lapply(rev(seq_col(centers)), function(i) centers[,i]))
    if (any(index != seq_row(centers))) {
      centers = centers[index,,drop = FALSE]
      gradients = gradients[index,,drop = FALSE]
    }
    # assertTRUE(nrow(centers) == length(unique(centers[,1])) * length(unique(centers[,2])))
    assertTRUE(nrow(centers) == prod(apply(as.matrix(centers), 2L, function(x) length(unique(x)))))
    assertNumber(prec.vector.length, lower = 0, finite = TRUE, null.ok = FALSE)
    assertNumber(prec.norm, lower = 0, finite = TRUE, null.ok = FALSE)
  }
  integrated = cumulateGradientsCPP(centers, gradients, local.efficient.ids, prec.vector.length, prec.norm, fix.diagonals, cumulate.gradient.length)
  dim(integrated$height) = c(length(integrated$height), 1)
  colnames(integrated$height) = c("height")
  
  return(integrated)
}
