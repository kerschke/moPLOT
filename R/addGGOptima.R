#' Add (Known) Local Optima of Smoof Functions to Gradient Field Heatmap.
#'
#' @description
#'   If the underlying single-objective functions are `smoof`-functions,
#'   extract their local optima and add them to the gradient field heatmap.
#'
#' @template arg_ggplot
#' @template arg_fni
#' @param fn3 [[function()]]\cr
#'   The third objective (if existing) used for computing the multi-objective gradient.
#' @param symbol1 [`[integer](1L) | [character](1L)`]\cr
#'   Symbol used for indicating the local optima of the first objective
#'   (default: `21L`, i.e., a filled circle).
#' @param symbol2 [`[integer](1L) | [character](1L)`]\cr
#'   Symbol used for indicating the local optima of the second objective
#'   (default: `22L`, i.e., a filled square).
#' @param symbol3 [`[integer](1L) | [character](1L)`]\cr
#'   Symbol used for indicating the local optima of the third objective
#'   (default: `24L`, i.e., a filled triangle).
#' @param ... [any]\cr
#'   Further arguments to be passed to the `geom_tile` function of `ggplot`.
#' @return [`ggplot`]\cr
#'   A `ggplot` object displaying the multi-objective gradient landscape.
#' @examples
#' 
#' @export
addGGOptima = function(g, fn1, fn2, fn3,
  symbol1 = 21L, symbol2 = 22L, symbol3 = 24L, ...) {

  assertClass(x = g, classes = "ggplot")
  df.optima = NULL
  if (!missing(fn1)) {
    assertFunction(fn1)
    assertClass(fn1, "smoof_function")
    x.opt = smoof::getLocalOptimum(fn1)$params
    if (is.null(x.opt)) {
      x.opt = smoof::getGlobalOptimum(fn1)$param
    }
    if (!is.null(x.opt)) {
      df.optima = rbind(df.optima, cbind(x.opt, obj = "fn1", shp = symbol1))
    }
  }
  if (!missing(fn2)) {
    assertFunction(fn2)
    assertClass(fn2, "smoof_function")
    x.opt = smoof::getLocalOptimum(fn2)$params
    if (is.null(x.opt)) {
      x.opt = smoof::getGlobalOptimum(fn2)$param
    }
    if (!is.null(x.opt)) {
      df.optima = rbind(df.optima, cbind(x.opt, obj = "fn2", shp = symbol2))
    }
  }
  if (!missing(fn3)) {
    assertFunction(fn3)
    assertClass(fn3, "smoof_function")
    x.opt = smoof::getLocalOptimum(fn3)$params
    if (is.null(x.opt)) {
      x.opt = smoof::getGlobalOptimum(fn3)$param
    }
    if (!is.null(x.opt)) {
      df.optima = rbind(df.optima, cbind(x.opt, obj = "fn3", shp = symbol3))
    }
  }
  if (!is.null(df.optima)) {
    g = g +
      geom_point(data = df.optima, mapping = aes_string(x = "x1", y = "x2", shape = "shp"), ...) +
      scale_shape_identity()
  }
  return(g)
}
