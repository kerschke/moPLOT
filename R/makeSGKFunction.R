#' Multi-Objective SGK Function
#'
#' @param dimensions `c(2, 3)`\cr
#' Number of (decision space) dimensions.
#'
#' @return Multi-objective smoof function
#' @export
#'
#' @examples
#' 
makeSGKFunction = function(dimensions = 2) {
  g = function(x, h, center) {
    h / (1 + 4 * (sum((x - center) ** 2)))
  }
  
  if (dimensions == 2) {
    unimodal_center = c(2 / 3, 1)
    unimodal_h = 1
    
    centers = list(
      c(0.5, 0),
      c(0.25, 2 / 3),
      c(1, 1)
    )
    hs = c(1.5, 2, 3)
  } else if (dimensions == 3) {
    unimodal_center = c(2 / 3, 1, 0)
    unimodal_h = 1
    
    centers = list(
      c(0.5, 0, 0.75),
      c(0.25, 2 / 3, 0.5),
      c(1, 1, 0)
    )
    hs = c(1.5, 2, 3)
  }
  
  lower = rep(-0.25, dimensions)
  upper = rep(1.25, dimensions)
  
  f = function(x) c(
    1 - g(x, unimodal_h, unimodal_center),
    1 - max(
      sapply(seq_along(centers), function(i) g(x, hs[i], centers[[i]]))
    )
  )
  
  smoof::makeMultiObjectiveFunction(
    name = "SGK Function", id = "sgk_function", description = "", fn = f,
    par.set = ParamHelpers::makeNumericParamSet(len = dimensions, lower = lower, upper = upper))
}
