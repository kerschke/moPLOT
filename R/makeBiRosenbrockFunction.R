#' Bi-Objective Rosenbrock Function
#' 
#' Creates one particular parametrization of a bi-objective Rosenbrock function.
#'
#' @return Multi-objective smoof function
#' @export
#'
#' @examples
#' 
makeBiRosenbrockFunction = function() {
  f1 <- function(x) {
    (1 - x[1]) ** 2 + 1 * (x[2] - x[1] ** 2) ** 2
  }
  
  f2 <- function(x) {
    (1 + x[1]) ** 2 + 1 * (-(x[2] - 3) - x[1] ** 2) ** 2
  }
  
  f <- function(x) {
    c(f1(x), f2(x))
  }
  
  smoof::makeMultiObjectiveFunction(
    name = "Bi-Rosenbrock Function", id = "bi_rosenbrock_function", description = "", fn = f,
    par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-1.5, 0), upper = c(1.5, 3)))
}

