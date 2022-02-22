#' Multi-Objective MinDist Functions
#'
#' Creates a multi-objective MinDist problem. Centers a defined for each
#' function with the `centers.fi` arguments. Box constraints are created at with
#' a distance of 1 to the extremal values of the centers.
#'
#' Decision space dimensionality is determined by the length of each center
#' value in the `centers.fi` lists. If `is.null(centers.f3)`, the problem is
#' bi-objective, otherwise it is tri-objective. For creating bi-objective
#' problems, you can also use `makeBiObjMinDistFunction()` directly.
#'
#' @param centers.f1 [[`list`]]\cr List of centers for the first objective.
#' @param centers.f2 [[`list`]]\cr List of centers for the second objective.
#' @param centers.f3 [[`list`]] (optional)\cr List of centers for the third
#'   objective, if function has three objectives.
#'
#' @return A multi-objective smoof function
#' @export
#'
#' @examples
#' 
makeMinDistFunction = function(centers.f1 = list(c(-2, -1), c(2, 1)),
                               centers.f2 = list(c(-2, 1), c(2, -1)),
                               centers.f3 = list(c(0, 1), c(0, -1))) {
  if (is.null(centers.f3)) {
    f <- function(x) {
      y1 <- min(sapply(centers.f1, function(center) sqrt(sum((x - center) ** 2))))
      y2 <- min(sapply(centers.f2, function(center) sqrt(sum((x - center) ** 2))))
      
      c(y1, y2)
    }
    
    all.centers <- Reduce(rbind, c(centers.f1, centers.f2))
  } else {
    f <- function(x) {
      y1 <- min(sapply(centers.f1, function(center) sqrt(sum((x - center) ** 2))))
      y2 <- min(sapply(centers.f2, function(center) sqrt(sum((x - center) ** 2))))
      y3 <- min(sapply(centers.f3, function(center) sqrt(sum((x - center) ** 2))))
      
      c(y1, y2, y3)
    }
    
    all.centers <- Reduce(rbind, c(centers.f1, centers.f2, centers.f3))
  }
  
  lower <- apply(all.centers, 2, min) - 1
  upper <- apply(all.centers, 2, max) + 1
  
  smoof::makeMultiObjectiveFunction(name = "MinDist Function", id = "mindist", description = "", fn = f,
                                    par.set = ParamHelpers::makeNumericParamSet(len = ncol(all.centers), lower = lower, upper = upper))
}

#' @rdname makeMinDistFunction
#' @export
makeBiObjMinDistFunction = function(centers.f1 = list(c(-2, -1), c(2, 1)),
                                    centers.f2 = list(c(-2, 1), c(2, -1))) {
  makeMinDistFunction(centers.f1, centers.f2, centers.f3 = NULL)
}
