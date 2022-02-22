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

#' Multi-Objective MPM2 Functions
#'
#' The `makeBiObjMPM2Function()` and `makeTriObjMPM2Function()` functions allow
#' the creation of multi-objective functions based on the MPM2 generator
#' implemented in [smoof::makeMPM2Function()].
#'
#' @param dimensions [[`numeric`]]\cr
#' Number of (decision space) dimensions of the MPM2 problems.
#' @param n.peaks.1 [[`numeric`]]\cr
#' Number of peaks for the first constituent MPM2 problem.
#' @param topology.1 `"random" || "funnel"`\cr
#' Topology of the first constituent MPM2 problem.
#' @param seed.1 [[`numeric`]]\cr
#' Random seed of the first constituent MPM2 problem.
#' @param n.peaks.2 [[`numeric`]]\cr
#' Number of peaks for the second constituent MPM2 problem.
#' @param topology.2 `"random" || "funnel"`\cr
#' Topology of the second constituent MPM2 problem.
#' @param seed.2 [[`numeric`]]\cr
#' Random seed of the second constituent MPM2 problem.
#' @param n.peaks.3 [[`numeric`]]\cr
#' Number of peaks for the third constituent MPM2 problem.
#' @param topology.3 `"random" || "funnel"`\cr
#' Topology of the third constituent MPM2 problem.
#' @param seed.3 [[`numeric`]]\cr
#' Random seed of the third constituent MPM2 problem.
#'
#' @return Multi-objective smoof function
#' @export
#'
#' @examples
#' 
makeBiObjMPM2Function = function(dimensions = 2, n.peaks.1 = 3, topology.1 = "random", seed.1 = 4,
                                 n.peaks.2 = 3, topology.2 = "random", seed.2 = 8) {
  f1 <- smoof::makeMPM2Function(n.peaks.1, dimensions, topology.1, seed.1)
  f2 <- smoof::makeMPM2Function(n.peaks.2, dimensions, topology.2, seed.2)
  
  smoof::makeMultiObjectiveFunction(
    name = paste0("Bi-MPM2 (", smoof::getName(f1), ", ", smoof::getName(f2), ")"),
    id = paste0("bi_mpm2_", smoof::getID(f1), "_", smoof::getID(f2)),
    fn = function(x) {
      drop(cbind(f1(x), f2(x)))
    },
    par.set = ParamHelpers::makeNumericParamSet("x", len = dimensions, lower = 0, upper = 1),
    vectorized = TRUE
  )
  
  # smoof::makeGOMOPFunction(dimensions = dimensions, funs = list(f1, f2))
}

#' @rdname makeBiObjMPM2Function
#' @export
makeTriObjMPM2Function = function(dimensions = 2, n.peaks.1 = 3, topology.1 = "random", seed.1 = 4,
                                  n.peaks.2 = 3, topology.2 = "random", seed.2 = 8,
                                  n.peaks.3 = 3, topology.3 = "random", seed.3 = 12) {
  f1 <- smoof::makeMPM2Function(n.peaks.1, dimensions, topology.1, seed.1)
  f2 <- smoof::makeMPM2Function(n.peaks.2, dimensions, topology.2, seed.2)
  f3 <- smoof::makeMPM2Function(n.peaks.3, dimensions, topology.3, seed.3)
  
  smoof::makeMultiObjectiveFunction(
    name = paste0("Tri-MPM2 (", smoof::getName(f1), ", ", smoof::getName(f2), ", ", smoof::getName(f3), ")"),
    id = paste0("tri_mpm2_", smoof::getID(f1), "_", smoof::getID(f2), "_", smoof::getID(f3)),
    fn = function(x) {
      drop(cbind(f1(x), f2(x), f3(x)))
    },
    par.set = ParamHelpers::makeNumericParamSet("x", len = dimensions, lower = 0, upper = 1),
    vectorized = TRUE
  )
  
  # smoof::makeGOMOPFunction(dimensions = dimensions, funs = list(f1, f2, f3))
}

# Aspar Function ====

f1_1 = function(x) (x[1]**4 - 2*x[1]**2 + x[2]**2 + 1)
f2_1 = function(x) ((x[1] + 0.5)**2 + (x[2]-2)**2)
f3_1 = function(x) ((x[1] + 0.25) ** 4 + 3 * (x[2] - 1) ** 2)
f_2d2d = function(x) c(f1_1(x), f2_1(x))
f_2d3d = function(x) c(f1_1(x), f2_1(x), f3_1(x))

f1_2 = function(x) (x[1]**4 - 2*x[1]**2 + x[2]**2 + 1 + x[3] ** 2)
f2_2 = function(x) ((x[1] + 0.5)**2 + (x[2]-2)**2 + (x[3] - 1) ** 4)
f3_2 = function(x) ((x[1] + 0.25) ** 4 + 3 * (x[2] - 1) ** 2 + (x[3] - 1) ** 2)
f_3d2d = function(x) c(f1_2(x), f2_2(x))
f_3d3d = function(x) c(f1_2(x), f2_2(x), f3_2(x))

#' Multi-Objective Aspar Function
#'
#' @param dimensions [[`numeric`]]\cr
#'   Number of dimensions
#' @param n.objectives [[`numeric`]]\cr
#'   Number of objectives
#'
#' @return Multi-objective smoof function
#' @export
#'
#' @examples
#' 
makeAsparFunction <- function(dimensions = 2, n.objectives = 2) {
  if (dimensions == 2 && n.objectives == 2) {
    smoof::makeMultiObjectiveFunction(name = "Aspar Function: 2D->2D", id = "aspar_2d2d", description = "", fn = f_2d2d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-2,-1), upper = c(2,3)))
  } else if (dimensions == 2 && n.objectives == 3) {
    smoof::makeMultiObjectiveFunction(name = "Aspar Function: 2D->3D", id = "aspar_2d3d", description = "", fn = f_2d3d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-2,-1), upper = c(2,3)))
  } else if (dimensions == 3 && n.objectives == 2) {
    smoof::makeMultiObjectiveFunction(name = "Aspar Function: 3D->2D", id = "aspar_3d2d", description = "", fn = f_3d2d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 3, lower = c(-2,-1,-2), upper = c(2,3,2)))
  } else if (dimensions == 3 && n.objectives == 3) {
    smoof::makeMultiObjectiveFunction(name = "Aspar Function: 3D->3D", id = "aspar_3d3d", description = "", fn = f_3d3d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 3, lower = c(-2,-1,-2), upper = c(2,3,2)))
  }
}

# SGK Function ====

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

# Bi-Rosenbrock Function ====

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

