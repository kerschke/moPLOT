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
#' @param rotated.1 [[`logical`]]\cr
#' Should peaks be rotated in the first constituent MPM2 problem?
#' @param peak.shape.1 `"ellipse" || "sphere"`\cr
#' The peak shape (sphere or ellipse) for the first constituent MPM2 problem
#' @param n.peaks.2 [[`numeric`]]\cr
#' Number of peaks for the second constituent MPM2 problem.
#' @param topology.2 `"random" || "funnel"`\cr
#' Topology of the second constituent MPM2 problem.
#' @param seed.2 [[`numeric`]]\cr
#' Random seed of the second constituent MPM2 problem.
#' @param rotated.2 [[`logical`]]\cr
#' Should peaks be rotated in the second constituent MPM2 problem?
#' @param peak.shape.2 `"ellipse" || "sphere"`\cr
#' The peak shape (sphere or ellipse) for the second constituent MPM2 problem
#' @param n.peaks.3 [[`numeric`]]\cr
#' Number of peaks for the third constituent MPM2 problem.
#' @param topology.3 `"random" || "funnel"`\cr
#' Topology of the third constituent MPM2 problem.
#' @param seed.3 [[`numeric`]]\cr
#' Random seed of the third constituent MPM2 problem.
#' @param rotated.3 [[`logical`]]\cr
#' Should peaks be rotated in the third constituent MPM2 problem?
#' @param peak.shape.3 `"ellipse" || "sphere"`\cr
#' The peak shape (sphere or ellipse) for the third constituent MPM2 problem
#'
#' @return Multi-objective smoof function
#' @export
#'
#' @examples
#' 
makeBiObjMPM2Function = function(
    dimensions = 2,
    n.peaks.1 = 3, topology.1 = "random", seed.1 = 4, rotated.1 = TRUE, peak.shape.1 = "ellipse",
    n.peaks.2 = 3, topology.2 = "random", seed.2 = 8, rotated.2 = TRUE, peak.shape.2 = "ellipse"
  ) {
  f1 <- smoof::makeMPM2Function(n.peaks.1, dimensions, topology.1, seed.1, rotated.1, peak.shape.1)
  f2 <- smoof::makeMPM2Function(n.peaks.2, dimensions, topology.2, seed.2, rotated.2, peak.shape.2)
  
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
makeTriObjMPM2Function = function(
    dimensions = 2,
    n.peaks.1 = 3, topology.1 = "random", seed.1 = 4, rotated.1 = TRUE, peak.shape.1 = "ellipse",
    n.peaks.2 = 3, topology.2 = "random", seed.2 = 8, rotated.2 = TRUE, peak.shape.2 = "ellipse",
    n.peaks.3 = 3, topology.3 = "random", seed.3 = 12, rotated.3 = TRUE, peak.shape.3 = "ellipse") {
  f1 <- smoof::makeMPM2Function(n.peaks.1, dimensions, topology.1, seed.1, rotated.1, peak.shape.1)
  f2 <- smoof::makeMPM2Function(n.peaks.2, dimensions, topology.2, seed.2, rotated.2, peak.shape.2)
  f3 <- smoof::makeMPM2Function(n.peaks.3, dimensions, topology.3, seed.3, rotated.3, peak.shape.3)
  
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
