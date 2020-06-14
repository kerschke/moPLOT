#' Performs a multi-objective local search.
#'
#' @description
#'   Given three initial points \code{x1}, \code{x2} and \code{x3}, of which is known that
#'   the angle between \code{x2 - x1} and \code{x3 - x2} is larger than 90 degree (i.e.,
#'   they point in opposite directions), this optimizer tries to find a local effient
#'   point that has to be located betwen \code{x1} and \code{x2}.
#'
#' @param x1 [\code{\link{numeric}(d)}]\cr
#'   d-dimensional individual located on one side of the (bi-objective) optimum.
#' @param x2 [\code{\link{numeric}(d)}]\cr
#'   d-dimensional individual located on the opposite side (w.r.t. \code{x1}) of
#'   the (bi-objective) optimum.
#' @param x3 [\code{\link{numeric}(d)}]\cr
#'   d-dimensional individual located such that the vector from \code{x2} to \code{x3}
#'   points in the opposite direction of the vector from \code{x1} to \code{x2}.
#' @template arg_fni
#' @template arg_precgrad
#' @template arg_precnorm
#' @template arg_precangle
#' @template arg_scalestep
#' @template arg_lower
#' @template arg_upper
#' @param max.steps [\code{\link{integer}(1L)}]
#'   Maximum number of local search steps to reach an optimum. The default is \code{1000L}.
#' @return [\code{\link{list}(5L)}]\cr
#'   List containing the matrix of points, which were visited in the course
#'   of the optimization, another matrix providing the number of performed
#'   function evaluations, two vectors providing the single-objective gradients
#'   of the last individual and a logical flag, indicating whether the optimizer
#'   found a local efficient point.
#' @examples
#' # Define two single-objective test problems:
#' fn1 = function(x) sum((x - c(2, 0))^2)
#' fn2 = function(x) sum((x - c(0, 1))^2)
#' fn = function(x) return(fn1(x), fn2(x))
#' 
#' # Visualize locally efficient set, i.e., the "area" where we ideally want to find a point:
#' plot(c(2, 0), c(0, 1), type = "o", pch = 19,
#'   xlab = expression(x[1]), ylab = expression(x[2]), las = 1, asp = 1)
#' text(2, 0, "Optimum of fn1", pos = 2, offset = 1.5)
#' text(0, 1, "Optimum of fn2", pos = 4, offset = 1.5)
#' 
#' # Place two points x1 and x2 on opposite sides of the bi-objective optimum:
#' x1 = c(1, 1)
#' x2 = c(0.5, 0)
#' x3 = c(0.8, 0.2)
#' points(rbind(x1, x2, x3), pch = 19, type = "o", lty = "dotted")
#' text(rbind(x1, x2, x3), labels = c("x1", "x2", "x3"), pos = 4)
#' 
#' # Optimize using weighted bisection optimization:
#' result = performMultiObjectiveLocalSearch(x1 = x1, x2 = x2, x3 = x3,
#'   fn = fn, lower = c(0, 0), upper = c(2, 1))
#' opt.path = result$opt.path
#' 
#' # Visualize the optimization path:
#' points(opt.path)
#' 
#' # Highlight the found local efficient point (= local optimum w.r.t. both objectives):
#' n = nrow(opt.path)
#' points(opt.path[n, 1], opt.path[n, 2], pch = 4, col = "red", cex = 2)
#' text(opt.path[n, 1], opt.path[n, 2], "Found Local Efficient Point", pos = 4, offset = 1.5)
#' @export
performMultiObjectiveLocalSearch = function(x1, x2, x3, fn,
  prec.grad = 1e-6, prec.norm = 1e-6, prec.angle = 1e-4,
  scale.step = 0.5, max.steps = 1000L) {

  d = getNumberOfParameters(fn)
  p = getNumberOfObjectives(fn)

  sp = seq_len(p)
  opt.path = matrix(x3, nrow = 1L)
  fn.evals = matrix(0L, nrow = 1L, ncol = 1L)
  gradient.mat = matrix(NA, nrow = p, ncol = d)
  
  while (nrow(opt.path) <= max.steps) {
    if (nrow(opt.path) > 1L) {
      i = nrow(opt.path)
      gradient.step = performGradientStep(ind = x2, fn = fn,
        gradient.mat = gradient.mat,
        scale.step = scale.step, prec.grad = prec.grad,
        prec.norm = prec.norm, prec.angle = prec.angle)
      x3 = gradient.step$offspring

      ## update function evaluations
      fn.evals[i,] = fn.evals[i,] + gradient.step$fn.evals[1L,]

      if (is.null(x3)) {
        ## if the x2 is a local efficient point, we can stop here without
        ## adding any new points to the optimization path
        gradient.mat = gradient.step$gradient.mat
        break
      } else if (identical(x2, x3)) {
        ## in the case, where x2 and x3 are identical, we are stuck in a trap;
        ## --> gradient.mat should be reset
        break
      } else {
        opt.path = rbind(opt.path, x3)
        fn.evals = rbind(fn.evals, 0L)
        if (nrow(opt.path) == max.steps + 1L) {
          break
        }
      }
    }

    ## compute angle between (x1, x2) and (x2, x3)
    v1 = x2 - x1
    v2 = x3 - x2
    angle = computeAngleCPP(vec1 = v1, vec2 = v2, prec = prec.norm)
    if (angle > 90) {
      ## if the vectors v1 and v2 show in oppsite directions, x2 must be
      ## located on the opposite side of the efficient set
      vl1 = computeVectorLengthCPP(v1)
      if (vl1 < (prec.norm / 100)) {
        ## in case x1 and x2 are extremely close to each other, the algorithm
        ## can not move forward --> abort
        break
      }
      ## reduce the step-size and place new point on connection of x1 and x2
      vl2 = computeVectorLengthCPP(v2)
      scale.step = (vl1 / (vl1 + vl2)) * scale.step
      x2 = x1 + (vl1 / (vl1 + vl2)) * v1
      opt.path = rbind(opt.path, x2)
      fn.evals = rbind(fn.evals, 0L)
    } else {
      ## if the angle is <= 90 degrees, rename x3 to x2 and x2 to x1 and proceed
      x1 = x2
      x2 = x3
    }
  }

  rownames(opt.path) = NULL
  rownames(fn.evals) = NULL
  return(list(
    opt.path = opt.path,
    fn.evals = fn.evals,
    gradient.mat = gradient.mat,
    found.optimum = is.null(x3)))
}
