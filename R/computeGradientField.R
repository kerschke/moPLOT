#' Compute the multi-objective gradient vector for a set of points.
#'
#' @description
#'   Computes the multi-objective gradients for a matrix of \code{points}.
#'
#' @param points [\code{\link{matrix}}]\cr
#'   Matrix of points, for which the multi-objective gradient should be computed.
#'   Each row of the matrix will be considered as a separate point, thus the number
#'   of rows corresponds to the number of observations and the number of columns
#'   to the dimensionality of the search space.
#' @template arg_fni
#' @param fn3 [\code{\link{function}}]\cr
#'   The third objective used for computing the multi-objective gradient. If not provided
#'   (\code{fn3 = NULL}), the gradient field will only consider fn1 and fn2.
#' @template arg_scalestep
#' @template arg_precgrad
#' @template arg_precnorm
#' @template arg_precangle
#' @template arg_lower
#' @template arg_upper
#' @param parallelize [\code{\link{logical}(1L)}]\cr
#'   Should the computation of the gradient vectors be parallelized (with \code{parallel::mclapply})?
#'   The default is \code{FALSE}.
#' @return [\code{\link{matrix}}]\cr
#'   Returns \code{matrix} of multi-objective gradients. The i-th row of the matrix
#'   contains the multi-objective gradient vector of the i-th observation (= row)
#'   of the input matrix \code{points}.
#' @examples
#' # Define bi-objective test problems:
#' fn = function(x) c(sum((x - c(0.2, 1))^2), sum(x))
#'
#' # Create a grid of points, for which the gradients should be computed:
#' points = expand.grid(x1 = seq(0, 1, 0.01), x2 = seq(0, 1, 0.05))
#' gradient.field = computeGradientField(points, fn)
#' @export
computeGradientField = function(points, fn, prec.grad = 1e-6,
  prec.norm = 1e-6, prec.angle = 1e-4, parallelize = FALSE) {
  if (parallelize) {
    r = parallel::mclapply(seq_row(points), function(i) {
      ind = as.numeric(points[i,])
      return(calcMOGradient(ind, fn, prec.grad, prec.norm, prec.angle))
    })
  } else {
    r = lapply(seq_row(points), function(i) {
      ind = as.numeric(points[i,])
      return(calcMOGradient(ind, fn, prec.grad, prec.norm, prec.angle))
    })
  }

  return(as.matrix(do.call(rbind, r)))
}

#' @export
computeGradientFieldGrid = function(grid, fn, prec.norm = 1e-6, prec.angle = 1e-4) {
  obj = smoof::getNumberOfObjectives(fn)
  d = ncol(points) # number of input dimensions
  n = nrow(points) # total number of points
  # calculate dimensions of given field of points

  if (is.null(grid$obj.space)) {
    cat("Evaluating grid of objective values ...\n")
    grid$obj.space = calculateObjectiveValues(grid$dec.space, fn)
  }

  cat("Estimating single-objective gradients ...\n")

  grad.mat.1 = -gridBasedGradientCPP(grid$obj.space[,1], grid$dims, grid$step.sizes, prec.norm, prec.angle)
  cat("Finished objective 1\n")
  grad.mat.2 = -gridBasedGradientCPP(grid$obj.space[,2], grid$dims, grid$step.sizes, prec.norm, prec.angle)
  cat("Finished objective 2\n")

  if (obj == 2) {
    cat("Estimating multi-objective gradients ...\n")

    mo.grad.mat = getBiObjGradientGridCPP(grad.mat.1, grad.mat.2, prec.norm, prec.angle)
    cat("Finished multiobjective gradients\n")
    return(list(
      multi.objective=mo.grad.mat,
      single.objective=list(grad.mat.1, grad.mat.2)
    ))
  }

  if (obj == 3) {
    grad.mat.3 = -gridBasedGradientCPP(grid$obj.space[,3], grid$dims, grid$step.sizes, prec.norm, prec.angle)
    cat("Finished objective 3\n")

    cat("Estimating multi-objective gradients ...\n")

    mo.grad.mat = getTriObjGradientGridCPP(grad.mat.1, grad.mat.2, grad.mat.3, prec.norm, prec.angle)
    cat("Finished multiobjective gradients\n")
    return(list(
      multi.objective=mo.grad.mat,
      single.objective=list(grad.mat.1, grad.mat.2, grad.mat.3)
    ))
  }

}

#' @export
computeDivergenceGrid = function(gradients, dims, step.sizes, prec.norm = 1e-6, prec.angle = 1e-4, normalize = F) {
  if (normalize) {
    gradients = normalizeMatrixRowsCPP(gradients, prec.norm)
  }

  l = lapply(1:ncol(gradients), function(i) {
    gridBasedGradientCPP(gradients[,i], dims, step.sizes, prec.norm, prec.angle)[,i]
  })

  Reduce('+', l)
}

#' @export
computeSecondOrderGrid = function(gradients, dims, step.sizes, prec.norm = 1e-6, prec.angle = 1e-4, normalize = F, epsilon=1e-2, get.matrix = T) {
  if (normalize) {
    gradients = normalizeMatrixRowsCPP(gradients, prec.norm)
  }

  if (length(dims) != 2) return()

  l = lapply(1:ncol(gradients), function(i) {
    gridBasedGradientCPP(gradients[,i], dims, rep(1, length(dims)), prec.norm, prec.angle)
  })

  if (get.matrix) {
    j = sapply(1:nrow(gradients), function(i) {
      matrix(c(l[[1]][i, 1], l[[1]][i, 2], l[[2]][i, 1], l[[2]][i, 2]), nrow = 2, ncol = 2, byrow=T)
    })
    return(t(j))
  } else {
    j = sapply(1:nrow(gradients), function(i) {
      H = matrix(c(l[[1]][i, 1], l[[1]][i, 2], l[[2]][i, 1], l[[2]][i, 2]), nrow = 2, ncol = 2, byrow=T)
      H = (t(H) + H) / 2

      Re(eigen(H, only.values = T)$values)
    })
    return(t(j))
  }
}


#' @export
genericMOGradient = function(G, prec.norm = 1e-6) {
  # usage: e.g. genericMOGradient(cbind(c(1,0,0), c(0,1,0), rep(sqrt(1/3), 3)))
  G = t(normalizeMatrixRowsCPP(t(G), prec.norm))
  pracma::qpspecial(G)$d
}

calcMOGradient = function(ind, fn, prec.grad, prec.norm, prec.angle) {
  g = -estimateGradientBothDirections(fn = fn, ind = ind, prec.grad = prec.grad, check.data = FALSE)
  if (nrow(g) < 3) {
    getBiObjGradientCPP(g[1,], g[2,], prec.norm, prec.angle)
  } else {
    getTriObjGradientCPP(g[1,], g[2,], g[3,], prec.norm, prec.angle)
  }
}
