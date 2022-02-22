#' Compute the multi-objective gradient vector for a set of points.
#'
#' @description
#'   Computes the multi-objective gradients for a matrix of `points`.
#'
#' @param points [[matrix]]\cr
#'   Matrix of points, for which the multi-objective gradient should be computed.
#'   Each row of the matrix will be considered as a separate point, thus the number
#'   of rows corresponds to the number of observations and the number of columns
#'   to the dimensionality of the search space.
#' @template arg_fn
#' @template arg_precgrad
#' @template arg_precnorm
#' @template arg_precangle
#' @param parallelize [[`logical`]]\cr
#'   Should the computation of the gradient vectors be parallelized (with `parallel::mclapply`)?
#'   The default is `FALSE`.
#' @param impute.boundary [[`logical`]]\cr
#'   Should values at the boundary be imputed? This is relevant for problems
#'   that have locally efficient points on the boundary box.
#' @template arg_lower
#' @template arg_upper
#' @return [[`matrix`]]\cr
#'   Returns `matrix` of multi-objective gradients. The i-th row of the matrix
#'   contains the multi-objective gradient vector of the i-th observation (= row)
#'   of the input matrix `points`.
#' @examples
#' 
#' @export
computeGradientField = function(points, fn, prec.grad = 1e-6,
  prec.norm = 1e-6, prec.angle = 1e-4, parallelize = FALSE, impute.boundary = TRUE, lower = NULL, upper = NULL) {

  cat("Estimating single-objective gradients ...\n")

  if (smoof::isSmoofFunction(fn)) {
    lower = smoof::getLowerBoxConstraints(fn)
    upper = smoof::getUpperBoxConstraints(fn)
  }

  if (parallelize) {
    gradients.list = parallel::mclapply(seq_row(points), function(i) {
      ind = as.numeric(points[i,])

      -estimateGradientBothDirections(fn = fn, ind = ind, prec.grad = prec.grad, check.data = FALSE, lower = lower, upper = upper)
    })
  } else {
    gradients.list = lapply(seq_row(points), function(i) {
      ind = as.numeric(points[i,])

      -estimateGradientBothDirections(fn = fn, ind = ind, prec.grad = prec.grad, check.data = FALSE, lower = lower, upper = upper)
    })
  }

  # Convert to matrix per objective (in a list)
  # With one row per gradient each
  single.objective = simplify2array(gradients.list)
  single.objective = asplit(single.objective, 1)
  single.objective = lapply(single.objective, t)

  obj = length(single.objective)

  cat("Estimating multi-objective gradients ...\n")
  if (obj == 2L) {
    multi.objective = getBiObjGradientGridCPP(gradMat1 = single.objective[[1L]], gradMat2 = single.objective[[2L]], precNorm = prec.norm, precAngle = prec.angle)
  } else if (obj == 3L) {
    multi.objective = getTriObjGradientGridCPP(gradMat1 = single.objective[[1L]], gradMat2 = single.objective[[2L]], gradMat3 = single.objective[[3]], precNorm = prec.norm, precAngle = prec.angle)
  } else {
    stop("Cannot cannot handle more than 3 objectives.")
  }

  cat("Finished multi-objective gradients\n")

  if (impute.boundary) {
    cat("Imputing multi-objective gradient at boundary\n")
    dims = apply(points, 2L, function(x) length(unique(x)))
    multi.objective = imputeBoundary(multi.objective, single.objective, dims)
  }

  return(list(
    multi.objective=multi.objective,
    single.objective=single.objective
  ))
}

#' @export
computeGradientFieldGrid = function(grid, prec.norm = 1e-6, prec.angle = 1e-4, impute.boundary = TRUE, normalized.scale = TRUE) {

  obj = ncol(grid$obj.space)
  
  assertList(grid, min.len = 3L, names = "named")
  assertSubset(c("dims", "step.sizes", "obj.space"), choices = names(grid))
  # currently need normalization, if number of objectives != 2
  assertTRUE(normalized.scale || obj == 2)
  
  cat("Estimating single-objective gradients ...\n")

  single.objective = lapply(seq_len(obj), function(i) {
    cat(paste("Differentiating objective", i, "\n"))
    -gridBasedGradientCPP(grid$obj.space[,i], grid$dims, grid$step.sizes, prec.norm, prec.angle)
  })
  
  cat("Estimating multi-objective gradients ...\n")
  
  if (obj == 2L) {
    multi.objective = getBiObjGradientGridCPP(gradMat1 = single.objective[[1L]], gradMat2 = single.objective[[2L]], precNorm = prec.norm, precAngle = prec.angle, normalized_scale = normalized.scale)
  } else if (obj == 3L) {
    multi.objective = getTriObjGradientGridCPP(gradMat1 = single.objective[[1L]], gradMat2 = single.objective[[2L]], gradMat3 = single.objective[[3L]], precNorm = prec.norm, precAngle = prec.angle)
  } else {
    stop("Cannot cannot handle more than 3 objectives.")
  }

  cat("Finished multi-objective gradients\n")
  
  if (impute.boundary) {
    cat("Imputing multi-objective gradient at boundary\n")
    
    multi.objective = imputeBoundary(multi.objective, single.objective, grid$dims, normalized_scale = normalized.scale)
  }
  
  return(list(
    multi.objective = multi.objective,
    single.objective = single.objective
  ))
}

#' @export
computeDivergenceGrid = function(gradients, dims, step.sizes, prec.norm = 1e-6, prec.angle = 1e-4, normalize = FALSE) {
  if (normalize) {
    gradients = normalizeMatrixRowsCPP(gradients, prec.norm)
  }
  
  d = length(dims)

  # TODO There has to be some extra handling to catch all cases on the decision boundary

  if (d == 2) {
    l = lapply(seq_len(ncol(gradients)), function(i) {
      gridBasedGradientCPP(gradients[, i], dims, step.sizes, prec.norm, prec.angle)[, i]
    })
    
    Reduce('+', l)
  } else if (d == 3) {
    l = lapply(seq_len(ncol(gradients)), function(i) {
      gridBasedGradientCPP(gradients[, i], dims, step.sizes, prec.norm, prec.angle)
    })
    
    sapply(seq_len(nrow(gradients)), function(i) {
      J = matrix(
        c(l[[1L]][i, 1L], l[[1L]][i, 2L], l[[1L]][i, 3L],
          l[[2L]][i, 1L], l[[2L]][i, 2L], l[[2L]][i, 3L],
          l[[3L]][i, 1L], l[[3L]][i, 2L], l[[3L]][i, 3L]),
        nrow = 3L, ncol = 3L, byrow = TRUE
      )
      
      real_parts <- Re(eigen(1/2 * (J + t(J)), symmetric = TRUE, only.values = TRUE)$values)
      
      pair_sums = c(
        real_parts[1] + real_parts[2],
        real_parts[1] + real_parts[3],
        real_parts[2] + real_parts[3]
      )
      
      if (all(pair_sums < 0)) {
        max(pair_sums)
      } else if (all(pair_sums > 0)) {
        min(pair_sums)
      } else {
        median(pair_sums)
      }
      
      # See: 
      # Chong, M.S., Perry, A.E., Cantwell, B.J.:
      # A general classification of three-dimensional flow fields.
      # Physics of Fluids A: Fluid Dynamics 2(5), 765–777 (1990)
      
      # (For the time) too inaccurate,
      # as we expect some eigenvalues to be approx. zero:
      
      # P = -(J[1,1] + J[2,2] + J[3,3])
      # Q = (J[1,1] * J[2,2] - J[2,1] * J[1,2]) +
      #     (J[1,1] * J[3,3] - J[3,1] * J[1,3]) +
      #     (J[2,2] * J[3,3] - J[3,2] * J[2,3])
      # R = -det(J)

      # if (P > 0) {
      #   if (R < 0) {
      #     -1 # stable
      #   } else if (R > 0 && Q <= R / P) {
      #     1 # unstable
      #   } else {
      #     -1 # stable
      #   }
      # } else if (P < 0) {
      #   if (R > 0) {
      #     1 # unstable
      #   } else if (R < 0 && Q <= R / P) {
      #     -1 # stable
      #   } else {
      #     1 # unstable
      #   }
      # } else {
      #   # P == 0
      #   if (R < 0) {
      #     -1 # stable
      #   } else if (R > 0) {
      #     1 # unstable
      #   } else {
      #     -1 # stable
      #   }
      # }
    })
  }
}

#' @export
computeSecondOrderGrid = function(gradients, dims, step.sizes, prec.norm = 1e-6, prec.angle = 1e-4, normalize = FALSE, epsilon=1e-2, get.matrix = TRUE) {
  if (normalize) {
    gradients = normalizeMatrixRowsCPP(gradients, prec.norm)
  }

  if (length(dims) != 2L)
    return(NULL)

  l = lapply(seq_len(ncol(gradients)), function(i) {
    gridBasedGradientCPP(gradients[,i], dims, rep(1L, length(dims)), prec.norm, prec.angle)
  })

  if (get.matrix) {
    j = sapply(seq_len(nrow(gradients)), function(i) {
      matrix(c(l[[1L]][i, 1L], l[[1L]][i, 2L], l[[2L]][i, 1L], l[[2L]][i, 2L]), nrow = 2L, ncol = 2L, byrow = TRUE)
    })
    return(t(j))
  } else {
    j = sapply(seq_len(nrow(gradients)), function(i) {
      H = matrix(c(l[[1L]][i, 1L], l[[1L]][i, 2L], l[[2L]][i, 1L], l[[2L]][i, 2L]), nrow = 2L, ncol = 2L, byrow = TRUE)
      H = (t(H) + H) / 2

      Re(eigen(H, only.values = TRUE)$values)
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
