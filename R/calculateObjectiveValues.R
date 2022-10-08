#' @export
calculateObjectiveValues = function(
    points, fn, parallelize = FALSE, parallel.cores = (parallel::detectCores() - 1),
    vectorized.evaluation = TRUE, convert.to.minimization = TRUE, verbose = TRUE) {
  n = smoof::getNumberOfObjectives(fn)
  names = paste0("y", 1:n)
  
  if (verbose) cat("Evaluating grid of objective values ...\n")
  
  if (vectorized.evaluation && smoof::isVectorized(fn)) {
    # Experimental feature!
    obj.space = fn(t(points))
  } else if (parallelize) {
    cuts = cut(1:nrow(points), parallel.cores)
    r = parallel::mclapply(levels(cuts), function(i) {
      t(apply(points[cuts == i,], 1, fn))
    }, mc.cores = parallel.cores)
    
    # TODO: Which one is really faster?
    # r = parallel::mclapply(seq_row(points), function(i) {
    #   as.numeric(fn(points[i,]))
    # }, mc.cores = parallel.cores)
    
    obj.space = as.matrix(do.call(rbind, r))
  } else {
    obj.space = t(apply(points, 1, fn))
  }
  
  if (convert.to.minimization) {
    minimize_flag = smoof::shouldBeMinimized(fn)
    for (i in seq_along(minimize_flag)) {
      obj.space[,i] <- obj.space[,i] * (if (minimize_flag[i]) 1 else -1)
    }
  }
  
  colnames(obj.space) = names
  
  return(obj.space)
}
