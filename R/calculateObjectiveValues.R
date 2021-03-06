#' @export
calculateObjectiveValues = function(points, fn, parallelize = FALSE, parallel.cores = (parallel::detectCores() - 1)) {
  n = smoof::getNumberOfObjectives(fn)
  names = paste0("y", 1:n)
  
  cat("Evaluating grid of objective values ...\n")
  
  if (parallelize) {
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
  
  colnames(obj.space) = names
  
  return(obj.space)
}
