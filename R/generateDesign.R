#' @export
generateDesign = function(fn, points.total = 1e6, points.per.dimension=NULL) {
  upper = smoof::getUpperBoxConstraints(fn)
  lower = smoof::getLowerBoxConstraints(fn)
  p = smoof::getNumberOfParameters(fn)
  
  if (!is.null(points.total) & !is.null(points.per.dimension)) {
    warning('points.per.dimension is set and will overwrite points.total')
  }
  
  l = list()
  step.sizes = c()
  
  if (!is.null(points.total)) {
    # Use points.total to calculate points.per.dimension
    ranges = upper - lower
    base.length = (points.total ** (1/p))
    factors = ranges / (prod(ranges) ** (1/p))
    points.per.dimension = round(base.length * factors)
  }
  
  if (length(points.per.dimension) == 1L) {
    points.per.dimension = rep(points.per.dimension, p)
  }
  
  for (i in 1:p) {
    x.i = c(paste0("x", i))
    l[[x.i]] = seq(lower[i], upper[i], length.out = points.per.dimension[i])
    step.sizes = c(step.sizes, (upper[i] - lower[i])/(points.per.dimension[i] - 1))
  }
  
  grid = list()
  grid$dec.space = as.matrix(expand.grid(l))
  grid$dims = points.per.dimension
  grid$step.sizes = step.sizes
  
  return(grid)
}
