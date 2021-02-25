#' @export
generateDesign = function(fn, points.total = NULL, points.per.dimension = NULL, upper = NULL, lower = NULL) {
  if (is.null(upper) && is.null(lower)) {
    upper = smoof::getUpperBoxConstraints(fn)
    lower = smoof::getLowerBoxConstraints(fn)
  }
  
  p = smoof::getNumberOfParameters(fn)
  
  if (!is.null(points.total) & !is.null(points.per.dimension)) {
    stop('Both points.per.dimension and points.total are set! Please choose only one parameter.')
  }
  
  if (is.null(points.total) & is.null(points.per.dimension)) {
    warning('Neither points.per.dimension nor points.total are set! Will use points.total = 1e6.')
    points.total = 1e6
  }
  
  l = list()
  step.sizes = c()
  
  if (!is.null(points.total) & is.null(points.per.dimension)) {
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
