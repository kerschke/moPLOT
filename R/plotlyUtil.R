plotlyColorscale = function(color.palette = NULL) {
  if (missing(color.palette)) {
    color.palette = viridisLite::turbo(500L, alpha = 1, begin = 0.05, end = 0.95, direction = 1)
  }
  
  l = length(color.palette)
  
  mapply(list, seq(0, 1, length.out = l), color.palette, SIMPLIFY = F)
}

imputeZero = function(height) {
  ## impute heights of zero for log-scale visualizations
  min.height = min(height[height != 0])
  height[height == 0] = min.height / 2
  return(height)
}

plotlyMarker = function(height, colorscale = plotlyColorscale()) {
  list(
    color = ~log(height),
    colorscale = colorscale,
    cmin = log(min(height)),
    cmax = log(max(height))
  )
}
