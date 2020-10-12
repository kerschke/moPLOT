plotlyColorscale = function(color.palette = NULL) {
  if (missing(color.palette)) {
    ## if no information on the color palette is provided, this
    ## function tries to sequentially tries to use tim.colors,
    ## viridis or at last the terrain.colors
    inst.pkgs = rownames(installed.packages())
    if ("fields" %in% inst.pkgs) {
      color.palette = fields::tim.colors(500L)
    } else if ("viridisLite" %in% inst.pkgs) {
      color.palette = viridisLite::viridis(500L,
                                           alpha = 1, begin = 0, end = 1, direction = 1, option = "D")
    } else {
      color.palette = terrain.colors(500L)
    }
  }
  
  l = length(color.palette)
  
  mapply(list, seq(0,1,length.out=l), color.palette, SIMPLIFY = F)
}

imputeZero = function(height) {
  ## impute heights of zero for log-scale visualizations
  min.height = min(height[height != 0])
  height[height == 0] = min.height / 2
  return(height)
}

plotlyMarker = function(height, colorscale = plotlyColorscale()) {
  list(
    color=~log(height),
    colorscale=colorscale,
    cmin=log(min(height)),
    cmax=log(max(height))
  )
}
