#' @export
ggplotPLOT = function(dec.space, obj.space, sinks, height, check.data = TRUE) {
  if (check.data) {
    assertMatrix(dec.space, ncols = 2, col.names = "unique") # 2D only
    assertMatrix(obj.space, min.cols = 2, nrows = nrow(dec.space)) # at least 2-objective
    
    assertInteger(sinks, lower = 1, upper = nrow(dec.space)) # sinks are valid indices for observations
    assertNumeric(height, len = nrow(dec.space)) # one height value per observation
  }
  
  # Compute coloring of locally efficient points
  
  nds = ecr::doNondominatedSorting(t(obj.space[sinks,]))
  
  # Impute Gradient Field Heatmap height
  
  height[height == 0] = min(height[height > 0]) / 2
  
  height.df = cbind.data.frame(dec.space, height=height)
  
  sinks.df = cbind.data.frame(dec.space[sinks,], rank=nds$dom.counter+1)
  sinks.df = sinks.df[order(sinks.df$rank, decreasing=T),] # Ensure more efficient points are plotted later, i.e. on top
  
  g = ggplotHeatmap(height.df, color.palette = gray.colors(500L, start = 0, end = 1, gamma = 0.5)) +
    geom_point(mapping = aes(x1, x2, color=rank), data = sinks.df) +
    scale_color_gradientn(colors = fields::tim.colors(500L), na.value="black", trans = "log") +
    theme_minimal() +
    theme(legend.position = "none", axis.title = element_blank())
  
  return(g)
}
