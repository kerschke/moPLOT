#' @export
ggplotPLOTObjSpace = function(obj.space, sinks, height, check.data = TRUE) {
  if (check.data) {
    assertMatrix(obj.space, ncols = 2) # 2-objective only
    
    assertInteger(sinks, lower = 1, upper = nrow(obj.space)) # Sinks are valid indices for observations
    assertNumeric(height, len = nrow(obj.space)) # One height value per observation
  }
  
  # Compute coloring of locally efficient points
  
  nds = ecr::doNondominatedSorting(t(obj.space[sinks,]))
  
  # Impute Gradient Field Heatmap height
  
  height[height == 0] = min(height[height > 0]) / 2
  
  height.obj.df = cbind.data.frame(obj.space, height=height)
  height.obj.df = height.obj.df[order(height.obj.df$height, decreasing=T),]
  
  sinks.obj.df = cbind.data.frame(obj.space[sinks,], height=nds$dom.counter+1)
  sinks.obj.df = sinks.obj.df[order(sinks.obj.df$height, decreasing=T),]
  
  var1 = colnames(obj.space)[1]
  var2 = colnames(obj.space)[2]
  
  g = ggplot() +
    geom_point(mapping = aes_string(var1, var2, color="height"), data = height.obj.df) +
    geom_point(mapping = aes_string(var1, var2, fill="height", color=NA), data = sinks.obj.df, shape=21) +
    scale_color_gradientn(colors = gray.colors(500L, start = 0, end = 1, gamma = 0.5), na.value="transparent", trans="log") +
    scale_fill_gradientn(colors = fields::tim.colors(500L), trans = "log") +
    xlab(variable.as.expression(var1)) +
    ylab(variable.as.expression(var2)) +
    theme_minimal() +
    theme(legend.position = "none")
  
  return(g)
}
