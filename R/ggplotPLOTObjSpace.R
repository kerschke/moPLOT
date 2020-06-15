#' @export
ggplotPLOTObjSpace = function(dec.space, obj.space, sinks, height) {
  
  nds = ecr::doNondominatedSorting(t(obj.space[sinks,]))
  
  height.obj.df = cbind.data.frame(obj.space, height=height)
  sinks.obj.df = cbind.data.frame(obj.space[sinks,], height=nds$dom.counter+1)
  height.obj.df = height.obj.df[order(height.obj.df$height, decreasing=T),]
  sinks.obj.df = sinks.obj.df[order(sinks.obj.df$height, decreasing=T),]
  
  g = ggplot() +
    geom_point(mapping = aes(y1, y2, color=height), data = height.obj.df) +
    geom_point(mapping = aes(y1, y2, fill=height, color=NA), data = sinks.obj.df, shape=21) +
    scale_color_gradientn(colors = gray.colors(500L, start = 0, end = 1, gamma = 0.5), na.value="transparent", trans="log") +
    scale_fill_gradientn(colors = fields::tim.colors(500L), trans = "log") +
    xlab(expression(y[1])) +
    ylab(expression(y[2])) +
    theme_minimal() +
    theme(legend.position = "none")  
  
  return(g)
}
