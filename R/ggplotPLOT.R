#' @export
ggplotPLOT = function(dec.space, obj.space, sinks, height) {
  
  nds = ecr::doNondominatedSorting(t(obj.space[sinks,]))
  
  height[height == 0] = min(height[height > 0]) / 2
  
  height.df = cbind.data.frame(dec.space, height=height)
  height.df = height.df[order(height.df$height, decreasing=T),]
  
  sinks.df = cbind.data.frame(dec.space[sinks,], rank=nds$dom.counter+1)
  sinks.df = sinks.df[order(sinks.df$rank, decreasing=T),]
  
  g = ggplotHeatmap(height.df, color.palette = gray.colors(500L, start = 0, end = 1, gamma = 0.5)) +
    geom_point(mapping = aes(x1, x2, color=rank), data = sinks.df) +
    scale_color_gradientn(colors = fields::tim.colors(500L), na.value="black", trans = "log") +
    theme_minimal() +
    theme(legend.position = "none", axis.title = element_blank())
  
  return(g)
}
