#' Create a PLOT Visualization
#' 
#' Returns a PLOT visualization as a [ggplot2::ggplot] object for the given visualization data.
#' 
#' @param dec.space 
#' Numeric [base::matrix] that defines the evaluated points in decision space.
#' One column per dimension.
#' @param obj.space 
#' Numeric [base::matrix] that defines the evaluated points in objective space.
#' One column per objective.
#' @param sinks
#' Integer [base::vector] of (row) indices that identify the _locally efficient_ points.
#' @param height 
#' Numeric [base::vector] that assigns a height value to each evaluated point.
#' @template arg_checkdata
#'
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
  
  var1 = colnames(dec.space)[1]
  var2 = colnames(dec.space)[2]
  
  colorscale.efficient = fields::tim.colors(500L)
  colorscale.heatmap = gray.colors(500L, start = 0, end = 1, gamma = 0.5)
  
  if (all(sinks.df$rank == sinks.df$rank[1])) {
    # If all "rank" values are identical, they are all globally efficient.
    # Change the colorscale so that they are plotted as the "lowest" color
    # value rather than the middle one (which is the default).
    colorscale.efficient = colorscale.efficient[1]
  }
  
  g = ggplotHeatmap(height.df, color.palette = colorscale.heatmap,
                    var1 = var1, var2 = var2) +
    geom_point(mapping = aes_string(var1, var2, color="rank"), data = sinks.df) +
    scale_color_gradientn(colors = colorscale.efficient, na.value="black", trans = "log") +
    theme(legend.position = "none")
  
  return(g)
}
