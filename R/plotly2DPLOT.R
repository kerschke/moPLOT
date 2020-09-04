#' @export
plotly2DPLOT = function(dec.space, obj.space, sinks, height, fn, mode = "decision.space", impute.zero = T) {
  # fn: smoof function, 2 dimensional decision space
  
  n = smoof::getNumberOfObjectives(fn)
  lower = smoof::getLowerBoxConstraints(fn)
  upper = smoof::getUpperBoxConstraints(fn)
  
  if (impute.zero) {
    height = imputeZero(height)
  }
  
  domination.counts = ecr::doNondominatedSorting(t(obj.space[sinks,]))$dom.counter + 1
  
  x = cbind.data.frame(dec.space, height, obj.space)
  
  x.locally.efficient = x[sinks,]
  x.locally.efficient$domination.counts = domination.counts
  
  x.order = order(x$height, decreasing = TRUE)
  x = x[x.order,] # relevant for obj.space
  x.locally.efficient = x.locally.efficient[order(x.locally.efficient$domination.counts, decreasing = TRUE),]
  
  if (n == 3) {
    objective.scene = list(
      aspectmode='cube',
      xaxis = list(range = c(min(x[,'y1']),max(x[,'y1'])), title='y₁'),
      yaxis = list(range = c(min(x[,'y2']),max(x[,'y2'])), title='y₂'),
      zaxis = list(range = c(min(x[,'y3']),max(x[,'y3'])), title='y₃')
    )
  }
  
  marker.heatmap = plotlyMarker(height, plotlyColorscale(gray.colorscale))
  marker.locally.efficient = list(
    color=~log(domination.counts),
    colorscale=plotlyColorscale(),
    cmin=log(min(domination.counts)),
    cmax=log(max(domination.counts))
  )
  
  if (mode == "both") {
    x.le.shared = highlight_key(x.locally.efficient)
    
    p.decision = plotly2DPLOTDecisionSpace(x, x.le.shared, fn, marker.locally.efficient)
    p.objective = plotly2DPLOTObjectiveSpace(x, x.le.shared, fn, marker.heatmap, marker.locally.efficient, scene="scene")
    
    domain.left = list(
      x=c(0,0.5),
      y=c(0,1)
    )
    decision.scene = list(domain=domain.left)
    
    domain.right = list(
      x=c(0.5,1),
      y=c(0,1)
    )
    if (n == 3) {
      objective.scene$domain = domain.right
    } else {
      objective.scene = list(domain=domain.right)
    }
    
    subplot(p.decision, p.objective) %>% layout(
      scene = objective.scene
    ) %>% highlight(
      on="plotly_click",
      off="plotly_deselect",
      opacityDim = 0.5,
      color = "red"
    ) %>% hide_guides()
  } else if (mode == "decision.space") {
    plotly2DPLOTDecisionSpace(x, x.locally.efficient, fn, marker.locally.efficient) %>%
      hide_guides()
  } else if (mode == "objective.space") {
    if (n == 3) {
      plotly2DPLOTObjectiveSpace(x, x.locally.efficient, fn, marker.heatmap, marker.locally.efficient) %>% layout(
        scene = objective.scene
      ) %>% hide_guides()
    } else {
      plotly2DPLOTObjectiveSpace(x, x.locally.efficient, fn, marker.heatmap, marker.locally.efficient) %>%
        hide_guides()
    }
  }
  
}

plotly2DPLOTObjectiveSpace = function(x, x.locally.efficient, fn, marker.heatmap, marker.locally.efficient, scene="scene") {
  n = smoof::getNumberOfObjectives(fn)
  
  if (n == 2) {
    plot_ly(
      data = x,
      type="scatter",
      x = ~y1, y = ~y2,
      mode = "markers",
      marker = marker.heatmap
    ) %>% add_trace(
      data = x.locally.efficient,
      type = "scatter",
      x = ~y1, y = ~y2,
      mode = "markers",
      marker = marker.locally.efficient,
      inherit = FALSE
    ) %>% layout(
      xaxis = list(
        title = "y₁",
        constrain = "domain"
      ),
      yaxis = list(
        title = "y₂",
        constrain = "domain"
      )
    )
  } else if (n == 3) {
    plot_ly(
      data = x,
      type="scatter3d",
      x = ~y1, y = ~y2, z = ~y3,
      scene = scene,
      mode = "markers",
      marker = marker.heatmap
    ) %>% add_trace(
      data = x.locally.efficient,
      type = "scatter3d",
      x = ~y1, y = ~y2, z = ~y3,
      mode = "markers",
      marker = marker.locally.efficient,
      inherit = FALSE
    )
  }
}

plotly2DPLOTDecisionSpace = function(x, x.locally.efficient, fn, marker.locally.efficient) {
  plot_ly(data = x,
          type = "heatmap",
          x = ~x1, y = ~x2, z = ~log(height),
          colorscale = plotlyColorscale(gray.colorscale)
  ) %>% add_trace(
    data = x.locally.efficient,
    type = "scattergl",
    x = ~x1, y = ~x2,
    mode = "markers",
    marker = marker.locally.efficient,
    inherit = FALSE
  ) %>% layout(
    xaxis = list(
      title = "x₁",
      constrain = "domain"
    ),
    yaxis = list(
      scaleanchor = "x",
      scaleratio = (max(x$x1) - min(x$x1)) / (max(x$x2) - min(x$x2)),
      title = "x₂",
      constrain = "domain"
    )
  )
}
