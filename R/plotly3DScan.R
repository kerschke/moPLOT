#' @export
plotly3DScan = function(grid, fn, sinks = NULL, mode = "decision.space", frame = "x3", impute.zero = T) {
  # grid: list of obj.space, dims, dec.space, step.sizes
  # fn: smoof function, 3 dimensional decision space
  
  n = smoof::getNumberOfObjectives(fn)
  lower = smoof::getLowerBoxConstraints(fn)
  upper = smoof::getUpperBoxConstraints(fn)
  
  if (impute.zero) {
    grid$height = imputeZero(grid$height)
  }
  
  if (!is.null(sinks)) {
    dom.counter = ecr::doNondominatedSorting(t(grid$obj.space[sinks,]))$dom.counter
  }
  
  decision.scene = list(
    aspectmode='cube',
    xaxis = list(range = c(lower[1],upper[1]), title='x₁'),
    yaxis = list(range = c(lower[2],upper[2]), title='x₂'),
    zaxis = list(range = c(lower[3],upper[3]), title='x₃')
  )
  
  x = cbind.data.frame(grid$dec.space, grid$obj.space, height = grid$height)
  
  if (n == 3) {
    objective.scene = list(
      aspectmode='cube',
      xaxis = list(range = c(min(x$y1),max(x$y1)), title='y₁'),
      yaxis = list(range = c(min(x$y2),max(x$y2)), title='y₂'),
      zaxis = list(range = c(min(x$y3),max(x$y3)), title='y₃')
    )
  } else {
    objective.scene = list()
  }
  
  if (!is.null(sinks)) {
    x.sinks = x[sinks,]
    x.sinks.shared = highlight_key(x.sinks)
    dom.height = log(dom.counter + 1)
    marker.sinks = list(
      color = dom.height,
      colorscale = plotlyColorscale(fields::tim.colors(500L)),
      cmin = min(dom.height),
      cmax = max(dom.height)
    )
    
    x.heatmap = x[-sinks,]
    
    heatmap.order = switch (
      frame,
      "x1" = order(x.heatmap$x1),
      "x2" = order(x.heatmap$x2),
      "x3" = order(x.heatmap$x3)
    )
    
    x.heatmap = x.heatmap[heatmap.order,]
    x.heatmap.shared = highlight_key(x.heatmap)
    marker.heatmap = plotlyMarker(x.heatmap$height, colorscale = plotlyColorscale(gray.colorscale))
  } else {
    x.sinks = NULL
    x.sinks.shared = NULL
    marker.sinks = NULL
    
    x.heatmap = x
    
    heatmap.order = switch (
      frame,
      "x1" = order(x.heatmap$x1),
      "x2" = order(x.heatmap$x2),
      "x3" = order(x.heatmap$x3)
    )
    
    x.heatmap = x.heatmap[heatmap.order,]
    x.heatmap.shared = highlight_key(x.heatmap)
    marker.heatmap = plotlyMarker(x.heatmap$height, colorscale = plotlyColorscale(fields::tim.colors(500L)))
  }
  
  if (mode == "both") {
    p.decision = plotly3DScanDecisionSpace(x.heatmap.shared, marker.heatmap, x.sinks.shared, marker.sinks, frame = frame, scene = "scene")
    
    p.objective = plotly3DScanObjectiveSpace(fn, x.heatmap.shared, marker.heatmap, x.sinks.shared, marker.sinks, frame = frame, scene = "scene2")
    
    domain.left = list(
      x=c(0,0.5),
      y=c(0,1)
    )
    decision.scene$domain = domain.left
    
    domain.right = list(
      x=c(0.5,1),
      y=c(0,1)
    )
    objective.scene$domain = domain.right
    
    p = subplot(p.decision, p.objective) %>% layout(
      scene = decision.scene,
      scene2 = objective.scene
    )
  } else if (mode == "decision.space") {
    p = plotly3DScanDecisionSpace(x.heatmap, marker.heatmap, x.sinks, marker.sinks, frame = frame) %>% layout(
      scene = decision.scene
    )
  } else if (mode == "objective.space") {
    if (n == 3) {
      p = plotly3DScanObjectiveSpace(fn, x.heatmap, marker.heatmap, x.sinks, marker.sinks, frame = frame) %>% layout(
        scene = objective.scene
      )
    } else {
      p = plotly3DScanObjectiveSpace(fn, x.heatmap, marker.heatmap, x.sinks, marker.sinks, frame = frame) %>% layout(
        xaxis = list(range = c(min(x$y1),max(x$y1)), title='y₁'),
        yaxis = list(range = c(min(x$y2),max(x$y2)), title='y₂')
      )
    }
  }
  
  p %>% animation_opts(
    frame = 1000,
    transition = 0
  ) %>% hide_guides()
}

plotly3DScanObjectiveSpace = function(fn, x.heatmap, marker.heatmap, x.sinks = NULL, marker.sinks = NULL, frame="x3", scene="scene") {
  n = smoof::getNumberOfObjectives(fn)
  
  if (frame == "x1") {
    frame = ~x1
  } else if (frame == "x2") {
    frame = ~x2
  } else if (frame == "x3") {
    frame = ~x3
  }
  
  if (n == 2) {
    p = plot_ly() %>% add_markers(
      data = x.heatmap,
      type = "scattergl",
      x = ~y1, y = ~y2,
      frame = frame,
      mode = "markers",
      marker = marker.heatmap
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
    
    if (!is.null(x.sinks)) {
      p = p %>% add_markers(
        type="scattergl",
        data = x.sinks,
        x = ~y1, y = ~y2,
        mode = "markers",
        marker = marker.sinks
      )
    }
  } else if (n == 3) {
    p = plot_ly(
      scene = scene
    ) %>% add_markers(
      data = x.heatmap,
      type = "scatter3d",
      x = ~y1, y = ~y2, z = ~y3,
      frame = frame,
      mode = "markers",
      marker = marker.heatmap
    )
    
    if (!is.null(x.sinks)) {
      p = p %>% add_markers(
        type = "scatter3d",
        x = ~y1, y = ~y2, z = ~y3,
        data = x.sinks,
        mode = "markers",
        marker = marker.sinks
      )
    }
  }
    
  p
}

plotly3DScanDecisionSpace = function(x.heatmap, marker.heatmap, x.sinks = NULL, marker.sinks = NULL, frame="x3", scene="scene") {
  if (frame == "x1") {
    frame = ~x1
  } else if (frame == "x2") {
    frame = ~x2
  } else if (frame == "x3") {
    frame = ~x3
  }
  
  p = plot_ly(
    scene = scene
  ) %>% add_markers(
    data = x.heatmap,
    type = "scatter3d",
    x = ~x1, y = ~x2, z = ~x3,
    frame = frame,
    mode = "markers",
    marker = marker.heatmap
  )
  
  if (!is.null(x.sinks)) {
    p = p %>% add_markers(
      data = x.sinks,
      type = "scatter3d",
      mode = "markers",
      x = ~x1, y = ~x2, z = ~x3,
      marker = marker.sinks
    )
  }
  
  p
}

