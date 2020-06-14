#' @export
plotly3DScan = function(grid, fn, mode = "decision.space", impute.zero = T) {
  # grid: list of obj.space, dims, dec.space, step.sizes
  # fn: smoof function, 3 dimensional decision space
  
  n = smoof::getNumberOfObjectives(fn)
  lower = smoof::getLowerBoxConstraints(fn)
  upper = smoof::getUpperBoxConstraints(fn)
  
  if (impute.zero) {
    grid$height = imputeZero(grid$height)
  }
  
  decision.scene = list(
    aspectmode='cube',
    xaxis = list(range = c(lower[1],upper[1]), title='x₁'),
    yaxis = list(range = c(lower[2],upper[2]), title='x₂'),
    zaxis = list(range = c(lower[3],upper[3]), title='x₃')
  )
  
  x = cbind.data.frame(grid$dec.space, grid$height, grid$obj.space)
  
  if (n == 3) {
    objective.scene = list(
      aspectmode='cube',
      xaxis = list(range = c(min(x$y1),max(x$y1)), title='y₁'),
      yaxis = list(range = c(min(x$y2),max(x$y2)), title='y₂'),
      zaxis = list(range = c(min(x$y3),max(x$y3)), title='y₃')
    )
  }
  
  marker = plotlyMarker(grid)
  
  if (mode == "both") {
    x.shared = highlight_key(x)
    p.decision = plotly3DScanDecisionSpace(x.shared, fn, marker, scene="scene") %>% layout(
      xaxis = list(range = c(min(x$y1),max(x$y1)), title='y₁'),
      yaxis = list(range = c(min(x$y2),max(x$y2)), title='y₂')
    )
    p.objective = plotly3DScanObjectiveSpace(x.shared, fn, marker, scene="scene2")
    
    domain.left = list(
      x=c(0,0.5),
      y=c(0,1)
    )
    decision.scene$domain = domain.left
    
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
      title = paste("Decision and Objective Space of", smoof::getName(fn)),
      scene = decision.scene,
      scene2 = objective.scene
    ) %>% highlight(
      on="plotly_click",
      off="plotly_deselect",
      opacityDim = 0.5,
      color = "red"
    ) %>% hide_guides() %>% animation_opts(
      # does not work?
      frame = 5000
    )
  } else if (mode == "decision.space") {
    plotly3DScanDecisionSpace(x,fn, marker) %>% layout(
      scene = decision.scene
    )
  } else if (mode == "objective.space") {
    if (n == 3) {
      plotly3DScanObjectiveSpace(x,fn, marker) %>% layout(
        scene = objective.scene
      )
    } else {
      plotly3DScanObjectiveSpace(x,fn, marker) %>% layout(
        xaxis = list(range = c(min(x$y1),max(x$y1)), title='y₁'),
        yaxis = list(range = c(min(x$y2),max(x$y2)), title='y₂')
      )
    }
  }
  
}

plotly3DScanObjectiveSpace = function(x, fn, marker.style, frame="x3", scene="scene") {
  n = smoof::getNumberOfObjectives(fn)
  
  if (frame == "x1") {
    frame = ~x1
    ids = ~paste(x2,x3)
  } else if (frame == "x2") {
    frame = ~x2
    ids = ~paste(x1,x3)
  } else if (frame == "x3") {
    frame = ~x3
    ids = ~paste(x1,x2)
  }
  
  if (n == 2) {
    plot_ly(data = x,
            type="scatter",
            x=~y1,y=~y2,
            ids=ids,
            frame=frame,
            mode = "markers",
            marker=marker.style
    ) %>% animation_opts(
      frame = 1000,
      transition = 0
    )
  } else if (n == 3) {
    plot_ly(data = x,
            type="scatter3d",
            x=~y1,y=~y2,z=~y3,
            ids = ids,
            frame = frame,
            scene = scene,
            mode = "markers",
            marker=marker.style
    ) %>% animation_opts(
      frame = 1000,
      transition = 0
    )
  }
}

plotly3DScanDecisionSpace = function(x, fn, marker.style, frame="x3", scene="scene") {
  if (frame == "x1") {
    frame = ~x1
    ids = ~paste(x2,x3)
  } else if (frame == "x2") {
    frame = ~x2
    ids = ~paste(x1,x3)
  } else if (frame == "x3") {
    frame = ~x3
    ids = ~paste(x1,x2)
  }
  
  plot_ly(data = x,
          type = "scatter3d",
          x=~x1,y=~x2,z=~x3,
          frame = frame,
          ids = ids,
          scene = scene,
          mode = "markers",
          marker=marker.style
  ) %>% animation_opts(
    frame = 1000
  )
}

