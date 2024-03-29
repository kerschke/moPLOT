#' @export
plotly3DLayers = function(grid, fn, sinks = NULL, mode = "decision.space", no.steps = 20L, impute.zero = TRUE,
                          colorscale.sinks = plotlyColorscale(), colorscale.heatmap = plotlyColorscale(gray.colorscale)) {
  # grid: list of obj.space, dims, dec.space, step.sizes
  # fn: smoof function, 3 dimensional decision space

  n = smoof::getNumberOfObjectives(fn)
  lower = smoof::getLowerBoxConstraints(fn)
  upper = smoof::getUpperBoxConstraints(fn)

  if (impute.zero) {
    grid$height = imputeZero(grid$height)
  }
  
  if (!is.null(sinks)) {
    x.sinks = cbind.data.frame(grid$dec.space[sinks,], grid$obj.space[sinks,])
    dom.counter = ecr::doNondominatedSorting(t(grid$obj.space[sinks,]))$dom.counter
  } else {
    x.sinks = NULL
    dom.counter = NULL
  }

  maxh = calculateMaxDisplayHeightCPP(grid$height, grid$dims, F)
  min.height = min(maxh)
  max.height = max(maxh[grid$height <= min.height])

  ecdf.height = stats::ecdf(grid$height)
  q.min = ecdf.height(min.height)
  q.max = ecdf.height(max.height)

  x.boundaries = c()

  quantiles = seq(q.min, q.max, length.out = no.steps)
  heights = stats::quantile(grid$height, probs=quantiles)

  for (height in heights) {ids = which(grid$height <= height & maxh >= height)
    boundary = as.data.frame(cbind(
      grid$dec.space[ids,,drop=F],
      grid$height[ids,,drop=F],
      grid$obj.space[ids,,drop=F]
    ))

    boundary$frame = height
    x.boundaries = rbind(x.boundaries, boundary)
  }

  x.boundaries = as.data.frame(x.boundaries)

  decision.scene = list(
    aspectmode='cube',
    xaxis = list(range = c(lower[1], upper[1]), title='x<sub>1</sub>'),
    yaxis = list(range = c(lower[2], upper[2]), title='x<sub>2</sub>'),
    zaxis = list(range = c(lower[3], upper[3]), title='x<sub>3</sub>')
  )

  if (n == 3) {
    objective.scene = list(
      aspectmode='cube',
      xaxis = list(range = c(min(x.boundaries$y1),max(x.boundaries$y1)), title='y<sub>1</sub>'),
      yaxis = list(range = c(min(x.boundaries$y2),max(x.boundaries$y2)), title='y<sub>2</sub>'),
      zaxis = list(range = c(min(x.boundaries$y3),max(x.boundaries$y3)), title='y<sub>3</sub>')
    )
  } else {
    objective.scene = list()
  }

  if (mode == "both") {
    p.decision = plotly3DLayersDecisionSpace(x.boundaries, fn, x.sinks, dom.counter, scene="scene",
                                             colorscale.sinks = colorscale.sinks, colorscale.heatmap = colorscale.heatmap)
    p.objective = plotly3DLayersObjectiveSpace(x.boundaries, fn, x.sinks, dom.counter, scene="scene2",
                                               colorscale.sinks = colorscale.sinks, colorscale.heatmap = colorscale.heatmap)

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
    ) %>% hide_guides()
  } else if (mode == "decision.space") {
    p = plotly3DLayersDecisionSpace(x.boundaries, fn, x.sinks, dom.counter,
                                    colorscale.sinks = colorscale.sinks, colorscale.heatmap = colorscale.heatmap) %>% layout(
      scene = decision.scene
    )
  } else if (mode == "objective.space") {
    if (n == 3) {
      p = plotly3DLayersObjectiveSpace(x.boundaries, fn, x.sinks, dom.counter,
                                       colorscale.sinks = colorscale.sinks, colorscale.heatmap = colorscale.heatmap) %>% layout(
        scene = objective.scene
      )
    } else {
      p = plotly3DLayersObjectiveSpace(x.boundaries, fn, x.sinks, dom.counter,
                                       colorscale.sinks = colorscale.sinks, colorscale.heatmap = colorscale.heatmap)
    }
  }

  p %>% animation_opts(
    frame = 1000,
    transition = 0
  ) %>% hide_guides()
}

plotly3DLayersObjectiveSpace = function(x, fn, x.sinks = NULL, dom.counter = NULL, scene = "scene",
                                        colorscale.sinks, colorscale.heatmap) {
  n.obj = smoof::getNumberOfObjectives(fn)
  
  if (!is.null(x.sinks)) {
    marker.sinks = plotlyMarker(dom.counter + 1, colorscale = colorscale.sinks)
    marker.heatmap = plotlyMarker(x$height, colorscale = colorscale.heatmap)
  } else {
    marker.heatmap = plotlyMarker(x$height, colorscale = colorscale.heatmap)
  }

  if (n.obj == 2) {
    p = plot_ly() %>% add_markers(
      data = x,
      type = "scattergl",
      x = ~y1, y = ~y2,
      frame = ~frame,
      mode = "markers",
      marker = marker.heatmap
    ) %>% layout(
      xaxis = list(
        title = "y<sub>1</sub>",
        constrain = "domain"
      ),
      yaxis = list(
        title = "y<sub>2</sub>",
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
  } else if (n.obj == 3) {
    p = plot_ly(
      scene = scene
    ) %>% add_markers(
      data = x,
      type = "scatter3d",
      x=~y1, y=~y2, z=~y3,
      frame = ~frame,
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

plotly3DLayersDecisionSpace = function(x, fn, x.sinks = NULL, dom.counter = NULL, scene="scene",
                                       colorscale.sinks, colorscale.heatmap) {
  if (!is.null(x.sinks)) {
    marker.sinks = plotlyMarker(dom.counter + 1, colorscale = colorscale.sinks)
    marker.heatmap = plotlyMarker(x$height, colorscale = colorscale.heatmap)
  } else {
    marker.heatmap = plotlyMarker(x$height, colorscale = colorscale.sinks)
  }
  
  p = plot_ly(
    scene = scene
  ) %>% add_markers(
    data = x,
    type = "scatter3d",
    x = ~x1, y = ~x2, z = ~x3,
    frame = ~frame,
    mode = "markers",
    marker = marker.heatmap
  )
  
  if (!is.null(x.sinks)) {
    p = p %>% add_markers(
      data = x.sinks,
      type = "scatter3d",
      x = ~x1, y = ~x2, z = ~x3,
      mode = "markers",
      marker = marker.sinks
    )
  }
  
  p
}
