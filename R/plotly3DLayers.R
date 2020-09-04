#' @export
plotly3DLayers = function(grid, fn, mode = "decision.space", no.steps = 20, impute.zero = T) {
  # grid: list of obj.space, dims, dec.space, step.sizes
  # fn: smoof function, 3 dimensional decision space

  n = smoof::getNumberOfObjectives(fn)
  lower = smoof::getLowerBoxConstraints(fn)
  upper = smoof::getUpperBoxConstraints(fn)

  if (impute.zero) {
    grid$height = imputeZero(grid$height)
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
  # print(head(x.boundaries))

  decision.scene = list(
    aspectmode='cube',
    xaxis = list(range = c(lower[1],upper[1]), title='x₁'),
    yaxis = list(range = c(lower[2],upper[2]), title='x₂'),
    zaxis = list(range = c(lower[3],upper[3]), title='x₃')
  )

  if (n == 3) {
    objective.scene = list(
      aspectmode='cube',
      xaxis = list(range = c(min(x.boundaries$y1),max(x.boundaries$y1)), title='y₁'),
      yaxis = list(range = c(min(x.boundaries$y2),max(x.boundaries$y2)), title='y₂'),
      zaxis = list(range = c(min(x.boundaries$y3),max(x.boundaries$y3)), title='y₃')
    )
  }

  marker = plotlyMarker(grid$height)

  if (mode == "both") {
    x.shared = highlight_key(x.boundaries)
    p.decision = plotly3DLayersDecisionSpace(x.shared, fn, marker, scene="scene")
    p.objective = plotly3DLayersObjectiveSpace(x.shared, fn, marker, scene="scene2")

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
    plotly3DLayersDecisionSpace(x.boundaries, fn, marker) %>% layout(
      scene = decision.scene
    )
  } else if (mode == "objective.space") {
    if (n == 3) {
      plotly3DLayersObjectiveSpace(x.boundaries, fn, marker) %>% layout(
        scene = objective.scene
      )
    } else {
      plotly3DLayersObjectiveSpace(x.boundaries, fn, marker)
    }
  }

}

plotly3DLayersObjectiveSpace = function(x, fn, marker.style, scene="scene") {
  p = smoof::getNumberOfObjectives(fn)

  if (p == 2) {
    plot_ly(data = x,
            type="scattergl",
            x=~y1,y=~y2,
            frame=~frame,
            ids=~paste(x1,x2,x3),
            mode = "markers",
            marker=marker.style
    ) %>% animation_opts(
      frame = 1000,
      transition = 0
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
  } else if (p == 3) {
    plot_ly(data = x,
            type="scatter3d",
            x=~y1,y=~y2,z=~y3,
            frame=~frame,
            scene=scene,
            ids=~paste(x1,x2,x3),
            mode = "markers",
            marker=marker.style
    ) %>% animation_opts(
      frame = 1000,
      transition = 0
    )
  }
}

plotly3DLayersDecisionSpace = function(x, fn, marker.style, scene="scene") {
  plot_ly(data = x,
          type="scatter3d",
          x=~x1,y=~x2,z=~x3,
          scene=scene,
          frame = ~frame,
          ids=~paste(x1,x2,x3),
          mode = "markers",
          marker=marker.style
  ) %>% animation_opts(
    frame = 1000,
    transition = 0
  )
}
