#' @export
plotly3DPareto = function(grid, fn, mode = "decision.space", impute.zero = T) {
  # grid: list of obj.space, dims, dec.space, step.sizes
  # fn: smoof function, 3 dimensional decision space
  
  n = smoof::getNumberOfObjectives(fn)
  lower = smoof::getLowerBoxConstraints(fn)
  upper = smoof::getUpperBoxConstraints(fn)
  
  if (impute.zero) {
    grid$height = imputeZero(grid$height)
  }
  
  nondomIndices = nondominated(grid$obj.space, grid$dims)
  
  x.nondom = cbind(
    grid$dec.space[nondomIndices,,drop=F],
    grid$height[nondomIndices,,drop=F],
    grid$obj.space[nondomIndices,,drop=F]
  )
  x.nondom = as.data.frame(x.nondom)
  
  decision.scene = list(
    aspectmode='cube',
    xaxis = list(range = c(lower[1],upper[1]), title='x<sub>1</sub>'),
    yaxis = list(range = c(lower[2],upper[2]), title='x<sub>2</sub>'),
    zaxis = list(range = c(lower[3],upper[3]), title='x<sub>3</sub>')
  )
  
  if (n == 3) {
    objective.scene = list(
      aspectmode='cube',
      xaxis = list(range = c(min(x.nondom[,'y1']),max(x.nondom[,'y1'])), title='y<sub>1</sub>'),
      yaxis = list(range = c(min(x.nondom[,'y2']),max(x.nondom[,'y2'])), title='y<sub>2</sub>'),
      zaxis = list(range = c(min(x.nondom[,'y3']),max(x.nondom[,'y3'])), title='y<sub>3</sub>')
    )
  }
  
  marker = plotlyMarker(grid$height)
  
  if (mode == "both") {
    x.shared = highlight_key(x.nondom)
    p.decision = plotly3DParetoDecisionSpace(x.shared, fn, marker, scene="scene")
    p.objective = plotly3DParetoObjectiveSpace(x.shared, fn, marker, scene="scene2")
    
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
    ) %>% hide_guides()
  } else if (mode == "decision.space") {
    plotly3DParetoDecisionSpace(x.nondom, fn, marker) %>% layout(
      scene = decision.scene
    )
  } else if (mode == "objective.space") {
    if (n == 3) {
      plotly3DParetoObjectiveSpace(x.nondom, fn, marker) %>% layout(
        scene = objective.scene
      )
    } else {
      plotly3DParetoObjectiveSpace(x.nondom, fn, marker)
    }
  }
  
}

#' @export
nondominated = function(obj.space, dims) {
  locallyNondom = locallyNondominatedCPP(obj.space, dims, TRUE)
  nondom = ecr::nondominated(t(obj.space[locallyNondom,]))
  
  locallyNondom[nondom]
}

plotly3DParetoObjectiveSpace = function(x, fn, marker.style, scene="scene") {
  n = smoof::getNumberOfObjectives(fn)
  
  if (n == 2) {
    plot_ly(data = x,
            type="scattergl",
            x=~y1,y=~y2,
            mode = "markers",
            marker = marker.style
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
  } else if (n == 3) {
    plot_ly(data = x,
            type="scatter3d",
            x=~y1,y=~y2,z=~y3,
            scene = scene,
            mode = "markers",
            marker = marker.style
    )
  }
}

plotly3DParetoDecisionSpace = function(x, fn, marker.style, scene="scene") {
  plot_ly(data = x,
          type="scatter3d",
          x=~x1,y=~x2,z=~x3,
          scene = scene,
          mode = "markers",
          marker = marker.style
  )
}
