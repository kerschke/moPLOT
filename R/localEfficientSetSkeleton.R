#' @export
localEfficientSetSkeleton = function(design, gradients, divergence, integration="fast", with.basins = F) {
  less = list()

  cat('Finding critical points ...\n')

  # lnd = lapply(1:ncol(grid$obj.space), function(i) {
  #   locallyNondominatedCPP(matrix(grid$obj.space[,i]), grid$dims, T)
  # })
  
  lnd = locallyNondominatedCPP(design$obj.space, design$dims, TRUE)

  critical = getCriticalPointsCellCPP(gradients$multi.objective, gradients$single.objective, divergence, lnd, design$dims, FALSE)

  sinks = critical$sinks

  # less$so.min = list()
  # less$so.min[[1]] = intersect(critical.1$sinks, sinks)
  # less$so.min[[2]] = intersect(critical.2$sinks, sinks)

  cat('Integrating vector field ...\n')
  if (integration == "fast") {
    integrated = computeCumulatedPathLengths(design$dec.space, gradients$multi.objective, sinks, fix.diagonals = T, prec.vector.length = 0, prec.norm = 0)
  } else {
    integrated = integrateVectorField(gradients$multi.objective, design$dims, sinks)
  }

  less$height = integrated$height
  # TODO  how can this even happen?!
  less$height[is.na(less$height)] = 0

  if (with.basins) {
    integration.sinks = sort(unique(integrated$last.visited))

    ccs = connectedComponentsGrid(integration.sinks, grid$dims)
    valid.ccs.ids = (ccs != 0 & !(ccs %in% as.numeric(names(table(ccs)))[table(ccs) < 4]))
    valid.ccs = ccs[valid.ccs.ids]
    valid.sinks = integration.sinks[valid.ccs.ids]

    less$sinks = valid.sinks

    cat('3/3 Calculating basins of attraction ...\n')
    sink.to.basin = rep(-1, prod(grid$dims))
    sink.to.basin[valid.sinks] = valid.ccs

    basins = sapply(1:length(integrated$last.visited), function(i) {
      ifelse(integrated$last.visited[i] != -1, sink.to.basin[integrated$last.visited[i]], -1)
    })

    less$basins = basins

    ridges = changeOfBasin(basins, grid$dims)
    ridges = union(ridges, which(basins == -1))

    less$ridges = ridges
  } else {
    less$sinks = sinks
  }

  return(less)
}
