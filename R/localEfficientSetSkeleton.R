#' @export
localEfficientSetSkeleton = function(design, gradients, divergence, integration = "fast", with.basins = FALSE) {
  less = list()
  
  cat('Finding critical points ...\n')
  
  lnd = locallyNondominatedCPP(design$obj.space, design$dims, TRUE)
  critical = getCriticalPointsCellCPP(gradients$multi.objective, gradients$single.objective, divergence, lnd, design$dims, sinks_only = TRUE)

  sinks = critical$sinks

  cat('Integrating vector field ...\n')
  if (integration == "fast") {
    integrated = computeCumulatedPathLengths(design$dec.space, gradients$multi.objective, sinks, fix.diagonals = T, prec.vector.length = 0, prec.norm = 0)
  } else {
    integrated = integrateVectorField(gradients$multi.objective, design$dims, sinks)
    dim(integrated$height) = c(length(integrated$height), 1)
    colnames(integrated$height) = c("height")
  }

  less$height = integrated$height
  # TODO  how can this even happen?!
  less$height[is.na(less$height)] = 0

  if (with.basins) {
    # integration.sinks = sort(unique(integrated$last.visited))
    integration.sinks = sinks

    ccs = connectedComponentsGrid(integration.sinks, design$dims)
    valid.ccs.ids = (ccs != 0 & !(ccs %in% as.numeric(names(table(ccs)))[table(ccs) < 4]))
    valid.ccs = ccs[valid.ccs.ids]
    valid.sinks = integration.sinks[valid.ccs.ids]

    less$sinks = valid.sinks

    cat('Calculating basins of attraction ...\n')
    sink.to.basin = rep(-1, prod(design$dims))
    sink.to.basin[valid.sinks] = valid.ccs

    basins = sapply(1:length(integrated$last.visited), function(i) {
      ifelse(integrated$last.visited[i] != -1, sink.to.basin[integrated$last.visited[i]], -1)
    })
    
    # Update Basin ID to be (-1,) 1, 2, 3, ...
    
    old_basins <- sort(unique(basins[basins != -1]))
    new_basins <- c(1:(length(old_basins)))
    
    basin_map <- rep(0, max(basins))
    basin_map[old_basins] <- new_basins
    
    basins <- sapply(basins, function(b) {
      ifelse(b == -1, -1, basin_map[b])
    })

    less$basins = basins

    basins_list = changeOfBasin(basins, design$dims, sinks)
    
    ridges = union(basins_list$ridges, which(basins == -1))
    less$ridges = ridges
    
    less$set_transitions = basins_list$set_transitions
    
  } else {
    less$sinks = sinks
  }

  return(less)
}
