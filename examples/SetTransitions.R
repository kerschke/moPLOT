fn <- makeExtendedBiObjBBOBFunction(2, 10, 1)

# Generate a design in the (rectangular) decision space of fn
design = generateDesign(fn, 500**2)
design$obj.space = calculateObjectiveValues(design$dec.space, fn, parallelize = T)

# Calculate single-objective and multi-objective gradients
gradients = computeGradientFieldGrid(design, prec.angle = 0, normalized.scale = FALSE)

# Calculate divergence of MOG
divergence = computeDivergenceGrid(gradients$multi.objective, design$dims, design$step.sizes)

# Calculate locally efficient points
less = localEfficientSetSkeleton(design, gradients, divergence, integration="fast", with.basins = TRUE, use.integration.sinks = FALSE)

basins_df <- cbind.data.frame(
  design$dec.space,
  basins = less$basins,
  is_sink = 1:nrow(design$dec.space) %in% less$sinks
)

basins_df$basins[basins_df$basins == -1] <- NA

ggplot(data = basins_df, aes(x1, x2)) +
  geom_raster(aes(fill = factor(basins))) +
  geom_raster(fill = ifelse(basins_df$is_sink, "black", "transparent"))



basins <- less$basins

n_descent_target <- table(basins[basins != -1])
transition_sinks <- intersect(which(less$set_transitions != -1), less$sinks)

set_transitions <- t(sapply(transition_sinks, function(s) {
  c(from = basins[s], to = less$set_transitions[s])
}))
set_transitions <- set_transitions[!duplicated(set_transitions),,drop=FALSE]

sets <- lapply(sort(unique(basins[basins != -1])), function(b) {
  set_ids <- intersect(which(basins == b), less$sinks)
  list(
    dec_space = design$dec.space[set_ids,],
    obj_space = design$obj.space[set_ids,]
  )
})

# === Set Visualization ===

compute_reach_proportions <- function(graph, weights) {
  sapply(igraph::V(graph), function(v) {
    reachable_v <- igraph::subcomponent(graph, v, "in")
    sum(weights[reachable_v]) / sum(weights)
  })
}

compute_nondominated_sets <- function(sets) {
  obj_space <- Reduce(rbind, lapply(sets, function(set) set$obj_space))
  
  # When do we change between different sets of the trace?
  set_change <- cumsum(sapply(sets, function(set) nrow(set$obj_space)))
  
  nd <- obj_space %>% t %>% ecr::nondominated()
  
  sapply(seq_along(set_change), function(i) {
    if (i == 1) {
      lower = 1
    } else {
      lower = set_change[i - 1] + 1
    }
    
    sum(nd[lower:set_change[i]])
  })
}

library(tidygraph)
library(ggraph)

if (length(set_transitions) != 0) {
  tbl_transitions <- as_tbl_graph(set_transitions)
  prop <- compute_reach_proportions(tbl_transitions, n_descent_target)
} else {
  tbl_transitions <- NULL
  prop <- n_descent_target / sum(n_descent_target)
}

set_nd_counts <- compute_nondominated_sets(sets)
node_color <- ifelse(set_nd_counts > 0, "green", "red")

(1 / prop[set_nd_counts > 0]) %>% sort(decreasing = TRUE)

ggraph(tbl_transitions, layout = "stress") +
  geom_edge_fan(arrow = arrow(length = unit(4, "mm")),
                end_cap = circle(4, "mm")) +
  geom_node_point(aes(size = prop), color = node_color) +
  scale_radius(limits = c(0,1))
