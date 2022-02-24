#' Set Transition Visualization in ggplot2
#'
#' @param design [[`list`]] \cr
#' @param less [[`list`]] \cr
#' @param node_size `"reachability" || "basin_size` \cr
#' @param layout `NULL || "first_dimensions" || "stress` \cr
#' @template arg_checkdata
#'
#' @return A [[ggplot2::ggplot]] object
#' @export
#'
#' @examples
#' 
ggplotSetTransitions <- function(design, less, node_size = c("reachability", "basin_size"),
                                 layout = NULL, check.data = TRUE) {
  if (!requireNamespace("tidygraph", quietly = TRUE) ||
      !requireNamespace("ggraph", quietly = TRUE)) {
    stop(
      "Packages \"tidygraph\" and \"ggraph\" must be installed to use ggplotSetTransitions().",
      call. = FALSE
    )
  }
  
  if (length(node_size) > 1L) {
    node_size <- node_size[[1L]]
  }
  
  if (check.data) {
    assertChoice(node_size, c("reachability", "basin_size"))
    assertChoice(layout, c("first_dimensions", "stress"), null.ok = TRUE)
    assertList(design)
    assertList(less)
  }
  
  dimensions <- length(design$dims)
  
  if (is.null(layout)) {
    layout = if (dimensions == 2L) "first_dimensions" else "stress"
  }
  
  # Prepare basins, sets, set transitions ====
  
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

  if (length(set_transitions) > 0) {
    tbl_transitions <- tidygraph::as_tbl_graph(igraph::graph(t(set_transitions), n = max(basins)))
  } else {
    tbl_transitions <- tidygraph::as_tbl_graph(igraph::graph(NULL, n = max(basins)))
  }
  
  if (length(set_transitions) != 0 && node_size == "reachability") {
    prop <- compute_reach_proportions(tbl_transitions, n_descent_target)
  } else {
    prop <- n_descent_target / sum(n_descent_target)
  }
  
  set_nd_counts <- compute_nondominated_sets(sets)
  node_color <- ifelse(set_nd_counts > 0, "green", "red")

  if (layout == "stress") {
    g <- ggraph::ggraph(tbl_transitions, layout = "stress") +
      ggraph::geom_node_point(aes(size = 1), color = "black", shape = 21) +
      ggraph::geom_node_point(aes(size = prop), color = node_color) +
      ggraph::geom_edge_fan(arrow = arrow(length = unit(4, "mm")),
                            end_cap = ggraph::circle(4, "mm")) +
      scale_size_area(limits = c(0,1)) +
      theme(legend.position = "none",
            panel.background = element_rect(fill = NA, size = 0),
            plot.background = element_rect(fill = NA, size = 0))
  } else {
    node_pos <- t(sapply(sets, function(set) {
      dec_space <- set$dec_space
      dec_space[ceiling(nrow(dec_space) / 2),]
    }))
    
    lower <- design$lower
    upper <- design$upper
    
    g <- ggraph::ggraph(tbl_transitions, layout = "manual", x = node_pos[,1], y = node_pos[,2]) +
      ggraph::geom_node_point(aes(size = 1), color = "black", shape = 21) +
      ggraph::geom_node_point(aes(size = prop), color = node_color) +
      ggraph::geom_edge_fan(arrow = arrow(length = unit(4, "mm")),
                            end_cap = ggraph::circle(4, "mm")) +
      scale_size_area(limits = c(0,1)) +
      theme_minimal() +
      coord_fixed(xlim = c(lower[1], upper[1]), ylim = c(lower[2], upper[2])) +
      theme(legend.position = "none",
            panel.background = element_rect(fill = NA, size = 0),
            plot.background = element_rect(fill = NA, size = 0)) +
      labs(x = expression(x[1]),
           y = expression(x[2]))
  }
  
  return(g)
}

# Helper Functions ====

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
  
  nd <- ecr::nondominated(t(obj_space))
  
  sapply(seq_along(set_change), function(i) {
    if (i == 1) {
      lower <- 1
    } else {
      lower <- set_change[i - 1] + 1
    }
    
    sum(nd[lower:set_change[i]])
  })
}
