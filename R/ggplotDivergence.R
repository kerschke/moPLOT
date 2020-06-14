#' @export
ggplotDivergence = function(grid) {
  positive.divergence = as.data.frame(grid$dec.space[which(grid$div > 0),])

  ggplotHeatmap(
    cbind.data.frame(grid$dec.space, grid$height),
    minimalistic.image = T
  ) + geom_tile(aes(x = x1, y = x2),
                data = positive.divergence
  )
}

#' @export
ggplotSecondOrder = function(grid, with.heatmap = F) {
  divergence = cbind.data.frame(grid$dec.space, div=grid$div)

  if (with.heatmap) {
    ggplotHeatmap(
      cbind.data.frame(grid$dec.space, grid$height),
      minimalistic.image = T
    ) + geom_tile(aes(x = x1, y = x2, fill = div, alpha = ifelse(is.na(div), 0, 1), width = grid$step.sizes[1], height=grid$step.sizes[2]),
                  data = divergence
    )
  } else {
    mins = apply(grid$dec.space, 2, min)
    maxs = apply(grid$dec.space, 2, max)
    ggplot(divergence) +
      xlim(mins[1], maxs[1]) +
      ylim(mins[2], maxs[2]) +
      geom_tile(aes(x = x1, y = x2, fill=div, alpha=ifelse(is.na(div), 0, 1), width = grid$step.sizes[1], height=grid$step.sizes[2])) +
      theme_minimal() +
      theme(legend.position = "none",
            axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank()
      )
  }
}

#' @export
ggplotHesse = function(grid) {
  locallyNondominatedCPP()
  divergence = cbind.data.frame(grid$dec.space[grid$div != 0,], div=grid$div[grid$div != 0])
}
