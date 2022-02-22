#' @export
plotly2DContours <- function(design, show.nondominated = TRUE) {
  z1 <- design$obj.space[,1L]
  dim(z1) <- design$dims
  z1 <- z1 %>% t
  
  z2 <- design$obj.space[,2L]
  dim(z2) <- design$dims
  z2 <- z2 %>% t
  
  p <- plot_ly() %>% 
    add_contour(x = unique(design$dec.space[,1]), y = unique(design$dec.space[,2]), z = z1,
                contours = list(
                  coloring = "none",
                  showlabels = TRUE),
                line = list(color = "blue")) %>% 
    add_contour(x = unique(design$dec.space[,1]), y = unique(design$dec.space[,2]), z = z2,
                contours = list(
                  coloring = "none",
                  showlabels = TRUE),
                line = list(color = "red"))
  
  if (ncol(design$obj.space) == 3L) {
    z3 <- design$obj.space[,3L]
    dim(z3) <- design$dims
    z3 <- z3 %>% t
    
    p <- p %>%
      add_contour(x = unique(design$dec.space[,1]), y = unique(design$dec.space[,2]), z = z3,
                  contours = list(
                    coloring = "none",
                    showlabels = TRUE),
                  line = list(color = "orange"))
    
  }
  
  if (show.nondominated) {
    nondom <- nondominated(design$obj.space, design$dims)
    # nondom <- which(ecr::nondominated(t(design$obj.space)))
    
    z.nondom <- rep(NA, nrow(design$dec.space))
    z.nondom[nondom] <- 1
    
    p <- p %>%
      add_heatmap(data = design$dec.space %>% as.data.frame, x = ~x1, y = ~x2,
                  z = z.nondom,
                  colorscale = plotlyColorscale(c("#000000", "#000000")),
                  showscale = FALSE)
      # add_markers(data = design$dec.space[nondom,] %>% as.data.frame, x = ~x1, y = ~x2,
      #             marker = list(color = "black"))
  }
  
  p %>%
    hide_legend() %>%
    layout(
      xaxis = list(
        title = "x<sub>1</sub>",
        constrain = "domain"
      ),
      yaxis = list(
        scaleanchor = "x",
        scaleratio = (diff(range(design$dec.space[,1]))) / diff(range(design$dec.space[,2])),
        title = "x<sub>2</sub>",
        constrain = "domain"
      )
    )
}
