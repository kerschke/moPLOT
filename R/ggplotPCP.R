#' Local Multi-Objective PCP Visualizations
#'
#' Returns a **Local PCP** (**P**arallel **C**oordinate **P**lot) visualization
#' as a [ggplot2::ggplot] object for the given visualization data. Local PCPs
#' highlight the different locally efficient sets.
#'
#' @param obj.space [[`matrix`]]\cr
#' @param basins [[`numeric`]]\cr
#' @param scale [[`character`]]\cr
#' @param alphaLines [[`numeric`]]\cr
#' @template arg_checkdata
#'
#' @export
ggplotLocalPCP <- function(obj.space, basins, scale = "uniminmax",
                           alphaLines = 0.1, check.data = TRUE) {
  if (check.data) {
    assertIntegerish(basins)
    assertNumeric(obj.space)
    assertNumeric(alphaLines, lower = 0, upper = 1)

    assertMatrix(obj.space, min.cols = 2L, nrows = length(basins))
  }

  objectives <- ncol(obj.space)

  combine_df <- cbind.data.frame(obj.space, basin = as.factor(basins))
  combine_df <- combine_df[order(combine_df$basin), ]

  ggparcoord(
    combine_df,
    columns = 1L:objectives,
    groupColumn = (objectives + 1L),
    alphaLines = alphaLines,
    scale = scale
  ) +
    labs(
      x = element_blank(),
      y = "Normalized Value",
      color = "Basin"
    ) + scale_x_discrete(
      labels = c(
        "y1" = expression(y[1]),
        "y2" = expression(y[2]),
        "y3" = expression(y[3])
      ),
      expand = expansion(mult = 0.1)
    ) + scale_color_viridis_d() +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    theme_minimal()
}

#' Global Multi-Objective PCP Visualizations
#'
#' Returns a **Global PCP** (**P**arallel **C**oordinate **P**lot) visualization
#' as a [ggplot2::ggplot] object for the given visualization data. Global PCPs
#' only show the Pareto front, i.e., the globally efficient points.
#'
#' @param obj.space [[`matrix`]]\cr
#' @param scale [[`character`]]\cr
#' @param alphaLines [[`numeric`]]\cr
#' @template arg_checkdata
#'
#' @export
ggplotGlobalPCP <- function(obj.space, scale = "uniminmax",
                            alphaLines = 0.2, check.data = TRUE) {
  if (check.data) {
    assertNumeric(obj.space)
    assertNumeric(alphaLines, lower = 0, upper = 1)
  }

  objectives <- ncol(obj.space)

  ggparcoord(
    obj.space,
    columns = 1L:objectives,
    alphaLines = alphaLines,
    scale = scale
  ) +
    labs(
      x = element_blank(),
      y = "Normalized Value"
    ) + scale_x_discrete(
      labels = c(
        "y1" = expression(y[1]),
        "y2" = expression(y[2]),
        "y3" = expression(y[3])
      ),
      expand = expansion(mult = 0.1)
    ) +
    theme_minimal()
}
