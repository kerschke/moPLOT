#' Local Multi-Objective PCP Visualizations
#'
#' Returns a **Local PCP** (**P**arallel **C**oordinate **P**lot) visualization
#' as a [ggplot2::ggplot] object for the given visualization data. Local PCPs
#' highlight the different locally efficient sets.
#'
#' @param design [[`list`]]\cr
#' @param less [[`list`]]\cr
#' @param space `"objective" || "decision" || "both"`
#' @param scale [[`character`]]\cr
#' @param alphaLines [[`numeric`]]\cr
#' @template arg_checkdata
#'
#' @export
ggplotLocalPCP <- function(design, less, space = c("objective", "decision", "both"),
                           scale = "uniminmax", alphaLines = 0.1, check.data = TRUE) {
  if (length(space) > 1) {
    space <- space[[1]]
  }
  
  if (check.data) {
    assertList(less)
    assertList(design)
    assertNumeric(alphaLines, lower = 0, upper = 1)
  }

  objectives <- ncol(design$obj.space)
  dimensions <- ncol(design$dec.space)
  
  combine_df <- switch(space,
                       "objective" = design$obj.space,
                       "decision" = design$dec.space,
                       "both" = cbind.data.frame(design$dec.space, design$obj.space)
  )
  
  combine_df <- cbind.data.frame(combine_df, basin = as.factor(less$basins))
  combine_df <- combine_df[less$sinks,]
  combine_df <- combine_df[order(combine_df$basin), ]

  ggparcoord(
    combine_df,
    columns = 1L:(ncol(combine_df) - 1L),
    groupColumn = "basin",
    alphaLines = alphaLines,
    scale = scale
  ) +
    labs(
      x = element_blank(),
      y = "Normalized Value",
      color = "Basin"
    ) + scale_x_discrete(
      labels = c(
        "x1" = expression(x[1]),
        "x2" = expression(x[2]),
        "x3" = expression(x[3]),
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
#' @param design [[`list`]]\cr
#' @param space `"objective" || "decision" || "both"`
#' @param scale [[`character`]]\cr
#' @param alphaLines [[`numeric`]]\cr
#' @param nondominated_ids `NULL ||`[[`numeric`]]\cr
#' @template arg_checkdata
#'
#' @export
ggplotGlobalPCP <- function(design, space = c("objective", "decision", "both"),
                            scale = "uniminmax", alphaLines = 0.2,
                            nondominated_ids = NULL, check.data = TRUE) {
  if (length(space) > 1) {
    space <- space[[1]]
  }

  if (check.data) {
    assertList(design)
    assertNumeric(alphaLines, lower = 0, upper = 1)
    assertIntegerish(nondominated_ids, lower = 1, null.ok = TRUE)
    assertChoice(space, c("objective", "decision", "both"))
  }

  dimensions <- ncol(design$dec.space)
  objectives <- ncol(design$obj.space)

  if (is.null(nondominated_ids)) {
    nondominated_ids <- nondominated(design$obj.space, design$dims)
  }

  combine_df <- switch(space,
    "objective" = design$obj.space,
    "decision" = design$dec.space,
    "both" = cbind.data.frame(design$dec.space, design$obj.space)
  )

  combine_df <- combine_df[nondominated_ids, ]

  ggparcoord(
    combine_df,
    alphaLines = alphaLines,
    scale = scale
  ) +
    labs(
      x = element_blank(),
      y = "Normalized Value"
    ) + scale_x_discrete(
      labels = c(
        "x1" = expression(x[1]),
        "x2" = expression(x[2]),
        "x3" = expression(x[3]),
        "y1" = expression(y[1]),
        "y2" = expression(y[2]),
        "y3" = expression(y[3])
      ),
      expand = expansion(mult = 0.1)
    ) +
    theme_minimal()
}
