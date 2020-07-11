#' @export
inferDesign = function(dec.space, obj.space) {
  upper = apply(dec.space, 2, max)
  lower = apply(dec.space, 2, min)
  p = ncol(dec.space)
  
  # Check input data is sensible
  
  assertTRUE(nrow(dec.space) == nrow(obj.space))
  
  # Data is expected to be ordered "in reverse"
  
  if (p == 2) {
    ordered.indices = order(dec.space[,2], dec.space[,1])
  } else if (p == 3) {
    ordered.indices = order(dec.space[,3], dec.space[,2], dec.space[,1])
  } else {
    stop(paste0("Cannot visualize data with ", p, " dimensions! Please use 2 or 3 dimensions only."))
  }
  
  # Order data by order of the decision space and convert to matrix (more efficient)
  
  dec.space = as.matrix(dec.space[ordered.indices,])
  obj.space = as.matrix(obj.space[ordered.indices,])
  
  # Determine grid points for each dimension
  
  grid.points = lapply(seq_col(dec.space), function(i) {
    unique(dec.space[,i])
  })
  
  # Assemble data to design object
  
  design = list()
  
  design$dec.space = dec.space
  design$obj.space = obj.space
  
  design$dims = sapply(grid.points, length)
  
  design$step.sizes = sapply(grid.points, function(grid.points.i) {
    steps = diff(sort(grid.points.i))
    
    # Check step sizes to four significant digits (relevant due to possible floating point errors)
    four_signif = ceiling(4 - log(min(steps), base=10))
    
    step.size = unique(round(steps, digits=four_signif))
    
    if (length(step.size) > 1) {
      stop(paste("Found multiple step sizes:", collapse(step.sizes)))
    } else {
      step.size
    }
  })
  
  # Check dimensions and number of observations match
  
  if (prod(design$dims) != nrow(design$dec.space)) {
    stop(
      paste0(
        "Extracted grid dimensions (",
        paste0(design$dims, collapse=", "),
        ") do not match number of observations (",
        paste0(nrow(design$dec.space)),
        ")!"
      )
    )
  }
  
  return(design)
}
