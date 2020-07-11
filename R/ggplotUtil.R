as.minimalistic.image = function(g) {
  # Create a minimalist version of the given plot
  g = g +
    theme(legend.position = "none",
          axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank()
    )
  
  return(g)
}

variable.as.expression = function(v) {
  # Convert x1-x3, y1-y3 to expressions
  # Otherwise return variable name
  if (v == "x1") {
    expression(x[1])
  } else if (v == "x2") {
    expression(x[2])
  } else if (v == "x3") {
    expression(x[3])
  } else if (v == "y1") {
    expression(y[1])
  } else if (v == "y2") {
    expression(y[2])
  } else if (v == "y3") {
    expression(y[3])
  } else {
    v
  }
}
