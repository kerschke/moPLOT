#' Rudolph 9-to-9
#' 
#' Creates the Rudolph 9-to-9 test function
#'
#' @return Multi-objective smoof function
#' @export
#'
#' @examples
#' 
makeRudolph9to9Function = function() {
  rudolph9to9 <- function(x) {
    a = 0.5
    b = 2
    c = 1
    
    t1 = sign(x[1]) * min(ceiling((abs(x[1]) - a - (c / 2)) / (2 * a + c)), 1)
    t2 = sign(x[2]) * min(ceiling((abs(x[2]) - b / 2) / b), 1)
    
    f1 = (x[1] - t1 * (c + 2*a) + a)^2 + (x[2] - t2 * b)^2
    f2 = (x[1] - t1 * (c + 2*a) - a)^2 + (x[2] - t2 * b)^2
    
    penalty = 0
    
    if (t1 == -1) penalty <- penalty + 0.03
    if (t1 ==  0) penalty <- penalty + 0.00
    if (t1 ==  1) penalty <- penalty + 0.06
    
    if (t2 == -1) penalty <- penalty + 0.09
    if (t2 ==  0) penalty <- penalty + 0.00
    if (t2 ==  1) penalty <- penalty + 0.18

    l <- 1
    
    f1 <- f1 + penalty * l
    f2 <- f2 + penalty * l
    
    c(f1, f2)
  }
  
  smoof::makeMultiObjectiveFunction(
    name = "Rudolph 9-to-9 Function", id = "rudolph9to9", description = "", fn = rudolph9to9,
    par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-3, -3), upper = c(3, 3)))
}

