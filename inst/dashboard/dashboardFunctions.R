benchmark_sets <- c(
  "Bi-objective BBOB" = "biobj_bbob",
  "DTLZ Functions" = "dtlz",
  "MinDist Functions" = "mindist",
  "MMF Functions" = "mmf",
  "MOP Functions" = "mop",
  "MPM2 Generator" = "mpm2",
  "ZDT Functions" = "zdt",
  "Other" = "other"
)

biobj_bbob_functions <- list(
  "Bi-objective BBOB" = smoof::makeBiObjBBOBFunction
)

dtlz_functions <- list(
  "DTLZ1" = smoof::makeDTLZ1Function,
  "DTLZ2" = smoof::makeDTLZ2Function,
  "DTLZ3" = smoof::makeDTLZ3Function,
  "DTLZ4" = smoof::makeDTLZ4Function,
  "DTLZ5" = smoof::makeDTLZ5Function,
  "DTLZ6" = smoof::makeDTLZ6Function,
  "DTLZ7" = smoof::makeDTLZ7Function
)

mmf_functions <- list(
  "MMF1" = smoof::makeMMF1Function,
  "MMF1e" = smoof::makeMMF1eFunction,
  "MMF1z" = smoof::makeMMF1zFunction,
  "MMF2" = smoof::makeMMF2Function,
  "MMF3" = smoof::makeMMF3Function,
  "MMF4" = smoof::makeMMF4Function,
  "MMF5" = smoof::makeMMF5Function,
  "MMF6" = smoof::makeMMF6Function,
  "MMF7" = smoof::makeMMF7Function,
  "MMF8" = smoof::makeMMF8Function,
  "MMF9" = smoof::makeMMF9Function,
  "MMF10" = smoof::makeMMF10Function,
  "MMF11" = smoof::makeMMF11Function,
  "MMF12" = smoof::makeMMF12Function,
  "MMF13" = smoof::makeMMF13Function,
  "MMF14" = smoof::makeMMF14Function,
  "MMF14a" = smoof::makeMMF14aFunction,
  "MMF15" = smoof::makeMMF15Function,
  "MMF15a" = smoof::makeMMF15aFunction,
  "SYMPART-simple" = smoof::makeSYMPARTsimpleFunction,
  "SYMPART-rotated" = smoof::makeSYMPARTrotatedFunction,
  "Omni test" = smoof::makeOmniTestFunction
)

zdt_functions = list(
  "ZDT1" = smoof::makeZDT1Function,
  "ZDT2" = smoof::makeZDT2Function,
  "ZDT3" = smoof::makeZDT3Function,
  "ZDT4" = smoof::makeZDT4Function,
  "ZDT6" = smoof::makeZDT6Function
)

mop_functions = list(
  "MOP1" = smoof::makeMOP1Function,
  "MOP2" = smoof::makeMOP2Function,
  "MOP3" = smoof::makeMOP3Function,
  "MOP4" = smoof::makeMOP4Function,
  "MOP5" = smoof::makeMOP5Function,
  "MOP6" = smoof::makeMOP6Function,
  "MOP7" = smoof::makeMOP7Function
)

makeMinDistFunction = function(centers.f1 = list(c(-2, -1), c(2, 1)),
                               centers.f2 = list(c(-2, 1), c(2, -1)),
                               centers.f3 = list(c(0, 1), c(0, -1))) {
  if (is.null(centers.f3)) {
    f <- function(x) {
      y1 <- min(sapply(centers.f1, function(center) sqrt(sum((x - center) ** 2))))
      y2 <- min(sapply(centers.f2, function(center) sqrt(sum((x - center) ** 2))))
      
      c(y1, y2)
    }
    
    all.centers <- Reduce(rbind, c(centers.f1, centers.f2))
  } else {
    f <- function(x) {
      y1 <- min(sapply(centers.f1, function(center) sqrt(sum((x - center) ** 2))))
      y2 <- min(sapply(centers.f2, function(center) sqrt(sum((x - center) ** 2))))
      y3 <- min(sapply(centers.f3, function(center) sqrt(sum((x - center) ** 2))))
      
      c(y1, y2, y3)
    }
    
    all.centers <- Reduce(rbind, c(centers.f1, centers.f2, centers.f3))
  }
  
  lower <- apply(all.centers, 2, min) - 1
  upper <- apply(all.centers, 2, max) + 1
  
  smoof::makeMultiObjectiveFunction(name = "MinDist Function", id = "mindist", description = "", fn = f,
                                    par.set = ParamHelpers::makeNumericParamSet(len = ncol(all.centers), lower = lower, upper = upper))
}

makeBiObjMinDistFunction = function(centers.f1 = list(c(-2, -1), c(2, 1)),
                                    centers.f2 = list(c(-2, 1), c(2, -1))) {
  makeMinDistFunction(centers.f1, centers.f2, centers.f3 = NULL)
}

mindist_functions = list(
  "Bi-objective MinDist" = makeBiObjMinDistFunction,
  "Tri-objective MinDist" = makeMinDistFunction
)

# MPM2 functions

makeBiObjMPM2Function = function(dimensions = 2, n.peaks.1 = 3, topology.1 = "random", seed.1 = 4,
                                 n.peaks.2 = 3, topology.2 = "random", seed.2 = 8) {
  f1 <- smoof::makeMPM2Function(n.peaks.1, dimensions, topology.1, seed.1)
  f2 <- smoof::makeMPM2Function(n.peaks.2, dimensions, topology.2, seed.2)
  
  smoof::makeGOMOPFunction(dimensions = dimensions, funs = list(f1, f2))
}

makeTriObjMPM2Function = function(dimensions = 2, n.peaks.1 = 3, topology.1 = "random", seed.1 = 4,
                                  n.peaks.2 = 3, topology.2 = "random", seed.2 = 8,
                                  n.peaks.3 = 3, topology.3 = "random", seed.3 = 12) {
  f1 <- smoof::makeMPM2Function(n.peaks.1, dimensions, topology.1, seed.1)
  f2 <- smoof::makeMPM2Function(n.peaks.2, dimensions, topology.2, seed.2)
  f3 <- smoof::makeMPM2Function(n.peaks.3, dimensions, topology.3, seed.3)
  
  smoof::makeGOMOPFunction(dimensions = dimensions, funs = list(f1, f2, f3))
}

mpm2_functions = list(
  "Bi-objective MPM2" = makeBiObjMPM2Function,
  "Tri-objective MPM2" = makeTriObjMPM2Function
)

# Aspar functions

f1_1 = function(x) (x[1]**4 - 2*x[1]**2 + x[2]**2 + 1)
f2_1 = function(x) ((x[1] + 0.5)**2 + (x[2]-2)**2)
f3_1 = function(x) ((x[1] + 0.25) ** 4 + 3 * (x[2] - 1) ** 2)
f_2d2d = function(x) c(f1_1(x), f2_1(x))
f_2d3d = function(x) c(f1_1(x), f2_1(x), f3_1(x))

f1_2 = function(x) (x[1]**4 - 2*x[1]**2 + x[2]**2 + 1 + x[3] ** 2)
f2_2 = function(x) ((x[1] + 0.5)**2 + (x[2]-2)**2 + (x[3] - 1) ** 4)
f3_2 = function(x) ((x[1] + 0.25) ** 4 + 3 * (x[2] - 1) ** 2 + (x[3] - 1) ** 2)
f_3d2d = function(x) c(f1_2(x), f2_2(x))
f_3d3d = function(x) c(f1_2(x), f2_2(x), f3_2(x))

makeAsparFunction <- function(dimensions = 2, n.objectives = 2) {
  if (dimensions == 2 && n.objectives == 2) {
    smoof::makeMultiObjectiveFunction(name = "2D->2D Test Function", id = "test_2d2d", description = "", fn = f_2d2d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-2,-1), upper = c(2,3)))
  } else if (dimensions == 2 && n.objectives == 3) {
    smoof::makeMultiObjectiveFunction(name = "2D->3D Test Function", id = "test_2d3d", description = "", fn = f_2d3d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-2,-1), upper = c(2,3)))
  } else if (dimensions == 3 && n.objectives == 2) {
    smoof::makeMultiObjectiveFunction(name = "3D->2D Test Function", id = "test_3d2d", description = "", fn = f_3d2d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 3, lower = c(-2,-1,-2), upper = c(2,3,2)))
  } else if (dimensions == 3 && n.objectives == 3) {
    smoof::makeMultiObjectiveFunction(name = "3D->3D Test Function", id = "test_3d3d", description = "", fn = f_3d3d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 3, lower = c(-2,-1,-2), upper = c(2,3,2)))
  }
}

other_functions = list(
  "Aspar" = makeAsparFunction,
  "BiSphere" = smoof::makeBiSphereFunction, # does not work as expected and is kinda boring
  "BK1" = smoof::makeBK1Function,
  "Dent" = smoof::makeDentFunction,
  "ED1" = smoof::makeED1Function,
  "ED2" = smoof::makeED2Function,
  "Kursawe" = smoof::makeKursaweFunction,
  "Viennet" = smoof::makeViennetFunction
)
