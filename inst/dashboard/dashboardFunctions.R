function_families <- c(
  "Bi-Objective BBOB" = "biobj_bbob",
  "DTLZ Functions" = "dtlz",
  "MMF Functions" = "mmf",
  "MOP Functions" = "mop",
  "MPM2 Generator (TODO)" = "mpm2",
  "ZDT Functions" = "zdt",
  "Other" = "other"
)

biobj_bbob_functions <- list(
  "Bi-Objective BBOB" = smoof::makeBiObjBBOBFunction
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

# Aspar functions

f1_1 = function(x) (x[1]**4 - 2*x[1]**2 + x[2]**2 + 1)
f2_1 = function(x) ((x[1] + 0.5)**2 + (x[2]-2)**2)
f3_1 = function(x) ((x[1] + 0.25) ** 4 + 3 * (x[2] - 1) ** 2)
f_2d2d = function(x) c(f1_1(x), f2_1(x))
f_2d3d = function(x) c(f1_1(x), f2_1(x), f3_1(x))

test.2d.2d = smoof::makeMultiObjectiveFunction(name = "2D->2D Test Function", id = "test_2d2d", description = "", fn = f_2d2d,
                                               par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-2,-1), upper = c(2,3)))

test.2d.3d = smoof::makeMultiObjectiveFunction(name = "2D->3D Test Function", id = "test_2d3d", description = "", fn = f_2d3d,
                                               par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-2,-1), upper = c(2,3)))

f1_2 = function(x) (x[1]**4 - 2*x[1]**2 + x[2]**2 + 1 + x[3] ** 2)
f2_2 = function(x) ((x[1] + 0.5)**2 + (x[2]-2)**2 + (x[3] - 1) ** 4)
f3_2 = function(x) ((x[1] + 0.25) ** 4 + 3 * (x[2] - 1) ** 2 + (x[3] - 1) ** 2)
f_3d2d = function(x) c(f1_2(x), f2_2(x))
f_3d3d = function(x) c(f1_2(x), f2_2(x), f3_2(x))

test.3d.2d = smoof::makeMultiObjectiveFunction(name = "3D->2D Test Function", id = "test_3d2d", description = "", fn = f_3d2d,
                                               par.set = ParamHelpers::makeNumericParamSet(len = 3, lower = c(-2,-1,-2), upper = c(2,3,2)))

test.3d.3d = smoof::makeMultiObjectiveFunction(name = "3D->3D Test Function", id = "test_3d3d", description = "", fn = f_3d3d,
                                               par.set = ParamHelpers::makeNumericParamSet(len = 3, lower = c(-2,-1,-2), upper = c(2,3,2)))

other_functions = list(
  "Viennet Function" = smoof::makeViennetFunction,
  "Kursawe" = smoof::makeKursaweFunction # TODO Dimensions
  # "2D->2D Test Function" = test.2d.2d, # TODO convert to generator function
  # "2D->3D Test Function" = test.2d.3d,
  # "3D->2D Test Function" = test.3d.2d,
  # "3D->3D Test Function" = test.3d.3d
)
