function_families <- c(
  "Bi-Objective BBOB" = "biobj_bbob",
  "DTLZ Functions" = "dtlz",
  "MMF Functions" = "mmf",
  "MPM2 Generator" = "mpm2",
  "ZDT Functions" = "zdt",
  "Other" = "other"
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
