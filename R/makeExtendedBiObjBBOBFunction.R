# FIDs of single-objective BBOB functions utilized in the
# regular bi-objective suite
fids <- c(1L, 2L, 6L, 8L, 13L, 14L, 15L, 17L, 20L, 21L)

# Regular suite: all (unique) combinations of aforementioned FIDs
fid_mapping <- expand.grid(fid_1 = fids, fid_2 = fids)
fid_mapping <- fid_mapping[fid_mapping[, 1L] <= fid_mapping[, 2L], ]
fid_mapping <- fid_mapping[order(fid_mapping[, 1L]),]

# Add all off-diagonal combinations per BBOB group (excluding FID 16)
# for the extended bi-objective suite
for (group in list(1L:5L, 6L:9L, 10L:14L, setdiff(15L:19L, 16L), 20L:24L)) {
  expanded <- expand.grid(fid_1 = group, fid_2 = group)
  expanded <- expanded[expanded[, 1L] < expanded[, 2L], ]
  expanded <- expanded[order(expanded[, 1L]),]
  
  fid_mapping <- rbind(fid_mapping, expanded)
  fid_mapping <- fid_mapping[!duplicated(fid_mapping),]
}

rownames(fid_mapping) <- 1:nrow(fid_mapping)

# Regularly, single objective IIDs are computed as
# IID_1 = 2 * IID + 1, and
# IID_2 = IID_1 + 1.
# 
# However, there are some exceptions:
iid_exceptions <- matrix(c(
  # IID, IID_1, IID_2
  1L, 2L, 4L,
  2L, 3L, 5L,
  9L, 19L, 21L,
  15L, 31L, 34L
), byrow = TRUE, ncol = 3L)

#' Extended Bi-Objective BBOB Functions with Metadata
#'
#' Smoof generator function for the extended bi-objective BBOB with validated
#' (corrected) instance mapping and returning useful metadata.
#'
#' @param dimensions `2L:40L`\cr
#'   Number of (decision space) dimensions
#' @param fid `1L:92L`\cr
#'   Function ID of the selected bi-objective BBOB problem.
#' @param iid `1L:10L`\cr
#'   Instance ID of the selected bi-objective BBOB problem.
#'
#' @export
makeExtendedBiObjBBOBFunction <- function(dimensions, fid, iid) {
  assert_int(dimensions, lower = 2L, upper = 40L)
  assert_int(fid, lower = 1L, upper = 92L)
  assert_int(iid, lower = 1L, upper = 15L)
  
  fid_1 <- fid_mapping[fid, 1]
  fid_2 <- fid_mapping[fid, 2]
  
  if (iid %in% iid_exceptions[,1]) {
    iid_1 <- iid_exceptions[iid == iid_exceptions[,1], 2]
    iid_2 <- iid_exceptions[iid == iid_exceptions[,1], 3]
  } else {
    iid_1 <- 2L * iid + 1
    iid_2 <- iid_1 + 1
  }
  
  fn_1 <- smoof::makeBBOBFunction(dimensions, fid_1, iid_1)
  fn_2 <- smoof::makeBBOBFunction(dimensions, fid_2, iid_2)
  
  param_set <- ParamHelpers::makeNumericParamSet("x", len = dimensions, lower = -5, upper = 5)
  
  smoof::makeMultiObjectiveFunction(
    name = sprintf("Extended Bi-Objective BBOB_%i_%i_%i", dimensions, fid, iid),
    id = paste0("ext_biobj_bbob_", fid, "_", iid, "_", dimensions),
    description = sprintf("%i-th noiseless Extended Bi-Objective BBOB function\n(FID: %i, IID: %i, DIMENSION: %i)",
                          fid, fid, iid, dimensions),
    fn = function(x) {
      drop(cbind(fn_1(x), fn_2(x)))
    },
    par.set = param_set,
    n.objectives = 2L,
    vectorized = TRUE
  )
}
