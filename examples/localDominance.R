library(moPLOT)
library(tidyverse)

fn <- smoof::makeBiObjBBOBFunction(dimensions = 2, fid = 10, iid = 10)
fn <- makeAsparFunction()
fn <- smoof::makeDTLZ1Function(2, 2)

design <- generateDesign(fn, 501**2)

ld_data <- computeLocalDominance(design$obj.space, design$dims)


basins <- sapply(ld_data$basins, function(v) {
  if (length(v) == 1) v
  else NA
})

chob <- changeOfBasin(basins, design$dims, ld_data$locally_efficient_ids)
basins[setdiff(chob$ridges, ld_data$locally_efficient_ids)] <- NA

display_height <- rep(0, length(ld_data$basins))
display_height[ld_data$locally_efficient_ids] <- -1
display_height[is.na(basins)] <- NA

ggplot(cbind.data.frame(design$dec.space, height = as.factor(display_height))) +
  coord_fixed() +
  geom_tile(aes(x1, x2, fill = height)) +
  scale_fill_manual(values = c("black", "gray"), na.value = "white") +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(cbind.data.frame(design$dec.space, height = as.factor(basins))) +
  coord_fixed() +
  geom_tile(aes(x1, x2, fill = height), alpha = (as.numeric(1:length(basins) %in% ld_data$locally_efficient_ids) * 0.5 + 0.5)) +
  scale_fill_viridis_d(na.value = NA) +
  theme_minimal() +
  theme(legend.position = "none")
