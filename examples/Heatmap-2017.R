library(moPLOT)

# Heatmaps as defined in 2017 EMO paper

# Define a smoof function fn
fn = smoof::makeDTLZ1Function(2,2)

# Generate a design in the (rectangular) decision space of fn
design = generateDesign(fn, 500**2)

# Calculate single-objective and multi-objective gradients
gradients = computeGradientField(design$dec.space, fn, impute.boundary = F, parallelize = T)

paths = computeCumulatedPathLengths(design$dec.space, gradients$multi.objective, dims = design$dims)

# Plot the results
ggplotHeatmap(cbind.data.frame(design$dec.space, height=paths$height))
ggplotObjectiveSpace(cbind.data.frame(design$obj.space, height=paths$height))
