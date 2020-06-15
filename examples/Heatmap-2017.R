library(moPLOT)

# Heatmaps as defined in 2017 EMO paper

# Define a smoof function fn
fn = smoof::makeDTLZ1Function(2,2)

# Generate a design in the (rectangular) decision space of fn
design = generateDesign(fn, 500**2)
design$obj.space = calculateObjectiveValues(design$dec.space, fn, parallelize = T)

# Calculate single-objective and multi-objective gradients
gradients = computeGradientField(design$dec.space, fn, impute.boundary = F, parallelize = T)

paths = computeCumulatedPathLengths(design$dec.space, gradients$multi.objective)

# Plot the results
ggplotHeatmap(cbind.data.frame(design$dec.space, height=paths$height))
ggplotObjectiveSpace(cbind.data.frame(design$obj.space, height=paths$height))
