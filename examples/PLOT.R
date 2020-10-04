library(moPLOT)

# PLOT as used in the 2020 PPSN paper

# Define a smoof function fn
fn = smoof::makeDTLZ1Function(2,2)

# Generate a design in the (rectangular) decision space of fn
design = generateDesign(fn, 500**2)
design$obj.space = calculateObjectiveValues(design$dec.space, fn, parallelize = T)

# Calculate single-objective and multi-objective gradients
gradients = computeGradientFieldGrid(design)

# Calculate divergence of MOG
divergence = computeDivergenceGrid(gradients$multi.objective, design$dims, design$step.sizes)

# Calculate locally efficient points
less = localEfficientSetSkeleton(design, gradients, divergence, integration="fast")

# Plot the results
ggplotPLOT(design$dec.space, design$obj.space, less$sinks, less$height)
ggplotPLOTObjSpace(design$obj.space, less$sinks, less$height)

