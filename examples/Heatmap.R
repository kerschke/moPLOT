library(moPLOT)

# Gradient field heatmap improved with explicitly calculated locally efficient points
# and imputation of the MOG at the boundary

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
ggplotHeatmap(cbind.data.frame(design$dec.space, height=less$height))
ggplotObjectiveSpace(cbind.data.frame(design$obj.space, height=less$height))
