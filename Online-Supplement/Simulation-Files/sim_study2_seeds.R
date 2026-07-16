#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
#-------------------------------------------------------------------------------

# OckhamSEM simulation study 2: All correlation matrices ----
# 24 variables with a varying number of factors
# Total data space: 1 million correlation matrices

# Set a seed
set.seed(20250107)

# Generate 1000 seeds
seeds <- round(runif(n = 1000, min = 1, max = 1e+6), 0)
models <- 1:48

pars <- expand.grid(models, seeds)
colnames(pars) <- c("model","seed")

# Save seeds
saveRDS(pars, file = "src/sim_study2_seeds.RDS")
