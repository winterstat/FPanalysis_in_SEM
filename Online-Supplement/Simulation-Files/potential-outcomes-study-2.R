#------------------------------------------------
# Author: Tive B.P. Khumalo and Sonja D. Winter
# Date: 07/09/2026
#------------------------------------------------


#-------------------------------------------------------------------
# Fitting propensity: Relationship between maximum iterations and FP
#-------------------------------------------------------------------


# loading packages
pkgs <- c('ockhamSEM', 'parallel', 'doParallel', 'lavaan', 'tidyverse')
lapply(pkgs, library, character.only = T)

# Set up data matrix
p <- 24
temp_mat <- diag(p)
colnames(temp_mat) <- rownames(temp_mat) <- paste0("x", seq(1:p))

# model specification
mod_2C <- " 
f1 =~ x1 + x2 + x3 + x4
f2 =~ x5 + x6 + x7 + x8
f3 =~ x9 + x10 + x11 + x12
f4 =~ x13 + x14 + x15 + x16
f5 =~ x17 + x18 + x19 + x20
f6 =~ x21 + x22 + x23 + x24"

# Conditions: complete data space, 10,000 correlation matrices, forced convergence
# reference-indicator, no variance bounds and max.iter = (150, 250, 750,
# 1500, 2500, 4000, 6000, 10,000)

# maximum number of iterations
max_iter <- c(150, 250, 750, 1500, 2500, 4000, 6000, 10000)

# using the a container (list)
model_results <- vector(mode = "list", length = length(max_iter))
names(model_results) <- paste0("max_", max_iter)

# number of correlation matrices
R <- 1e4

# parallel simulation
# cl <- makeCluster(8, type = "FORK")
cl <- makeCluster(8)
registerDoParallel(cl)

# on.exit(stopCluster(cl))

# Wrapping everything using lapply:
for(i in 1:length(max_iter)){

# Model Fitting  
mod_2C_fit <- cfa(mod_2C, sample.cov = temp_mat,
control = list(iter.max = max_iter[i]),
optim.force.converged = TRUE, bounds = "none",
sample.nobs = 1e3)

# Fitting Propensity
model_results[[i]] <- run.fitprop(mod_2C_fit,
fit.measure = c("logl", "cfi", "srmr", "rmsea"),
rmethod = "onion", reps = R, onlypos = FALSE, 
cluster = cl, seed = 3242, saveR = F, saveModel = F)

}

stopCluster(cl)

# saving the results
saveRDS(model_results, "model_2C_fp.rds")





