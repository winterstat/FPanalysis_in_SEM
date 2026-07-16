#------------------------------------------------
# Author: Tive B.P. Khumalo and Sonja D. Winter
# Date: 07/09/2026
#------------------------------------------------


#-------------------------------------------------
# Fitting propensity: Potential Outcomes
#-------------------------------------------------

# loading packages
pkgs <- c('ockhamSEM', 'parallel', 'doParallel', 'lavaan', 'tidyverse')
lapply(pkgs, library, character.only = T)

# Set up data matrix
p <- 24
temp_mat <- diag(p)
colnames(temp_mat) <- rownames(temp_mat) <- paste0("x", seq(1:p))

# model specification
mod_2A <- "
f1 =~ x1 + x2 + x3 + x4 +x5 + x6 + x7 + x8 +
x9 + x10 + x11 + x12 +x13 + x14 + x15 + x16 + 
x17 + x18 + x19 + x20 + x21 + x22 + x23 + x24
"

mod_2B <- "
f1 =~ x1 + x2 + x3 + x4 +x5 + x6 + x7 + x8
f2 =~ x9 + x10 + x11 + x12 +x13 + x14 + x15 + x16
f3 =~ x17 + x18 + x19 + x20 + x21 + x22 + x23 + x24
"

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
model_results_2A <- vector(mode = "list", length = length(max_iter))
names(model_results_2A) <- paste0("max_", max_iter)

model_results_2B <- vector(mode = "list", length = length(max_iter))
names(model_results_2B) <- paste0("max_", max_iter)

model_results_2C <- vector(mode = "list", length = length(max_iter))
names(model_results_2C) <- paste0("max_", max_iter)

# number of correlation matrices
R <- 1e4

# parallel simulation
ncores <- 64

cl <- makeCluster(ncores, type = "FORK")

registerDoParallel(cl)

on.exit(stopCluster(cl))

# Wrapping everything using lapply:
for(i in 1:length(max_iter)){

# Model Fitting  
mod_2A_fit <- cfa(mod_2A, sample.cov = temp_mat,
control = list(iter.max = max_iter[i]),
optim.force.converged = TRUE, bounds = "none",
sample.nobs = 1e3)  
  
mod_2B_fit <- cfa(mod_2B, sample.cov = temp_mat,
control = list(iter.max = max_iter[i]),
optim.force.converged = TRUE, bounds = "none",
sample.nobs = 1e3)
    
mod_2C_fit <- cfa(mod_2C, sample.cov = temp_mat,
control = list(iter.max = max_iter[i]),
optim.force.converged = TRUE, bounds = "none",
sample.nobs = 1e3)

# Fitting Propensity
model_results_2A[[i]] <- run.fitprop(mod_2A_fit,
fit.measure = c("logl", "cfi", "srmr", "rmsea"),
rmethod = "onion", reps = R, onlypos = FALSE, 
cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)


model_results_2B[[i]] <- run.fitprop(mod_2B_fit,
fit.measure = c("logl", "cfi", "srmr", "rmsea"),
rmethod = "onion", reps = R, onlypos = FALSE,
cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)

model_results_2C[[i]] <- run.fitprop(mod_2C_fit,
fit.measure = c("logl", "cfi", "srmr", "rmsea"),
rmethod = "onion", reps = R, onlypos = FALSE, 
cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)

}

# saving the results under conditions reference-indicator and no variance bounds
saveRDS(model_results_2A, "model_2A_fp_con_1.rds")
saveRDS(model_results_2B, "model_2B_fp_con_1.rds")
saveRDS(model_results_2C, "model_2C_fp_con_1.rds")

#-------------------------------------------------------------------------------
# Results under conditions reference-indicator and positive variance

# using the a container (list)
model_results_22A <- vector(mode = "list", length = length(max_iter))
names(model_results_22A) <- paste0("max_", max_iter)

model_results_22B <- vector(mode = "list", length = length(max_iter))
names(model_results_22B) <- paste0("max_", max_iter)

model_results_22C <- vector(mode = "list", length = length(max_iter))
names(model_results_22C) <- paste0("max_", max_iter)



# Wrapping everything using lapply:
for(i in 1:length(max_iter)){
  
  # Model Fitting  
  mod_22A_fit <- cfa(mod_2A, sample.cov = temp_mat,
                    control = list(iter.max = max_iter[i]),
                    optim.force.converged = TRUE, bounds = "pos.var",
                    sample.nobs = 1e3)  
  
  mod_22B_fit <- cfa(mod_2B, sample.cov = temp_mat,
                    control = list(iter.max = max_iter[i]),
                    optim.force.converged = TRUE, bounds = "pos.var",
                    sample.nobs = 1e3)
  
  mod_22C_fit <- cfa(mod_2C, sample.cov = temp_mat,
                    control = list(iter.max = max_iter[i]),
                    optim.force.converged = TRUE, bounds = "pos.var",
                    sample.nobs = 1e3)
  
  # Fitting Propensity
  model_results_22A[[i]] <- run.fitprop(mod_22A_fit,
  fit.measure = c("logl", "cfi", "srmr", "rmsea"),
  rmethod = "onion", reps = R, onlypos = FALSE, 
  cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
  
  model_results_22B[[i]] <- run.fitprop(mod_22B_fit,
  fit.measure = c("logl", "cfi", "srmr", "rmsea"),
  rmethod = "onion", reps = R, onlypos = FALSE,
  cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
  model_results_22C[[i]] <- run.fitprop(mod_22C_fit,
  fit.measure = c("logl", "cfi", "srmr", "rmsea"),
  rmethod = "onion", reps = R, onlypos = FALSE, 
  cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
}


# saving the results under conditions reference-indicator and positive variance
saveRDS(model_results_22A, "model_2A_fp_con_2.rds")
saveRDS(model_results_22B, "model_2B_fp_con_2.rds")
saveRDS(model_results_22C, "model_2C_fp_con_2.rds")
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Results under conditions unit variance and no variance bounds

# using the a container (list)
model_results_222A <- vector(mode = "list", length = length(max_iter))
names(model_results_222A) <- paste0("max_", max_iter)

model_results_222B <- vector(mode = "list", length = length(max_iter))
names(model_results_222B) <- paste0("max_", max_iter)

model_results_222C <- vector(mode = "list", length = length(max_iter))
names(model_results_222C) <- paste0("max_", max_iter)



# Wrapping everything using lapply:
for(i in 1:length(max_iter)){
  
  # Model Fitting  
  mod_222A_fit <- cfa(mod_2A, sample.cov = temp_mat,
                     control = list(iter.max = max_iter[i]),
                     optim.force.converged = TRUE, bounds = "none",
                     sample.nobs = 1e3, std.lv = TRUE)  
  
  mod_222B_fit <- cfa(mod_2B, sample.cov = temp_mat,
                     control = list(iter.max = max_iter[i]),
                     optim.force.converged = TRUE, bounds = "none",
                     sample.nobs = 1e3, std.lv = TRUE)
  
  mod_222C_fit <- cfa(mod_2C, sample.cov = temp_mat,
                     control = list(iter.max = max_iter[i]),
                     optim.force.converged = TRUE, bounds = "none",
                     sample.nobs = 1e3, std.lv = TRUE)
  
  # Fitting Propensity
  model_results_222A[[i]] <- run.fitprop(mod_222A_fit,
  fit.measure = c("logl", "cfi", "srmr", "rmsea"),
  rmethod = "onion", reps = R, onlypos = FALSE, 
  cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
  
  model_results_222B[[i]] <- run.fitprop(mod_222B_fit,
  fit.measure = c("logl", "cfi", "srmr", "rmsea"),
  rmethod = "onion", reps = R, onlypos = FALSE,
  cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
  model_results_222C[[i]] <- run.fitprop(mod_222C_fit,
  fit.measure = c("logl", "cfi", "srmr", "rmsea"),
  rmethod = "onion", reps = R, onlypos = FALSE, 
  cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
}


# saving the results under conditions unit variance and no variance bounds
saveRDS(model_results_222A, "model_2A_fp_con_3.rds")
saveRDS(model_results_222B, "model_2B_fp_con_3.rds")
saveRDS(model_results_222C, "model_2C_fp_con_3.rds")
#-------------------------------------------------------------------------------

# Results under conditions unit variance and positive variance

# using the a container (list)
model_results_2222A <- vector(mode = "list", length = length(max_iter))
names(model_results_2222A) <- paste0("max_", max_iter)

model_results_2222B <- vector(mode = "list", length = length(max_iter))
names(model_results_2222B) <- paste0("max_", max_iter)

model_results_2222C <- vector(mode = "list", length = length(max_iter))
names(model_results_2222C) <- paste0("max_", max_iter)



# Wrapping everything using lapply:
for(i in 1:length(max_iter)){
  
  # Model Fitting  
  mod_2222A_fit <- cfa(mod_2A, sample.cov = temp_mat,
                      control = list(iter.max = max_iter[i]),
                      optim.force.converged = TRUE, bounds = "pos.var",
                      sample.nobs = 1e3, std.lv = TRUE)  
  
  mod_2222B_fit <- cfa(mod_2B, sample.cov = temp_mat,
                      control = list(iter.max = max_iter[i]),
                      optim.force.converged = TRUE, bounds = "pos.var",
                      sample.nobs = 1e3, std.lv = TRUE)
  
  mod_2222C_fit <- cfa(mod_2C, sample.cov = temp_mat,
                      control = list(iter.max = max_iter[i]),
                      optim.force.converged = TRUE, bounds = "pos.var",
                      sample.nobs = 1e3, std.lv = TRUE)
  
  # Fitting Propensity
  model_results_2222A[[i]] <- run.fitprop(mod_2222A_fit,
  fit.measure = c("logl", "cfi", "srmr", "rmsea"),
  rmethod = "onion", reps = R, onlypos = FALSE, 
  cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
  
  model_results_2222B[[i]] <- run.fitprop(mod_2222B_fit,
  fit.measure = c("logl", "cfi", "srmr", "rmsea"),
  rmethod = "onion", reps = R, onlypos = FALSE,
  cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
  model_results_2222C[[i]] <- run.fitprop(mod_2222C_fit,
  fit.measure = c("logl", "cfi", "srmr", "rmsea"),
  rmethod = "onion", reps = R, onlypos = FALSE, 
  cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
}


# saving the results under conditions unit variance and positive variance
saveRDS(model_results_2222A, "model_2A_fp_con_4.rds")
saveRDS(model_results_2222B, "model_2B_fp_con_4.rds")
saveRDS(model_results_2222C, "model_2C_fp_con_4.rds")