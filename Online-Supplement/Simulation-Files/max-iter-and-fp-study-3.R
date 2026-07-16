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
p <- 9
temp_mat_a <- diag(p)
colnames(temp_mat_a) <- rownames(temp_mat_a) <- paste0("x", seq(1:p))

p <- 27
temp_mat_b <- diag(p)
colnames(temp_mat_b) <- rownames(temp_mat_b) <- paste0("x", seq(1:p))

p <- 45
temp_mat_c <- diag(p)
colnames(temp_mat_c) <- rownames(temp_mat_c) <- paste0("x", seq(1:p))


# model specification
mod_3A <- "
f1 =~ x1 + x2 + x3
f2 =~ x4 + x5 + x6
f3 =~ x7 + x8 + x9
"

mod_3B <- "
f1 =~ x1 + x2 + x3 + x4 +x5 + x6 + x7 + x8 + x9
f2 =~ x10 + x11 + x12 +x13 + x14 + x15 + x16 + x17 + x18
f3 =~ x19 + x20 + x21 + x22 + x23 + x24 + x25 + x26 + x27
"

mod_3C <- " 
f1 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13 + x14 + x15 
f2 =~ x16 +  x17 + x18 + x19 + x20 + x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + x30
f3 =~ x31 + x32 + x33 + x34 + x35 + x36 + x37 + x38 + x39 + x40 + x41 + x42 + x43 + x44 + x45"

# Conditions: complete data space, 10,000 correlation matrices, forced convergence
# reference-indicator, no variance bounds and max.iter = (150, 250, 750,
# 1500, 2500, 4000, 6000, 10,000)

# maximum number of iterations
max_iter <- c(150, 250, 750, 1500, 2500, 4000, 6000, 10000)

# using the a container (list)
model_results_3A <- vector(mode = "list", length = length(max_iter))
names(model_results_3A) <- paste0("max_", max_iter)

model_results_3B <- vector(mode = "list", length = length(max_iter))
names(model_results_3B) <- paste0("max_", max_iter)

model_results_3C <- vector(mode = "list", length = length(max_iter))
names(model_results_3C) <- paste0("max_", max_iter)

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
  mod_3A_fit <- cfa(mod_3A, sample.cov = temp_mat_a,
                    control = list(iter.max = max_iter[i]),
                    optim.force.converged = TRUE, bounds = "none",
                    sample.nobs = 1e3)  
  
  mod_3B_fit <- cfa(mod_3B, sample.cov = temp_mat_b,
                    control = list(iter.max = max_iter[i]),
                    optim.force.converged = TRUE, bounds = "none",
                    sample.nobs = 1e3)
  
  mod_3C_fit <- cfa(mod_3C, sample.cov = temp_mat_c,
                    control = list(iter.max = max_iter[i]),
                    optim.force.converged = TRUE, bounds = "none",
                    sample.nobs = 1e3)
  
  # Fitting Propensity
  model_results_3A[[i]] <- run.fitprop(mod_3A_fit,
                                       fit.measure = c("logl", "cfi", "srmr", "rmsea"),
                                       rmethod = "onion", reps = R, onlypos = FALSE, 
                                       cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
  
  model_results_3B[[i]] <- run.fitprop(mod_3B_fit,
                                       fit.measure = c("logl", "cfi", "srmr", "rmsea"),
                                       rmethod = "onion", reps = R, onlypos = FALSE,
                                       cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
  model_results_3C[[i]] <- run.fitprop(mod_3C_fit,
                                       fit.measure = c("logl", "cfi", "srmr", "rmsea"),
                                       rmethod = "onion", reps = R, onlypos = FALSE, 
                                       cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
}

# saving the results under conditions reference-indicator and no variance bounds
saveRDS(model_results_3A, "model_3A_fp_con_1.rds")
saveRDS(model_results_3B, "model_3B_fp_con_1.rds")
saveRDS(model_results_3C, "model_3C_fp_con_1.rds")

#-------------------------------------------------------------------------------
# Results under conditions reference-indicator and positive variance

# using the a container (list)
model_results_33A <- vector(mode = "list", length = length(max_iter))
names(model_results_33A) <- paste0("max_", max_iter)

model_results_33B <- vector(mode = "list", length = length(max_iter))
names(model_results_33B) <- paste0("max_", max_iter)

model_results_33C <- vector(mode = "list", length = length(max_iter))
names(model_results_33C) <- paste0("max_", max_iter)



# Wrapping everything using lapply:
for(i in 1:length(max_iter)){
  
  # Model Fitting  
  mod_33A_fit <- cfa(mod_3A, sample.cov = temp_mat_a,
                     control = list(iter.max = max_iter[i]),
                     optim.force.converged = TRUE, bounds = "pos.var",
                     sample.nobs = 1e3)  
  
  mod_33B_fit <- cfa(mod_3B, sample.cov = temp_mat_b,
                     control = list(iter.max = max_iter[i]),
                     optim.force.converged = TRUE, bounds = "pos.var",
                     sample.nobs = 1e3)
  
  mod_33C_fit <- cfa(mod_3C, sample.cov = temp_mat_c,
                     control = list(iter.max = max_iter[i]),
                     optim.force.converged = TRUE, bounds = "pos.var",
                     sample.nobs = 1e3)
  
  # Fitting Propensity
  model_results_33A[[i]] <- run.fitprop(mod_33A_fit,
                                        fit.measure = c("logl", "cfi", "srmr", "rmsea"),
                                        rmethod = "onion", reps = R, onlypos = FALSE, 
                                        cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
  
  model_results_33B[[i]] <- run.fitprop(mod_33B_fit,
                                        fit.measure = c("logl", "cfi", "srmr", "rmsea"),
                                        rmethod = "onion", reps = R, onlypos = FALSE,
                                        cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
  model_results_33C[[i]] <- run.fitprop(mod_33C_fit,
                                        fit.measure = c("logl", "cfi", "srmr", "rmsea"),
                                        rmethod = "onion", reps = R, onlypos = FALSE, 
                                        cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
}


# saving the results under conditions reference-indicator and positive variance
saveRDS(model_results_33A, "model_3A_fp_con_2.rds")
saveRDS(model_results_33B, "model_3B_fp_con_2.rds")
saveRDS(model_results_33C, "model_3C_fp_con_2.rds")
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Results under conditions unit variance and no variance bounds

# using the a container (list)
model_results_333A <- vector(mode = "list", length = length(max_iter))
names(model_results_333A) <- paste0("max_", max_iter)

model_results_333B <- vector(mode = "list", length = length(max_iter))
names(model_results_333B) <- paste0("max_", max_iter)

model_results_333C <- vector(mode = "list", length = length(max_iter))
names(model_results_333C) <- paste0("max_", max_iter)



# Wrapping everything using lapply:
for(i in 1:length(max_iter)){
  
  # Model Fitting  
  mod_333A_fit <- cfa(mod_2A, sample.cov = temp_mat_a,
                      control = list(iter.max = max_iter[i]),
                      optim.force.converged = TRUE, bounds = "none",
                      sample.nobs = 1e3, std.lv = TRUE)  
  
  mod_333B_fit <- cfa(mod_2B, sample.cov = temp_mat_b,
                      control = list(iter.max = max_iter[i]),
                      optim.force.converged = TRUE, bounds = "none",
                      sample.nobs = 1e3, std.lv = TRUE)
  
  mod_333C_fit <- cfa(mod_2C, sample.cov = temp_mat_c,
                      control = list(iter.max = max_iter[i]),
                      optim.force.converged = TRUE, bounds = "none",
                      sample.nobs = 1e3, std.lv = TRUE)
  
  # Fitting Propensity
  model_results_333A[[i]] <- run.fitprop(mod_333A_fit,
                                         fit.measure = c("logl", "cfi", "srmr", "rmsea"),
                                         rmethod = "onion", reps = R, onlypos = FALSE, 
                                         cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
  
  model_results_333B[[i]] <- run.fitprop(mod_333B_fit,
                                         fit.measure = c("logl", "cfi", "srmr", "rmsea"),
                                         rmethod = "onion", reps = R, onlypos = FALSE,
                                         cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
  model_results_333C[[i]] <- run.fitprop(mod_333C_fit,
                                         fit.measure = c("logl", "cfi", "srmr", "rmsea"),
                                         rmethod = "onion", reps = R, onlypos = FALSE, 
                                         cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
}


# saving the results under conditions unit variance and no variance bounds
saveRDS(model_results_333A, "model_3A_fp_con_3.rds")
saveRDS(model_results_333B, "model_3B_fp_con_3.rds")
saveRDS(model_results_333C, "model_3C_fp_con_3.rds")
#-------------------------------------------------------------------------------

# Results under conditions unit variance and positive variance

# using the a container (list)
model_results_3333A <- vector(mode = "list", length = length(max_iter))
names(model_results_3333A) <- paste0("max_", max_iter)

model_results_3333B <- vector(mode = "list", length = length(max_iter))
names(model_results_3333B) <- paste0("max_", max_iter)

model_results_3333C <- vector(mode = "list", length = length(max_iter))
names(model_results_3333C) <- paste0("max_", max_iter)



# Wrapping everything using lapply:
for(i in 1:length(max_iter)){
  
  # Model Fitting  
  mod_2222A_fit <- cfa(mod_2A, sample.cov = temp_mat_a,
                       control = list(iter.max = max_iter[i]),
                       optim.force.converged = TRUE, bounds = "pos.var",
                       sample.nobs = 1e3, std.lv = TRUE)  
  
  mod_2222B_fit <- cfa(mod_2B, sample.cov = temp_mat_b,
                       control = list(iter.max = max_iter[i]),
                       optim.force.converged = TRUE, bounds = "pos.var",
                       sample.nobs = 1e3, std.lv = TRUE)
  
  mod_2222C_fit <- cfa(mod_2C, sample.cov = temp_mat_c,
                       control = list(iter.max = max_iter[i]),
                       optim.force.converged = TRUE, bounds = "pos.var",
                       sample.nobs = 1e3, std.lv = TRUE)
  
  # Fitting Propensity
  model_results_3333A[[i]] <- run.fitprop(mod_3333A_fit,
                                          fit.measure = c("logl", "cfi", "srmr", "rmsea"),
                                          rmethod = "onion", reps = R, onlypos = FALSE, 
                                          cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
  
  model_results_3333B[[i]] <- run.fitprop(mod_3333B_fit,
                                          fit.measure = c("logl", "cfi", "srmr", "rmsea"),
                                          rmethod = "onion", reps = R, onlypos = FALSE,
                                          cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
  model_results_3333C[[i]] <- run.fitprop(mod_3333C_fit,
                                          fit.measure = c("logl", "cfi", "srmr", "rmsea"),
                                          rmethod = "onion", reps = R, onlypos = FALSE, 
                                          cluster = cl, seed = 3242, saveR = FALSE, saveModel = FALSE)
  
}


# saving the results under conditions unit variance and positive variance
saveRDS(model_results_3333A, "model_3A_fp_con_4.rds")
saveRDS(model_results_3333B, "model_3B_fp_con_4.rds")
saveRDS(model_results_3333C, "model_3C_fp_con_4.rds")