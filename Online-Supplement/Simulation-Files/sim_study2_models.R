# OckhamSEM simulation study 2: All correlation matrices ----
# 24 variables with a varying number of factors
# Total data space: 1 million correlation matrices
# Conditions included: 
#           - space size (post simulation)
#           - Model (3 models)
#           - All correlations vs only positive (two different data spaces)
#           - ML iterations: 150 (default) or 1500 
#           - force convergence: F (default) or T
#           - factor identification: reference marker (default) or unit-variance
#           - variance bounds: "none" (default), or "pos.var"
#           - within poscor data space, we have: 3 (models) x 2 (MLit) x 2 (force conv) x 2 (ident) x 2 (bounds) = 48 conditions
# Fit indices exported
#           - -LL
#           - CFI
#           - SRMR
#           - RMSEA

# Load packages ----
library(ockhamSEM)

# Set up temp matrix ----
p <- 24
temp_mat <- diag(p)
colnames(temp_mat) <- rownames(temp_mat) <- paste0("x", seq(1:p))

# Set up models ----
mod1 <- '
f1 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + 
x9 + x10 + x11 + x12 + x13 + x14 + x15 + x16 +
x17 + x18 + x19 + x20 + x21 + x22 + x23 + x24
'

mod2 <- '
f1 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 
f2 =~ x9 + x10 + x11 + x12 + x13 + x14 + x15 + x16
f3 =~ x17 + x18 + x19 + x20 + x21 + x22 + x23 + x24
'

mod3 <- '
f1 =~ x1 + x2 + x3 + x4
f2 =~ x5 + x6 + x7 + x8
f3 =~ x9 + x10 + x11 + x12
f4 =~ x13 + x14 + x15 + x16
f5 =~ x17 + x18 + x19 + x20
f6 =~ x21 + x22 + x23 + x24'

# Set up different conditions ----

# A: force convergence = FALSE, ML iterations = 150, identification = reference marker, bounds = "none"
mod1a.fit <- cfa(mod1, sample.cov = temp_mat,     # sample covariance/correlation matrix
                 control = list(iter.max = 150),  # increases the number of maximum iterations
                 optim.force.converged = F,       # forces the model to converge if True else 
                 bounds = "none",                 # it restricts or constrains the observed variances(positive or negative or both)(none is the default case)
                 sample.nobs = 1000)              # number of observations in the correlation/covariance matrix               
# the same comments apply on the models below.

mod2a.fit <- cfa(mod2, sample.cov = temp_mat,    
                 control = list(iter.max = 150),
                 optim.force.converged = F,
                 bounds = "none",
                 sample.nobs = 1000)

mod3a.fit <- cfa(mod3, sample.cov = temp_mat, 
                 control = list(iter.max = 150),
                 optim.force.converged = F,
                 bounds = "none",
                 sample.nobs = 1000)

# B: force convergence = FALSE, ML iterations = 150, identification = reference marker, bounds = "pos.var"
mod1b.fit <- cfa(mod1, sample.cov = temp_mat, 
                 control = list(iter.max = 150),
                 optim.force.converged = F,
                 bounds = "pos.var",             # it restricts or constrains the observed variances to be positive
                 sample.nobs = 1000)

mod2b.fit <- cfa(mod2, sample.cov = temp_mat, 
                 control = list(iter.max = 150),
                 optim.force.converged = F,
                 bounds = "pos.var",
                 sample.nobs = 1000)

mod3b.fit <- cfa(mod3, sample.cov = temp_mat, 
                 control = list(iter.max = 150),
                 optim.force.converged = F,
                 bounds = "pos.var",
                 sample.nobs = 1000)

# C: force convergence = FALSE, ML iterations = 150, identification = unit variance, bounds = "none"
mod1c.fit <- cfa(mod1, sample.cov = temp_mat, 
                 std.lv = T,                     # identification = unit variance(the latent factor is now on a standard normal distribution)
                 control = list(iter.max = 150),
                 optim.force.converged = F,
                 bounds = "none",
                 sample.nobs = 1000)

mod2c.fit <- cfa(mod2, sample.cov = temp_mat, 
                 std.lv = T,
                 control = list(iter.max = 150),
                 optim.force.converged = F,
                 bounds = "none",
                 sample.nobs = 1000)

mod3c.fit <- cfa(mod3, sample.cov = temp_mat, 
                 std.lv = T,
                 control = list(iter.max = 150),
                 optim.force.converged = F,
                 bounds = "none",
                 sample.nobs = 1000)

# D: force convergence = FALSE, ML iterations = 150, identification = unit variance, bounds = "pos.var"
mod1d.fit <- cfa(mod1, sample.cov = temp_mat, 
                 std.lv = T,                    
                 control = list(iter.max = 150),
                 optim.force.converged = F,
                 bounds = "pos.var",
                 sample.nobs = 1000)

mod2d.fit <- cfa(mod2, sample.cov = temp_mat, 
                 std.lv = T,
                 control = list(iter.max = 150),
                 optim.force.converged = F,
                 bounds = "pos.var",
                 sample.nobs = 1000)

mod3d.fit <- cfa(mod3, sample.cov = temp_mat,
                 std.lv = T,
                 control = list(iter.max = 150),
                 optim.force.converged = F,
                 bounds = "pos.var",
                 sample.nobs = 1000)

# E: force convergence = FALSE, ML iterations = 1500, identification = reference marker, bounds = "none"
mod1e.fit <- cfa(mod1, sample.cov = temp_mat, 
                 control = list(iter.max = 1500),
                 optim.force.converged = F,
                 bounds = "none",
                 sample.nobs = 1000)

mod2e.fit <- cfa(mod2, sample.cov = temp_mat, 
                 control = list(iter.max = 1500),
                 optim.force.converged = F,
                 bounds = "none",
                 sample.nobs = 1000)

mod3e.fit <- cfa(mod3, sample.cov = temp_mat, 
                 control = list(iter.max = 1500),
                 optim.force.converged = F,
                 bounds = "none",
                 sample.nobs = 1000)

# F: force convergence = FALSE, ML iterations = 1500, identification = reference marker, bounds = "pos.var"
mod1f.fit <- cfa(mod1, sample.cov = temp_mat, 
                 control = list(iter.max = 1500),
                 optim.force.converged = F,
                 bounds = "pos.var",
                 sample.nobs = 1000)

mod2f.fit <- cfa(mod2, sample.cov = temp_mat, 
                 control = list(iter.max = 1500),
                 optim.force.converged = F,
                 bounds = "pos.var",
                 sample.nobs = 1000)

mod3f.fit <- cfa(mod3, sample.cov = temp_mat, 
                 control = list(iter.max = 1500),
                 optim.force.converged = F,
                 bounds = "pos.var",
                 sample.nobs = 1000)

# G: force convergence = FALSE, ML iterations = 1500, identification = unit variance, bounds = "none"
mod1g.fit <- cfa(mod1, sample.cov = temp_mat, 
                 std.lv = T,
                 control = list(iter.max = 1500),
                 optim.force.converged = F,
                 bounds = "none",
                 sample.nobs = 1000)

mod2g.fit <- cfa(mod2, sample.cov = temp_mat, 
                 std.lv = T,
                 control = list(iter.max = 1500),
                 optim.force.converged = F,
                 bounds = "none",
                 sample.nobs = 1000)

mod3g.fit <- cfa(mod3, sample.cov = temp_mat, 
                 std.lv = T,
                 control = list(iter.max = 1500),
                 optim.force.converged = F,
                 bounds = "none",
                 sample.nobs = 1000)

# H: force convergence = FALSE, ML iterations = 1500, identification = unit variance, bounds = "pos.var"
mod1h.fit <- cfa(mod1, sample.cov = temp_mat, 
                 std.lv = T,
                 control = list(iter.max = 1500),
                 optim.force.converged = F,
                 bounds = "pos.var",
                 sample.nobs = 1000)

mod2h.fit <- cfa(mod2, sample.cov = temp_mat, 
                 std.lv = T,
                 control = list(iter.max = 1500),
                 optim.force.converged = F,
                 bounds = "pos.var",
                 sample.nobs = 1000)

mod3h.fit <- cfa(mod3, sample.cov = temp_mat,
                 std.lv = T,
                 control = list(iter.max = 1500),
                 optim.force.converged = F,
                 bounds = "pos.var",
                 sample.nobs = 1000)

# I: force convergence = TRUE, ML iterations = 150, identification = reference marker, bounds = "none"
mod1i.fit <- cfa(mod1, sample.cov = temp_mat, 
                 control = list(iter.max = 150),
                 optim.force.converged = T,
                 bounds = "none",
                 sample.nobs = 1000)

mod2i.fit <- cfa(mod2, sample.cov = temp_mat, 
                 control = list(iter.max = 150),
                 optim.force.converged = T,
                 bounds = "none",
                 sample.nobs = 1000)

mod3i.fit <- cfa(mod3, sample.cov = temp_mat, 
                 control = list(iter.max = 150),
                 optim.force.converged = T,
                 bounds = "none",
                 sample.nobs = 1000)

# J: force convergence = TRUE, ML iterations = 150, identification = reference marker, bounds = "pos.var"
mod1j.fit <- cfa(mod1, sample.cov = temp_mat, 
                 control = list(iter.max = 150),
                 optim.force.converged = T,
                 bounds = "pos.var",
                 sample.nobs = 1000)

mod2j.fit <- cfa(mod2, sample.cov = temp_mat, 
                 control = list(iter.max = 150),
                 optim.force.converged = T,
                 bounds = "pos.var",
                 sample.nobs = 1000)

mod3j.fit <- cfa(mod3, sample.cov = temp_mat, 
                 control = list(iter.max = 150),
                 optim.force.converged = T,
                 bounds = "pos.var",
                 sample.nobs = 1000)

# K: force convergence = TRUE, ML iterations = 150, identification = unit variance, bounds = "none"
mod1k.fit <- cfa(mod1, sample.cov = temp_mat, 
                 std.lv = T,
                 control = list(iter.max = 150),
                 optim.force.converged = T,
                 bounds = "none",
                 sample.nobs = 1000)

mod2k.fit <- cfa(mod2, sample.cov = temp_mat, 
                 std.lv = T,
                 control = list(iter.max = 150),
                 optim.force.converged = T,
                 bounds = "none",
                 sample.nobs = 1000)

mod3k.fit <- cfa(mod3, sample.cov = temp_mat, 
                 std.lv = T,
                 control = list(iter.max = 150),
                 optim.force.converged = T,
                 bounds = "none",
                 sample.nobs = 1000)

# L: force convergence = TRUE, ML iterations = 150, identification = unit variance, bounds = "pos.var"
mod1l.fit <- cfa(mod1, sample.cov = temp_mat, 
                 std.lv = T,
                 control = list(iter.max = 150),
                 optim.force.converged = T,
                 bounds = "pos.var",
                 sample.nobs = 1000)

mod2l.fit <- cfa(mod2, sample.cov = temp_mat, 
                 std.lv = T,
                 control = list(iter.max = 150),
                 optim.force.converged = T,
                 bounds = "pos.var",
                 sample.nobs = 1000)

mod3l.fit <- cfa(mod3, sample.cov = temp_mat,
                 std.lv = T,
                 control = list(iter.max = 150),
                 optim.force.converged = T,
                 bounds = "pos.var",
                 sample.nobs = 1000)

# M: force convergence = TRUE, ML iterations = 1500, identification = reference marker, bounds = "none"
mod1m.fit <- cfa(mod1, sample.cov = temp_mat, 
                 control = list(iter.max = 1500),
                 optim.force.converged = T,
                 bounds = "none",
                 sample.nobs = 1000)

mod2m.fit <- cfa(mod2, sample.cov = temp_mat, 
                 control = list(iter.max = 1500),
                 optim.force.converged = T,
                 bounds = "none",
                 sample.nobs = 1000)

mod3m.fit <- cfa(mod3, sample.cov = temp_mat, 
                 control = list(iter.max = 1500),
                 optim.force.converged = T,
                 bounds = "none",
                 sample.nobs = 1000)

# N: force convergence = TRUE, ML iterations = 1500, identification = reference marker, bounds = "pos.var"
mod1n.fit <- cfa(mod1, sample.cov = temp_mat, 
                 control = list(iter.max = 1500),
                 optim.force.converged = T,
                 bounds = "pos.var",
                 sample.nobs = 1000)

mod2n.fit <- cfa(mod2, sample.cov = temp_mat, 
                 control = list(iter.max = 1500),
                 optim.force.converged = T,
                 bounds = "pos.var",
                 sample.nobs = 1000)

mod3n.fit <- cfa(mod3, sample.cov = temp_mat, 
                 control = list(iter.max = 1500),
                 optim.force.converged = T,
                 bounds = "pos.var",
                 sample.nobs = 1000)

# O: force convergence = TRUE, ML iterations = 1500, identification = unit variance, bounds = "none"
mod1o.fit <- cfa(mod1, sample.cov = temp_mat, 
                 std.lv = T,
                 control = list(iter.max = 1500),
                 optim.force.converged = T,
                 bounds = "none",
                 sample.nobs = 1000)

mod2o.fit <- cfa(mod2, sample.cov = temp_mat, 
                 std.lv = T,
                 control = list(iter.max = 1500),
                 optim.force.converged = T,
                 bounds = "none",
                 sample.nobs = 1000)

mod3o.fit <- cfa(mod3, sample.cov = temp_mat, 
                 std.lv = T,
                 control = list(iter.max = 1500),
                 optim.force.converged = T,
                 bounds = "none",
                 sample.nobs = 1000)

# P: force convergence = TRUE, ML iterations = 1500, identification = unit variance, bounds = "pos.var"
mod1p.fit <- cfa(mod1, sample.cov = temp_mat, 
                 std.lv = T,
                 control = list(iter.max = 1500),
                 optim.force.converged = T,
                 bounds = "pos.var",
                 sample.nobs = 1000)

mod2p.fit <- cfa(mod2, sample.cov = temp_mat, 
                 std.lv = T,
                 control = list(iter.max = 1500),
                 optim.force.converged = T,
                 bounds = "pos.var",
                 sample.nobs = 1000)

mod3p.fit <- cfa(mod3, sample.cov = temp_mat,
                 std.lv = T,
                 control = list(iter.max = 1500),
                 optim.force.converged = T,
                 bounds = "pos.var",
                 sample.nobs = 1000)

dir.create("src")

models <- list(mod1a.fit, mod1b.fit,mod1c.fit,mod1d.fit,mod1e.fit,mod1f.fit,mod1g.fit, mod1h.fit,
               mod1i.fit,mod1j.fit,mod1k.fit,mod1l.fit,mod1m.fit,mod1n.fit,mod1o.fit,mod1p.fit,
               mod2a.fit, mod2b.fit,mod2c.fit,mod2d.fit,mod2e.fit,mod2f.fit,mod2g.fit, mod2h.fit,
               mod2i.fit,mod2j.fit,mod2k.fit,mod2l.fit,mod2m.fit,mod2n.fit,mod2o.fit,mod2p.fit,
               mod3a.fit, mod3b.fit,mod3c.fit,mod3d.fit,mod3e.fit,mod3f.fit,mod3g.fit, mod3h.fit,
               mod3i.fit,mod3j.fit,mod3k.fit,mod3l.fit,mod3m.fit,mod3n.fit,mod3o.fit,mod3p.fit)

save(mod1a.fit, mod1b.fit,mod1c.fit,mod1d.fit,mod1e.fit,mod1f.fit,mod1g.fit, mod1h.fit,
     mod1i.fit,mod1j.fit,mod1k.fit,mod1l.fit,mod1m.fit,mod1n.fit,mod1o.fit,mod1p.fit,
     mod2a.fit, mod2b.fit,mod2c.fit,mod2d.fit,mod2e.fit,mod2f.fit,mod2g.fit, mod2h.fit,
     mod2i.fit,mod2j.fit,mod2k.fit,mod2l.fit,mod2m.fit,mod2n.fit,mod2o.fit,mod2p.fit,
     mod3a.fit, mod3b.fit,mod3c.fit,mod3d.fit,mod3e.fit,mod3f.fit,mod3g.fit, mod3h.fit,
     mod3i.fit,mod3j.fit,mod3k.fit,mod3l.fit,mod3m.fit,mod3n.fit,mod3o.fit,mod3p.fit,
     file ="src/sim_study2_models.RData")

save(models, file = "src/sim_study2_models_list.RData")
