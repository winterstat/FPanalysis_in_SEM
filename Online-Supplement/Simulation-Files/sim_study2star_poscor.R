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
#           - variance bounds: "none" (default), or "pos.var"(positive variances)
#           - within positive correlation data space, we have: 3 (models) x 2 (MLit) x 2 (force conv) x 2 (ident) x 2 (bounds) = 48 conditions
# Measures of fit used to assess the fitting propensity
#           - LL: Log-Likelihood
#           - CFI: Comparative Fit Index
#           - SRMR : Standardized Root Mean Squared Residual
#           - RMSEA: Root Mean Square Error of Approximation 

# Load packages ----
library(ockhamSEM)
library(parallel)
library(doParallel)
library(lavaan) # though it is a dependencies of the ockhamSEM package (load it too)

# Set a seed
set.seed(20250107)

# Generate 1000 seeds
seeds <- round(runif(n = 1000, min = 1, max = 1e+6), 0)
models <- 1:48

pars <- expand.grid(models, seeds)
colnames(pars) <- c("model","seed")


# Set up seed ----
seeds <- pars$seed

i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(i)

seedi <- seeds[i]

# Set up temp matrix ----
p <- 24 # number of variables
temp_mat <- diag(p) # diagonal matrix
colnames(temp_mat) <- rownames(temp_mat) <- paste0("x", seq(1:p)) # labeling the rows and columns

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

load("src/sim_study2_models.RData")

# Run fitting propensity analysis ----
# number of cores to use, detect Cores will give the total number of cores available
cores <- detectCores() - 2 
cl <- makeCluster(cores, type = 'PSOCK')
# registering the cluster
registerDoParallel(cl)

# number of correlation matrices
R <- 25 # 1000


# time
start <- Sys.time()


# Computing the fitting propensity 
res <- run.fitprop(mod1a.fit, mod1b.fit, mod1c.fit, mod1d.fit, mod1e.fit, mod1f.fit, 
                   
                   mod1g.fit, mod1h.fit, mod1i.fit, mod1j.fit, mod1k.fit, mod1l.fit, 
                   
                   mod1m.fit, mod1n.fit, mod1o.fit, mod1p.fit, mod2a.fit, mod2b.fit,
                   
                   mod2c.fit, mod2d.fit, mod2e.fit, mod2f.fit, mod2g.fit, mod2h.fit,
                   
                   mod2i.fit, mod2j.fit, mod2k.fit, mod2l.fit, mod2m.fit, mod2n.fit,
                   
                   mod2o.fit, mod2p.fit, mod3a.fit, mod3b.fit, mod3c.fit, mod3d.fit,
                   
                   mod3e.fit, mod3f.fit, mod3g.fit, mod3h.fit, mod3i.fit, mod3j.fit,
                   
                   mod3k.fit, mod3l.fit, mod3m.fit, mod3n.fit, mod3o.fit, mod3p.fit,
                   
                   fit.measure = c("logl", "cfi", "srmr", "rmsea"), rmethod = "onion", 
                   
                   reps = R, onlypos = TRUE,
                   cluster = cl,
                   seed = seedi)

# computation time
diff <- Sys.time() - start

# closing cluster
stopCluster(cl)

# save the results
saveRDS(res, file = paste0('results/study2_poscor_',i,'.RDS'))

