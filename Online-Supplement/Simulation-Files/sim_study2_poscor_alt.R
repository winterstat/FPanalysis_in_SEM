#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
#-------------------------------------------------------------------------------

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

# packages
pkgs <- c('ockhamSEM', 'parallel', 'doParallel', 'lavaan', 'rslurm')

# check if this packages have been installed
rlang::check_required(pkgs)

# loading them
lapply(pkgs, library, character.only = TRUE)

# Set up seed ----
pars <- readRDS("src/sim_study2_seeds.RDS")
pars <- pars[1:24,]

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

load("src/sim_study2_models_list.RData")

# Set up fitting propensity function

sim.fp <- function(model, seed, R = 1000) {
  # time
  start <- Sys.time()
  
  # FP
  res <- run.fitprop(models[[model]],
                     fit.measure = c("logl", "cfi", "srmr", "rmsea"), 
                     rmethod = "onion", 
                     reps = R, 
                     onlypos = TRUE,
                     seed = seed)
  
  # computation time
  diff <- Sys.time() - start
  
  # results
  result <- list(model, diff, res$fit_list)
  
  saveRDS(result, file = paste0("output/study2_poscor_", model, "_", seed, ".RDS"))
  print(paste0("Finished model ", model, ", seed ", seed))
  
  return(result)
}

sopt <- list(
  partition='general',
  "mem-per-cpu" = 4096, #4096, #12288, #8192, #4096 #2048
  "mail-type" = 'END,FAIL',
  "mail-user" = 'btk2xm@umsystem.edu')

sjob <- slurm_apply(sim.fp, pars, jobname = 'sim2_pos0',
                    #ntasks = 48,
                    nodes = 4, cpus_per_node = 4,
                    preschedule_cores = FALSE,
                    submit = FALSE,
                    slurm_options= sopt,
                    global_objects = c("mod1","mod2","mod3",
                                       "temp_mat",
                                       "models")
)
