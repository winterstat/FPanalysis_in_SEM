#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
#-------------------------------------------------------------------------------

# OckhamSEM simulation study 1: All correlation matrices ----
# Exploring Introduction example preacher (3 variables)
# Total data space: 1 million correlation matrices
# Conditions included: 
#           - space size (post simulation)
#           - Models (2 models)
#           - All correlations vs only positive (two different data spaces)
# Measures of fit used to assess the fitting propensity
#           - LL: Log-Likelihood
#           - CFI: Comparative Fit Index
#           - SRMR : Standardized Root Mean Squared Residual
#           - RMSEA: Root Mean Square Error of Approximation 

# Load packages ----
# packages
pkgs <- c('ockhamSEM', 'parallel', 'doParallel', 'lavaan')

# check if this packages have been installed
rlang::check_required(pkgs)

# loading them
lapply(pkgs, library, character.only = TRUE)

# Set up correlation/covariance matrix matrix ----
p <- 3 # number of variables
temp_mat <- diag(p)
colnames(temp_mat) <- rownames(temp_mat) <- paste0("v", seq(1:p))

# Set up models ----
mod1 <- '
v3 ~ v1 + v2
v1 ~~ 0*v2  # correlation between v1 and v2 is fixed to zero.
'

mod2 <- '
v3 ~ v1
v2 ~ v3
'

# Set up different conditions ----

# A: force convergence = FALSE, ML iterations = default (150)
mod1.fit <- sem(mod1, sample.cov = temp_mat, sample.nobs = 1000)
mod2.fit <- sem(mod2, sample.cov = temp_mat, sample.nobs = 1000)


# Run fitting propensity analysis ----
# number of cores to use, detect Cores will give the total number of cores available
cores <- detectCores() - 2 
cl <- makeCluster(cores, type = 'PSOCK')
# registering the cluster
registerDoParallel(cl)

# number of correlation matrices
R <- 1e6 # 1 million 

# time
start <- Sys.time()


# Computing the fitting propensity 
res <- run.fitprop(mod1.fit, mod2.fit, 
                   fit.measure = c("logl", "cfi", "srmr", "rmsea"),
                   rmethod = "onion", reps = R,
                   onlypos = FALSE,# relaxing the positive only correlations
                   cluster = cl,
                   seed = 3242, saveR = T, saveModel = T)

# computation time
Sys.time() - start

# closing cluster
stopCluster(cl)

# save the results
saveRDS(res, file = 'results/study1_allcor.RDS')

