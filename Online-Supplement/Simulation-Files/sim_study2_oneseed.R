#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
Sys.Date()
#-------------------------------------------------------------------------------


# OckhamSEM simulation study 2: All correlation matrices ----
# One of the seeds in the original simulation caused issues, 
# trying a different seed for those final 25k


# packages
pkgs <- c('ockhamSEM', 'parallel', 'doParallel', 'lavaan', 'rslurm')

# check if this packages have been installed
rlang::check_required(pkgs)

# loading them
lapply(pkgs, library, character.only = TRUE)
# Set up seed ----
pars <- readRDS("src/sim_study2_seeds.RDS")
pars <- pars[1921:1968,]

# Set up correlation/covariance matrix matrix ----
p <- 24 # number of variables
temp_mat <- diag(p) # diagonal matrix
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

load("src/sim_study2_models_list.RData")

# Run fitting propensity analysis ----

# Set up fitprop function
sim.fp <- function(model, seed, R = 25000) {
  
  tryCatch({
    suppressWarnings({
      # time
      start <- Sys.time()
      
      # FP
      res <- run.fitprop(models[[model]],
                         fit.measure = c("logl", "cfi", "srmr", "rmsea"), 
                         rmethod = "onion", 
                         reps = R, 
                         onlypos = FALSE,
                         seed = seed)
      
      # computation time
      diff <- Sys.time() - start
      
      result <- list(model, diff, res$fit_list)
      
      saveRDS(result, file = paste0("output/study2_allcor_", model, "_", seed, ".RDS"))
      print(paste0("Finished model ", model, ", seed ", seed))
      
      result
    })},
    error = function(e) {
      err <<- e
      err <- as.character(err)
    
    result <- list(model, err)
    
    saveRDS(result, file = paste0("output/study2_allcor_err_", model, "_", seed, ".RDS"))
    print(paste0("Did not finish model ", model, ", seed ", seed))
    
    result
    })
  
}

sopt <- list(
  time = '48:00:00', # 48hrs 0minutes
  partition='general',
  "mem-per-cpu" = 4096,   # other options #12288, #8192, #4096 #2048
  "mail-type" = 'END,FAIL',
  "mail-user" = 'sdwfvv@umsystem.edu')

sjob <- slurm_apply(sim.fp, pars, jobname = 'sim2_all_oneseed',
                    nodes = 4, cpus_per_node = 12,
                    preschedule_cores = FALSE,
                    submit = FALSE,
                    slurm_options= sopt,
                    global_objects = c("mod1","mod2","mod3",
                                       "temp_mat",
                                       "models")
)


# creating a directory
dir.create("_rslurm_sim2_all_oneseeed/output")

