# OckhamSEM simulation study 3: All correlation matrices ----
# Some of the seeds x models in the original simulation caused issues, 
# trying different seeds for those final 50k-200k

# next 5 lines of code can be used locally to get the model numbers that have missing 25ks
# files <- list.files("results/all_correlations_res/output_all_3/")
# 
# tab <- table(str_split_i(files, "_", 3))
# as.numeric(names(tab[which(tab == 38)]))
# as.numeric(names(tab[which(tab == 32)]))

#miss2 <- c(2, 5, 8,11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47) #all are model 2
miss8 <- c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48) # all are model 3


# packages
pkgs <- c('ockhamSEM', 'parallel', 'doParallel', 'lavaan', 'rslurm', 'tidyverse')

# check if this packages have been installed
rlang::check_required(pkgs)

# loading them
lapply(pkgs, library, character.only = TRUE)

# Set up seed ----
pars <- readRDS("src/sim_study2_seeds.RDS")
pars <- pars[2305:2784,]

# this should be 160 rows
pars <- pars %>% filter((model %in% miss8))
# Set up correlation/covariance matrix matrix ----

# Set up temp matrix ----
p <- 9
temp_mat.a <- diag(p)
colnames(temp_mat.a) <- rownames(temp_mat.a) <- paste0("x", seq(1:p))

p <- 27
temp_mat.b <- diag(p)
colnames(temp_mat.b) <- rownames(temp_mat.b) <- paste0("x", seq(1:p))

p <- 45
temp_mat.c <- diag(p)
colnames(temp_mat.c) <- rownames(temp_mat.c) <- paste0("x", seq(1:p))


# Set up models ----
mod1 <- '
f1 =~ x1 + x2 + x3 
f2 =~ x4 + x5 + x6 
f3 =~ x7 + x8 + x9
'

mod2 <- '
f1 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
f2 =~ x10 + x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18
f3 =~ x19 + x20 + x21 + x22 + x23 + x24 + x25 + x26 + x27
'

mod3 <- '
f1 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13 + x14 + x15
f2 =~ x16 + x17 + x18 + x19 + x20 + x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + x30
f3 =~ x31 + x32 + x33 + x34 + x35 + x36 + x37 + x38 + x39 + x40 + x41 + x42 + x43 + x44 + x45
'

# Set up different conditions ----

load("src/sim_study3_models_list.RData")

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
      
      saveRDS(result, file = paste0("output/study3_allcor_", model, "_", seed, ".RDS"))
      print(paste0("Finished model ", model, ", seed ", seed))
      
      result
    })},
    error = function(e) {
      err <<- e
      err <- as.character(err)
      
      result <- list(model, err)
      
      saveRDS(result, file = paste0("output/study3_allcor_err_", model, "_", seed, ".RDS"))
      print(paste0("Did not finish model ", model, ", seed ", seed))
      
      result
    })
  
}

sopt <- list(
  time = '48:00:00', # 48hrs 0minutes
  partition='general',
  "mem-per-cpu" = 8192,   # other options #12288, #8192, #4096 #2048 #6144,
  "mail-type" = 'END,FAIL',
  "mail-user" = 'sdwfvv@umsystem.edu')

sjob <- slurm_apply(sim.fp, pars, jobname = 'sim3_all_oneseed_c',
                    nodes = 20, cpus_per_node = 8,
                    preschedule_cores = FALSE,
                    submit = FALSE,
                    slurm_options= sopt,
                    global_objects = c("mod1","mod2","mod3",
                                       "temp_mat.a", "temp_mat.b", "temp_mat.c",
                                       "models")
)


# creating a directory
dir.create("_rslurm_sim3_all_oneseed_c/output")

