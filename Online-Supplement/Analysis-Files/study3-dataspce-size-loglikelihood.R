#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 09/23/2025
#-------------------------------------------------------------------------------
# Relationship between the data space and the fitting propensity of the models
# in study 3
#-------------------------------------------------------------------------------

## Load packages ---------------------------------------------------------------
pkgs <- c('tidyverse', 'ggdist', 'ggthemes', 'effectsize', 'dplyr', 'ggplot2')

# check if packages are available
# rlang::check_installed(pkgs)

# loading packages
lapply(pkgs, library, character.only = TRUE)

# # create a local package library (run once)
# dir.create("C:/Rlibs", showWarnings = FALSE)
# 
# # put the new folder first in the library search path for this session
# .libPaths("C:/Rlibs")
# 
# # now install into the local folder
# install.packages("xfun", lib = .libPaths()[1], dependencies = TRUE)

## Load full data---------------------------------------------------------------
load(here::here("Online-Supplement","Data-Files", "output_poscor_3_all_conditions.RData"))

# positive correlations
study3_pos_fits <- output_conditions
load(here::here("Online-Supplement","Data-Files", "output_all_3_conditions.RData"))

study3_all_fits <- output_conditions

# removing the data
rm(output_conditions)
#-------------------------------------------------------------------------------

## Select range of data space sizes ----
sizes <- c(50, 100, 500, 1000, 5000, 10000, 50000, 100000, 150000, 175000, 200000, 500000, 750000)

# -Cutoffs----------------------------------------------------------------------
# Study 3 for all:
#   M1: -11460, M2: -36999; M3: -62493
# 
# Positive:
#   M1: -10,306, M2: -30,479; M3: -50,076

#-----------------------------------------------------
# Case I: Where we replace the non-convergence with some badness of fit value
#-------------------------------------------------------------------------------
# Positive correlations
#-------------------------------------------------------------------------------

# subsets
model_1 <- lapply(study3_pos_fits, function(x) x %>% filter(model == 1))

model_2 <- lapply(study3_pos_fits, function(x) x %>% filter(model == 2))

model_3 <- lapply(study3_pos_fits, function(x) x %>% filter(model == 3))

# removing the empty data frames
model_1 <- model_1[sapply(model_1,nrow) > 0]

model_2 <- model_2[sapply(model_2,nrow) > 0]

model_3 <- model_3[sapply(model_3,nrow) > 0]

  
# Model 1
# function to get prop CFI >= .51
cfi.fp_1 <- function(fits, cutoff = .51) {
  
  fits[is.na(fits)] <- 0
    
  return((sum(fits >= cutoff) / length(fits)))
}


# function to get prop RMSEA < .41
rmsea.fp_1 <- function(fits, cutoff = .41) {
  
  fits[is.na(fits)] <- 1
  
  return((sum(fits < cutoff) / length(fits)))
}


# function to get prop SRMR < .17
srmr.fp_1 <- function(fits, cutoff = .17) {
  
  fits[is.na(fits)] <- 1
  
  return((sum(fits < cutoff) / length(fits)))
}

# function to get the log likelihood with > -10,306---------------------------

logl.fp_1 <- function(fits, cutoff = -10306){
  
  fits[is.na(fits)] <- -30000
  
  return((sum(fits >= cutoff) / length(fits)))
  
}


# functions to sample X fit values and compute FP
samfits <- function(fits_all, size = 100) {
  
  rows <- sample(1:nrow(fits_all[[1]]), size, replace = FALSE)
  
  # container
  results <- matrix(NA, length(fits_all), 5)
  
  # renaming columns
  colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl')
  
  for(m in seq_along(fits_all)){
    
    # sub-setting  
    fits = fits_all[[m]][rows,]
    
    # CFI
    cfi<- cfi.fp_1(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp_1(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp_1(fits[,"srmr"])
    
    # LOGL
    logl <- logl.fp_1(fits[,"logl"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi, rmsea, srmr, logl)
    
    # saving results
    results[m,] <- res
    
  }
  
  # output
  return(results)
  
  
}
#-------------------------------------------------------------------------------
# Model 2
# function to get prop CFI >= .431
cfi.fp_2 <- function(fits, cutoff = .431) {
  
  fits[is.na(fits)] <- 0
  
  return((sum(fits >= cutoff) / length(fits)))
}


# function to get prop RMSEA < .249
rmsea.fp_2 <- function(fits, cutoff = .249) {
  
  fits[is.na(fits)] <- 1
  
  return((sum(fits < cutoff) / length(fits)))
}


# function to get prop SRMR < .084
srmr.fp_2 <- function(fits, cutoff = .084) {
  
  fits[is.na(fits)] <- 1
  
  return((sum(fits < cutoff) / length(fits)))
}

# function to get the log likelihood with > -30,479------------------------------

logl.fp_2 <- function(fits, cutoff = -30479){
  
  fits[is.na(fits)] <- -35000
  
  return((sum(fits >= cutoff) / length(fits)))
  
}


# functions to sample X fit values and compute FP
samfits_2 <- function(fits_all, size = 100) {
  
  rows <- sample(1:nrow(fits_all[[1]]), size, replace = FALSE)
  
  # container
  results <- matrix(NA, length(fits_all), 5)
  
  # renaming columns
  colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl')
  
  for(m in seq_along(fits_all)){
    # sub-setting  
    fits = fits_all[[m]][rows,]
    
    # CFI
    cfi<- cfi.fp_2(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp_2(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp_2(fits[,"srmr"])
    
    # LOGL
    logl <- logl.fp_2(fits[,"logl"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi, rmsea, srmr, logl)
    
    # saving results
    results[m,] <- res
    
  }
  
  # output
  return(results)
}  
#-------------------------------------------------------------------------------

# Model 3
# function to get prop CFI >= .35
cfi.fp_3 <- function(fits, cutoff = .35) {
  
  fits[is.na(fits)] <- 0
  
  return((sum(fits >= cutoff) / length(fits)))
}


# function to get prop RMSEA < .17
rmsea.fp_3 <- function(fits, cutoff = .17) {
  
  fits[is.na(fits)] <- 1
  
  return((sum(fits < cutoff) / length(fits)))
}


# function to get prop SRMR < .065
srmr.fp_3 <- function(fits, cutoff = .065) {
  
  fits[is.na(fits)] <- 1
  
  return((sum(fits < cutoff) / length(fits)))
}

# function to get the log likelihood with > -50,076------------------------------

logl.fp_3 <- function(fits, cutoff = -50076){
  
  fits[is.na(fits)] <- -55000
  
  return((sum(fits >= cutoff) / length(fits)))
  
}


# functions to sample X fit values and compute FP
samfits_3 <- function(fits_all, size = 100) {
  
  rows <- sample(1:nrow(fits_all[[1]]), size, replace = FALSE)
  
  # container
  results <- matrix(NA, length(fits_all), 5)
  
  # renaming columns
  colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl')
  
  for(m in seq_along(fits_all)){
    # sub-setting  
    fits = fits_all[[m]][rows,]
    
    # CFI
    cfi<- cfi.fp_3(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp_3(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp_3(fits[,"srmr"])
    
    # LOGL
    logl <- logl.fp_3(fits[,"logl"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi, rmsea, srmr, logl)
    
    # saving results
    results[m,] <- res
    
  }
  
  # output
  return(results)
}  



#-------------------------------------------------------------------------------
## Get FP for different space sizes ----

# Model 1
# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits(model_1, size = size)))
  tmp <- (do.call(rbind, tmp))
  tmp <- cbind(tmp, rep(sizes, each = 16))
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl', 'size')
  
  if(i == 1) {
    fp.size.results_1 <- tmp
  } else {
    fp.size.results_1 <- rbind(fp.size.results_1, tmp)
  }
}


# Grid of conditions that we can use for post-processing ----
study3_conditions <- expand.grid(model = 1:3, 
                                 bounds = c("none", "pos.var"),
                                 ident = c("refmarker", "unitvar"),
                                 mlit = c(150, 1500),
                                 conv = c("F", "T"))

# Individual study conditions
study_con_model_1 <- study3_conditions %>%
  filter(model == 1)

# merging the data frame with the conditions
study_con_model_1$condition <- 1:16

res_pos_model_1 <- merge(fp.size.results_1, study_con_model_1)

res_pos_model_1$model <- factor(res_pos_model_1$model)

# saving the results
saveRDS(fp.size.results_1, file = here::here("Online-Supplement","Data-Files", "study3_poscor_model_1_dataspacesize.RDS"))

saveRDS(res_pos_model_1, file = here::here("Online-Supplement", "Data-Files", "results_poscor_study_3_model_1_dataspacesize.RDS"))

# Model 2

# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits_2(model_2, size = size)))
  tmp <- (do.call(rbind, tmp))
  tmp <- cbind(tmp, rep(sizes, each = 16))
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl', 'size')
  
  if(i == 1) {
    fp.size.results_2 <- tmp
  } else {
    fp.size.results_2 <- rbind(fp.size.results_2, tmp)
  }
}


# Individual study conditions
study_con_model_2 <- study3_conditions %>%
  filter(model == 2)

# merging the data frame with the conditions
study_con_model_2$condition <- 1:16

res_pos_model_2 <- merge(fp.size.results_2, study_con_model_2)

res_pos_model_2$model <- factor(res_pos_model_2$model)

# saving the results
saveRDS(fp.size.results_2, file = here::here("Online-Supplement", "Data-Files", "study3_poscor_model_2_dataspacesize.RDS"))

saveRDS(res_pos_model_2, file = here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_2_dataspacesize.RDS"))

# Model 3

# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits_3(model_3, size = size)))
  tmp <- (do.call(rbind, tmp))
  tmp <- cbind(tmp, rep(sizes, each = 16))
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr','logl', 'size')
  
  if(i == 1) {
    fp.size.results_3 <- tmp
  } else {
    fp.size.results_3 <- rbind(fp.size.results_3, tmp)
  }
}

# Individual study conditions
study_con_model_3 <- study3_conditions %>%
  filter(model == 3)

# merging the data frame with the conditions
study_con_model_3$condition <- 1:16

res_pos_model_3 <- merge(fp.size.results_3, study_con_model_3)

res_pos_model_3$model <- factor(res_pos_model_3$model)

# saving the results
saveRDS(fp.size.results_3, file = here::here("Online-Supplement", "Data-Files", "study3_poscor_model_3_dataspacesize.RDS"))
saveRDS(res_pos_model_3, file = here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_3_dataspacesize.RDS"))

#-------------------------------------------------------------------------------
# All Correlations
#-------------------------------------------------------------------------------


# subsets
model_1 <- lapply(study3_all_fits, function(x) x %>% filter(model == 1))

model_2 <- lapply(study3_all_fits, function(x) x %>% filter(model == 2))

model_3 <- lapply(study3_all_fits, function(x) x %>% filter(model == 3))

# removing the empty data frames
model_1 <- model_1[sapply(model_1,nrow) > 0]

model_2 <- model_2[sapply(model_2,nrow) > 0]

model_3 <- model_3[sapply(model_3,nrow) > 0]


# Model 1
# function to get prop CFI >= .51
cfi.fp_1 <- function(fits, cutoff = .51) {
  
  fits[is.na(fits)] <- 0
    
    return((sum(fits >= cutoff) / length(fits)))
}


# function to get prop RMSEA < .41
rmsea.fp_1 <- function(fits, cutoff = .41) {
  
  fits[is.na(fits)] <- 1
  
  return((sum(fits < cutoff) / length(fits)))
}


# function to get prop SRMR < .17
srmr.fp_1 <- function(fits, cutoff = .17) {
  
  fits[is.na(fits)] <- 1
  
  return((sum(fits < cutoff) / length(fits)))
}

# function to get the log likelihood with > -11,460---------------------------

logl.fp_1 <- function(fits, cutoff = -11460){
  
  fits[is.na(fits)] <- -30000
  
  return((sum(fits >= cutoff) / length(fits)))
  
}


# functions to sample X fit values and compute FP
samfits <- function(fits_all, size = 100) {
  
  rows <- sample(1:nrow(fits_all[[1]]), size, replace = FALSE)
  
  # container
  results <- matrix(NA, length(fits_all), 5)
  
  # renaming columns
  colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl')
  
  for(m in seq_along(fits_all)){
    
    # sub-setting  
    fits = fits_all[[m]][rows,]
    
    # CFI
    cfi<- cfi.fp_1(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp_1(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp_1(fits[,"srmr"])
    
    # LOGL
    logl <- logl.fp_1(fits[,"logl"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi, rmsea, srmr, logl)
    
    # saving results
    results[m,] <- res
    
  }
  
  # output
  return(results)
  
  
}
#-------------------------------------------------------------------------------
# Model 2
# function to get prop CFI >= .431
cfi.fp_2 <- function(fits, cutoff = .431) {
  
  fits[is.na(fits)] <- 0
  
  return((sum(fits >= cutoff) / length(fits)))
}


# function to get prop RMSEA < .249
rmsea.fp_2 <- function(fits, cutoff = .249) {
  
  fits[is.na(fits)] <- 1
  
  return((sum(fits < cutoff) / length(fits)))
}


# function to get prop SRMR < .09
srmr.fp_2 <- function(fits, cutoff = .084) {
  
  fits[is.na(fits)] <- 1
  
  return((sum(fits < cutoff) / length(fits)))
}

# function to get the log likelihood with > -36,999------------------------------

logl.fp_2 <- function(fits, cutoff = -36999){
  
  fits[is.na(fits)] <- -40000
  
  return((sum(fits >= cutoff) / length(fits)))
  
}


# functions to sample X fit values and compute FP
samfits_2 <- function(fits_all, size = 100) {
  
  rows <- sample(1:nrow(fits_all[[1]]), size, replace = FALSE)
  
  # container
  results <- matrix(NA, length(fits_all), 5)
  
  # renaming columns
  colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl')
  
  for(m in seq_along(fits_all)){
    # sub-setting  
    fits = fits_all[[m]][rows,]
    
    # CFI
    cfi<- cfi.fp_2(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp_2(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp_2(fits[,"srmr"])
    
    # LOGL
    logl <- logl.fp_2(fits[,"logl"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi, rmsea, srmr, logl)
    
    # saving results
    results[m,] <- res
    
  }
  
  # output
  return(results)
}  
#-------------------------------------------------------------------------------

# Model 3
# function to get prop CFI >= .35
cfi.fp_3 <- function(fits, cutoff = .35) {
  
  fits[is.na(fits)] <- 0
  
  return((sum(fits >= cutoff) / length(fits)))
}


# function to get prop RMSEA < .17
rmsea.fp_3 <- function(fits, cutoff = .17) {
  
  fits[is.na(fits)] <- 1
  
  return((sum(fits < cutoff) / length(fits)))
}


# function to get prop SRMR < .065
srmr.fp_3 <- function(fits, cutoff = .065) {
  
  fits[is.na(fits)] <- 1
  
  return((sum(fits < cutoff) / length(fits)))
}

# function to get the log likelihood with > -62,493------------------------------

logl.fp_3 <- function(fits, cutoff = -62493){
  
  fits[is.na(fits)] <- -70000
  
  return((sum(fits >= cutoff) / length(fits)))
  
}


# functions to sample X fit values and compute FP
samfits_3 <- function(fits_all, size = 100) {
  
  rows <- sample(1:nrow(fits_all[[1]]), size, replace = FALSE)
  
  # container
  results <- matrix(NA, length(fits_all), 5)
  
  # renaming columns
  colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl')
  
  for(m in seq_along(fits_all)){
    # sub-setting  
    fits = fits_all[[m]][rows,]
    
    # CFI
    cfi<- cfi.fp_3(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp_3(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp_3(fits[,"srmr"])
    
    # LOGL
    logl <- logl.fp_3(fits[,"logl"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi, rmsea, srmr, logl)
    
    # saving results
    results[m,] <- res
    
  }
  
  # output
  return(results)
}  



#-------------------------------------------------------------------------------
## Get FP for different space sizes ----

# Model 1
# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits(model_1, size = size)))
  tmp <- (do.call(rbind, tmp))
  tmp <- cbind(tmp, rep(sizes, each = 16))
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl', 'size')
  
  if(i == 1) {
    fp.size.results_1 <- tmp
  } else {
    fp.size.results_1 <- rbind(fp.size.results_1, tmp)
  }
}


# Grid of conditions that we can use for post-processing ----
study3_conditions <- expand.grid(model = 1:3, 
                                 bounds = c("none", "pos.var"),
                                 ident = c("refmarker", "unitvar"),
                                 mlit = c(150, 1500),
                                 conv = c("F", "T"))

# Individual study conditions
study_con_model_1 <- study3_conditions %>%
  filter(model == 1)

# merging the data frame with the conditions
study_con_model_1$condition <- 1:16

res_all_model_1 <- merge(fp.size.results_1, study_con_model_1)

res_all_model_1$model <- factor(res_all_model_1$model)

# saving the results
saveRDS(fp.size.results_1, file = here::here("Online-Supplement","Data-Files", "study3_allcor_model_1_dataspacesize.RDS"))
saveRDS(res_all_model_1, file = here::here("Online-Supplement","Data-Files", "results_allcor_study_3_model_1_dataspacesize.RDS"))

# Model 2

# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits_2(model_2, size = size)))
  tmp <- (do.call(rbind, tmp))
  tmp <- cbind(tmp, rep(sizes, each = 16))
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl', 'size')
  
  if(i == 1) {
    fp.size.results_2 <- tmp
  } else {
    fp.size.results_2 <- rbind(fp.size.results_2, tmp)
  }
}


# Individual study conditions
study_con_model_2 <- study3_conditions %>%
  filter(model == 2)

# merging the data frame with the conditions
study_con_model_2$condition <- 1:16

res_all_model_2 <- merge(fp.size.results_2, study_con_model_2)

res_all_model_2$model <- factor(res_all_model_2$model)

# saving the results
saveRDS(fp.size.results_2, file = here::here("Online-Supplement","Data-Files", "study3_allcor_model_2_dataspacesize.RDS"))
saveRDS(res_all_model_2, file = here::here("Online-Supplement","Data-Files", "results_allcor_study_3_model_2_dataspacesize.RDS"))

# Model 3

# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits_3(model_3, size = size)))
  tmp <- (do.call(rbind, tmp))
  tmp <- cbind(tmp, rep(sizes, each = 16))
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr','logl', 'size')
  
  if(i == 1) {
    fp.size.results_3 <- tmp
  } else {
    fp.size.results_3 <- rbind(fp.size.results_3, tmp)
  }
}

# Individual study conditions
study_con_model_3 <- study3_conditions %>%
  filter(model == 3)

# merging the data frame with the conditions
study_con_model_3$condition <- 1:16

res_all_model_3 <- merge(fp.size.results_3, study_con_model_3)

res_all_model_3$model <- factor(res_all_model_3$model)

# saving the results
saveRDS(fp.size.results_3, file = here::here("Online-Supplement","Data-Files", "study3_allcor_model_3_dataspacesize.RDS"))
saveRDS(res_all_model_3, file = here::here("Online-Supplement","Data-Files", "results_allcor_study_3_model_3_dataspacesize.RDS"))


# clean work space
rm(list = ls())

#-------------------------------------------------------------------------------
# Analysis:
# loading the datasets
res_pos_model_1 <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_1_dataspacesize.RDS"))
res_all_model_1 <- readRDS(here::here("Online-Supplement","Data-Files", "results_allcor_study_3_model_1_dataspacesize.RDS"))

res_pos_model_2<- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_2_dataspacesize.RDS"))
res_all_model_2 <- readRDS(here::here("Online-Supplement","Data-Files", "results_allcor_study_3_model_2_dataspacesize.RDS"))

res_pos_model_3 <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_3_dataspacesize.RDS"))
res_all_model_3 <- readRDS(here::here("Online-Supplement","Data-Files", "results_allcor_study_3_model_3_dataspacesize.RDS"))
#-------------------------------------------------------------------------------


# new variable in each data frame-----------------------------------------------
# Model 1
results_full_1 <- bind_rows(
  
  res_pos_model_1 %>% mutate(cor = 'pos') %>% mutate(cor = factor(cor)),
  
  res_all_model_1 %>% mutate(cor = 'all') %>% mutate(cor = factor(cor))
)

# Model 2
results_full_2 <- bind_rows(
  
  res_pos_model_2 %>% mutate(cor = 'pos') %>% mutate(cor = factor(cor)),
  
  res_all_model_2 %>% mutate(cor = 'all') %>% mutate(cor = factor(cor))
)

# Model 3
results_full_3 <- bind_rows(
  
  res_pos_model_3 %>% mutate(cor = 'pos') %>% mutate(cor = factor(cor)),
  
  res_all_model_3 %>% mutate(cor = 'all') %>% mutate(cor = factor(cor))
)

# checking for any missing data so we can confirm that our merging strategy worked
apply(results_full_1, 2, function(x) 100* mean(is.na(x)))

apply(results_full_2, 2, function(x) 100* mean(is.na(x)))

apply(results_full_3, 2, function(x) 100* mean(is.na(x)))

# Fitting propensity------------------------------------------------------------
results_full_1 %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(size == 1e4 & cor == "all") %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size, index, bounds,ident, mlit, conv) %>%
  summarize(m = mean(value),
            se = sd(value), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>%
  pivot_wider(names_from = c("cor", "model","conv", "mlit"), values_from = result) -> table_11

results_full_2 %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(size == 1e4 & cor == "all") %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size, index, bounds,ident, mlit, conv) %>%
  summarize(m = mean(value),
            se = sd(value), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>%
  pivot_wider(names_from = c("cor", "model","conv", "mlit"), values_from = result) -> table_22

results_full_3 %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(size == 1e4 & cor == "all") %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size, index, bounds,ident, mlit, conv) %>%
  summarize(m = mean(value),
            se = sd(value), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>%
  pivot_wider(names_from = c("cor", "model","conv", "mlit"), values_from = result) -> table_33
#-------------------------------------------------------------------------------
## Total fitting propensity for each fit measure
results_full_1 %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(size == 1e4 & cor == "all") %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size, index, bounds,ident, mlit, conv) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop")%>% 
  pivot_wider(names_from = c("cor", "model","conv", "mlit"), values_from = total) -> total_11

results_full_2 %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(size == 1e4 & cor == "all") %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size, index, bounds,ident, mlit, conv) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop")%>% 
  pivot_wider(names_from = c("cor", "model","conv", "mlit"), values_from = total) -> total_22

results_full_3 %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(size == 1e4 & cor == "all") %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size, index, bounds,ident, mlit, conv) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop")%>% 
  pivot_wider(names_from = c("cor", "model","conv", "mlit"), values_from = total) -> total_33
#-------------------------------------------------------------------------------
# Saving the results above
saveRDS(table_11, file = here::here("Online-Supplement","Data-Files", "table_11.RDS"))
saveRDS(table_22, file = here::here("Online-Supplement","Data-Files", "table_22.RDS"))
saveRDS(table_33, file = here::here("Online-Supplement","Data-Files", "table_33.RDS"))

saveRDS(total_11, file = here::here("Online-Supplement","Data-Files", "total_11.RDS"))
saveRDS(total_22, file = here::here("Online-Supplement","Data-Files", "total_22.RDS"))
saveRDS(total_33, file = here::here("Online-Supplement","Data-Files", "total_33.RDS"))

#-------------------------------------------------------------------------------
# Plots
#-------------------------------------------------------------------------------

# setting the theme
my_theme <- theme_grey() + theme(panel.background = element_blank(),
            legend.position = 'bottom',
            axis.title = element_text(face = "bold"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), strip.text = element_text(face = 'bold', size = 6.0),
            strip.background = element_rect(colour = 'grey79'), 
            axis.line = element_line(colour = 'grey79', linewidth = .35),
            axis.text = element_text(size = 6.0)) 

theme_set(my_theme)

#--------------------------------------------------------------------------------
results_full_1 %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  filter(size <= 10000) %>% 
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c('Model 1C', 'Model 2C', 'Model 3C')),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr", "logl"), labels = c("CFI", "RMSEA", "SRMR", "LL")), 
         
         bounds = factor(bounds, levels = c('none', 'pos.var'), labels = c('Unrestricted', 'Positive')),
         
         ident = factor(ident, levels = c('refmarker', 'unitvar'), labels = c('Reference', 'Unit Variance')),
         
         mlit = factor(mlit, levels = c('150', '1500')),
         
         conv = factor(conv, levels = c('F', 'T'), labels = c('Default', 'Forced')),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == 'LL' & cor == 'Positive Subspace') %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  linewidth = .1, alpha = 1) +
  facet_grid(rows = vars(mlit),
             cols =vars(ident, conv),
             scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (proportion  > cutoff)",
       fill = 'Model', fill_ramp = 'Confidence Level', linetype = 'Model') +
  scale_fill_manual(name = 'Model',
                   values = c('Model 1C' = "#276419",
                              'Model 2C' = "#8e0152", 
                              'Model 3C' = "#b2182b")) +
  scale_linetype_manual(
    name = "Model",
    values = c("Model 1C" = "solid", "Model 2C" = "dashed", "Model 3C" = "dotted")) +
  scale_x_continuous(sec.axis = sec_axis(~ . , 
                                         name = "Positive Subspace", breaks = NULL, labels = NULL)) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Maximum Iterations", 
                                         breaks = NULL, labels = NULL)) -> plot_1

# saving the graph
ggsave(plot = plot_1, path = here::here("Online-Supplement", "Plots", "Study Three"), 
filename = 'Fit-positive-LL-dataspace-study-3-model-1.png', width = 10, height = 8,
       units = 'in')


#---
results_full_1 %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  filter(size <= 10000) %>% 
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c('Model 1C', 'Model 2C', 'Model 3C')),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr", "logl"), labels = c("CFI", "RMSEA", "SRMR", "LL")), 
         
         bounds = factor(bounds, levels = c('none', 'pos.var'), labels = c('Unrestricted', 'Positive')),
         
         ident = factor(ident, levels = c('refmarker', 'unitvar'), labels = c('Reference', 'Unit Variance')),
         
         mlit = factor(mlit, levels = c('150', '1500')),
         
         conv = factor(conv, levels = c('F', 'T'), labels = c('Default', 'Forced')),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == 'LL' & cor == 'Complete Space') %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  linewidth = .1, alpha = 1) +
  facet_grid(rows = vars(mlit),
             cols =vars(ident, conv),
             scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (proportion  > cutoff)",
       fill = 'Model', fill_ramp = 'Confidence Level', linetype = 'Model') +
  scale_fill_manual(name = 'Model',
                   values = c('Model 1C' = "#276419",
                              'Model 2C' = "#8e0152", 
                              'Model 3C' = "#b2182b")) +
  scale_linetype_manual(
    name = "Model",
    values = c("Model 1C" = "solid", "Model 2C" = "dashed", "Model 3C" = "dotted")) +
  scale_x_continuous(sec.axis = sec_axis(~ . , 
                                         name = "Complete Space", breaks = NULL, labels = NULL)) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Maximum Iterations", 
                                         breaks = NULL, labels = NULL)) -> plot_11

# saving the graph
ggsave(plot = plot_11, path = here::here("Online-Supplement","Plots", "Study Three"),
filename = 'Fit-all-LL-dataspace-study-3-model-1.png', width = 10, height = 8,
       units = 'in')

#--------------------------------------------------------------------------------
results_full_2 %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  filter(size <= 10000) %>% # selecting fit values less that the value specified
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c('Model 1C', 'Model 2C', 'Model 3C')),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr", "logl"), labels = c("CFI", "RMSEA", "SRMR", "LL")), 
         
         bounds = factor(bounds, levels = c('none', 'pos.var'), labels = c('Unrestricted', 'Positive')),
         
         ident = factor(ident, levels = c('refmarker', 'unitvar'), labels = c('Reference', 'Unit Variance')),
         
         mlit = factor(mlit, levels = c('150', '1500')),
         
         conv = factor(conv, levels = c('F', 'T'), labels = c('Default', 'Forced')),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == 'LL' & cor == 'Positive Subspace') %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  linewidth = .1, alpha = 1) +
  facet_grid(rows = vars(mlit),
             cols =vars(ident, conv),
             scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (proportion  > cutoff)",
       fill = 'Model', fill_ramp = 'Confidence Level', linetype = 'Model') +
  scale_fill_manual(name = 'Model',
                   values = c('Model 1C' = "#276419",
                              'Model 2C' = "#8e0152", 
                              'Model 3C' = "#b2182b")) +
  scale_linetype_manual(
    name = "Model",
    values = c("Model 1C" = "solid", "Model 2C" = "dashed", "Model 3C" = "dotted")) +
  scale_x_continuous(sec.axis = sec_axis(~ . , 
                                         name = "Positive Subspace", breaks = NULL, labels = NULL)) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Maximum Iterations", 
                                         breaks = NULL, labels = NULL)) -> plot_2

# saving the graph
ggsave(plot = plot_2, path = here::here("Online-Supplement","Plots", "Study Three"), 
filename = 'Fit-positive-LL-dataspace-study-3-model-2.png', width = 10, height = 8,
units = 'in')


#---
results_full_2 %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  filter(size <= 10000) %>% # selecting fit values less that the value specified
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c('Model 1C', 'Model 2C', 'Model 3C')),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr", "logl"), labels = c("CFI", "RMSEA", "SRMR", "LL")), 
         
         bounds = factor(bounds, levels = c('none', 'pos.var'), labels = c('Unrestricted', 'Positive')),
         
         ident = factor(ident, levels = c('refmarker', 'unitvar'), labels = c('Reference', 'Unit Variance')),
         
         mlit = factor(mlit, levels = c('150', '1500')),
         
         conv = factor(conv, levels = c('F', 'T'), labels = c('Default', 'Forced')),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Subspace"))) %>%
  filter(index == 'LL' & cor == 'Complete Subspace') %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  linewidth = .1, alpha = 1) +
  facet_grid(rows = vars(mlit),
             cols =vars(ident, conv),
             scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (proportion  > cutoff)",
       fill = 'Model', fill_ramp = 'Confidence Level', linetype = 'Model') +
  scale_fill_manual(name = 'Model',
                   values = c('Model 1C' = "#276419",
                              'Model 2C' = "#8e0152", 
                              'Model 3C' = "#b2182b")) +
  scale_linetype_manual(
    name = "Model",
    values = c("Model 1C" = "solid", "Model 2C" = "dashed", "Model 3C" = "dotted")) +
  scale_x_continuous(sec.axis = sec_axis(~ . , 
                                         name = "Complete Space", breaks = NULL, labels = NULL)) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Maximum Iterations", 
                                         breaks = NULL, labels = NULL)) -> plot_22

# saving the graph
ggsave(plot = plot_22, path = here::here("Online-Supplement","Plots", "Study Three"), 
filename = 'Fit-all-LL-dataspace-study-3-model-2.png', width = 10, height = 8,
units = 'in')

#--------------------------------------------------------------------------------
results_full_3 %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  filter(size <= 10000) %>% # selecting fit values less that the value specified
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c('Model 1C', 'Model 2C', 'Model 3C')),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr", "logl"), labels = c("CFI", "RMSEA", "SRMR", "LL")), 
         
         bounds = factor(bounds, levels = c('none', 'pos.var'), labels = c('Unrestricted', 'Positive')),
         
         ident = factor(ident, levels = c('refmarker', 'unitvar'), labels = c('Reference', 'Unit Variance')),
         
         mlit = factor(mlit, levels = c('150', '1500')),
         
         conv = factor(conv, levels = c('F', 'T'), labels = c('Default', 'Forced')),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == 'LL' & cor == 'Positive Subspace') %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  linewidth = .1, alpha = 1) +
  facet_grid(rows = vars(mlit),
             cols =vars(ident, conv),
             scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (proportion  > cutoff)",
       fill = 'Model', fill_ramp = 'Confidence Level', linetype = 'Model') +
  scale_fill_manual(name = 'Model',
                   values = c('Model 1C' = "#276419",
                              'Model 2C' = "#8e0152", 
                              'Model 3C' = "#b2182b")) +
  scale_linetype_manual(
    name = "Model",
    values = c("Model 1C" = "solid", "Model 2C" = "dashed", "Model 3C" = "dotted")) +
  scale_x_continuous(sec.axis = sec_axis(~ . , 
                                         name = "Positive Subspace", breaks = NULL, labels = NULL)) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Maximum Iterations", 
                                         breaks = NULL, labels = NULL)) -> plot_3

# saving the graph
ggsave(plot = plot_3, path = here::here("Online-Supplement","Plots", "Study Three"), 
filename = 'Fit-positive-LL-dataspace-study-3-model-3.png', width = 10, height = 8,
units = 'in')


#---
results_full_3 %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  filter(size <= 10000) %>% # selecting fit values less that the value specified
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c('Model 1C', 'Model 2C', 'Model 3C')),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr", "logl"), labels = c("CFI", "RMSEA", "SRMR", "LL")), 
         
         bounds = factor(bounds, levels = c('none', 'pos.var'), labels = c('Unrestricted', 'Positive')),
         
         ident = factor(ident, levels = c('refmarker', 'unitvar'), labels = c('Reference', 'Unit Variance')),
         
         mlit = factor(mlit, levels = c('150', '1500')),
         
         conv = factor(conv, levels = c('F', 'T'), labels = c('Default', 'Forced')),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == 'LL' & cor == 'Complete Space') %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  linewidth = .1, alpha = 1) +
  facet_grid(rows = vars(mlit),
             cols =vars(ident, conv),
             scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (proportion  > cutoff)",
       fill = 'Model', fill_ramp = 'Confidence Level', linetype = 'Model') +
  scale_fill_manual(name = 'Model',
                   values = c('Model 1C' = "#276419",
                              'Model 2C' = "#8e0152", 
                              'Model 3C' = "#b2182b")) +
  scale_linetype_manual(
    name = "Model",
    values = c("Model 1C" = "solid", "Model 2C" = "dashed", "Model 3C" = "dotted")) +
  scale_x_continuous(sec.axis = sec_axis(~ . , 
                                         name = " Complete Space", breaks = NULL, labels = NULL)) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Maximum Iterations",
                                         breaks = NULL, labels = NULL)) -> plot_33

# saving the graph
ggsave(plot = plot_33, path = here::here("Online-Supplement","Plots", "Study Three"),
filename = 'Fit-all-LL-dataspace-study-3-model-3.png', width = 10, height = 8,
units = 'in')

#-------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
# cleaning the environment
rm(list = ls())
#-------------------------------------------------------------------------------

# Case II: Using list-wise deletion---------------------------------------------

## Load full data---------------------------------------------------------------
load("Online-Supplement","Data-Files", "output_poscor_3_all_conditions.RData")

# positive correlations
study3_pos_fits <- output_conditions

load("Online-Supplement","Data-Files", "output_all_3_conditions.RData")

study3_all_fits <- output_conditions

# removing the data
rm(output_conditions)

# Select range of data space sizes ----
sizes <- c(50, 100, 500, 1000, 5000, 10000, 50000, 100000, 150000, 175000, 200000, 500000, 750000)
#-------------------------------------------------------------------------------
# Positive correlations
#-------------------------------------------------------------------------------

# subsets
model_1 <- lapply(study3_pos_fits, function(x) x %>% filter(model == 1))

model_2 <- lapply(study3_pos_fits, function(x) x %>% filter(model == 2))

model_3 <- lapply(study3_pos_fits, function(x) x %>% filter(model == 3))

# removing the empty data frames
model_1 <- model_1[sapply(model_1,nrow) > 0]

model_2 <- model_2[sapply(model_2,nrow) > 0]

model_3 <- model_3[sapply(model_3,nrow) > 0]

# Model 1
# function to get prop CFI >= .51
cfi.fp_1 <- function(fits, cutoff = .51) {
    
    return((sum(fits >= cutoff,  na.rm = TRUE) / length(fits)))
}


# function to get prop RMSEA < .41
rmsea.fp_1 <- function(fits, cutoff = .41) {

  return((sum(fits < cutoff,na.rm = TRUE) / length(fits)))
}


# function to get prop SRMR < .17
srmr.fp_1 <- function(fits, cutoff = .17) {
  
  return((sum(fits < cutoff, na.rm = TRUE) / length(fits)))
}

# function to get the log likelihood with > -10,306---------------------------

logl.fp_1 <- function(fits, cutoff = -10306){
  
  return((sum(fits >= cutoff, na.rm = TRUE) / length(fits)))
  
}


# functions to sample X fit values and compute FP
samfits <- function(fits_all, size = 100) {
  
  rows <- sample(1:nrow(fits_all[[1]]), size, replace = FALSE)
  
  # container
  results <- matrix(NA, length(fits_all), 5)
  
  # renaming columns
  colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl')
  
  for(m in seq_along(fits_all)){
    
    # sub-setting  
    fits = fits_all[[m]][rows,]
    
    # CFI
    cfi<- cfi.fp_1(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp_1(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp_1(fits[,"srmr"])
    
    # LOGL
    logl <- logl.fp_1(fits[,"logl"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi, rmsea, srmr, logl)
    
    # saving results
    results[m,] <- res
    
  }
  
  # output
  return(results)
  
  
}
#-------------------------------------------------------------------------------
# Model 2
# function to get prop CFI >= .431
cfi.fp_2 <- function(fits, cutoff = .431) {
  
  return((sum(fits >= cutoff,  na.rm = TRUE) / length(fits)))
}


# function to get prop RMSEA < .249
rmsea.fp_2 <- function(fits, cutoff = .249) {
  
  return((sum(fits < cutoff,  na.rm = TRUE) / length(fits)))
}


# function to get prop SRMR < .09
srmr.fp_2 <- function(fits, cutoff = .084) {
  
  return((sum(fits < cutoff,  na.rm = TRUE) / length(fits)))
}

# function to get the log likelihood with > -30,479------------------------------

logl.fp_2 <- function(fits, cutoff = -30479){
  
  return((sum(fits >= cutoff, na.rm = TRUE) / length(fits)))
  
}


# functions to sample X fit values and compute FP
samfits_2 <- function(fits_all, size = 100) {
  
  rows <- sample(1:nrow(fits_all[[1]]), size, replace = FALSE)
  
  # container
  results <- matrix(NA, length(fits_all), 5)
  
  # renaming columns
  colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl')
  
  for(m in seq_along(fits_all)){
    # sub-setting  
    fits = fits_all[[m]][rows,]
    
    # CFI
    cfi<- cfi.fp_2(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp_2(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp_2(fits[,"srmr"])
    
    # LOGL
    logl <- logl.fp_2(fits[,"logl"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi, rmsea, srmr, logl)
    
    # saving results
    results[m,] <- res
    
  }
  
  # output
  return(results)
}  
#-------------------------------------------------------------------------------

# Model 3
# function to get prop CFI >= .35
cfi.fp_3 <- function(fits, cutoff = .35) {
  
  
  return((sum(fits >= cutoff, na.rm = TRUE) / length(fits)))
}


# function to get prop RMSEA < .17
rmsea.fp_3 <- function(fits, cutoff = .17) {
  
  return((sum(fits < cutoff, na.rm = TRUE) / length(fits)))
}


# function to get prop SRMR < .065
srmr.fp_3 <- function(fits, cutoff = .065) {

  
  return((sum(fits < cutoff, na.rm = TRUE) / length(fits)))
}

# function to get the log likelihood with > -50,076------------------------------

logl.fp_3 <- function(fits, cutoff = -50076){
  
  return((sum(fits >= cutoff, na.rm = TRUE) / length(fits)))
  
}


# functions to sample X fit values and compute FP
samfits_3 <- function(fits_all, size = 100) {
  
  rows <- sample(1:nrow(fits_all[[1]]), size, replace = FALSE)
  
  # container
  results <- matrix(NA, length(fits_all), 5)
  
  # renaming columns
  colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl')
  
  for(m in seq_along(fits_all)){
    # sub-setting  
    fits = fits_all[[m]][rows,]
    
    # CFI
    cfi<- cfi.fp_3(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp_3(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp_3(fits[,"srmr"])
    
    # LOGL
    logl <- logl.fp_3(fits[,"logl"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi, rmsea, srmr, logl)
    
    # saving results
    results[m,] <- res
  }
  
  # output
  return(results)
}  



#-------------------------------------------------------------------------------
# Get FP for different space sizes ---------------------------------------------

# Model 1
# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits(model_1, size = size)))
  tmp <- (do.call(rbind, tmp))
  tmp <- cbind(tmp, rep(sizes, each = 16))
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl', 'size')
  
  if(i == 1) {
    fp.size.results_1 <- tmp
  } else {
    fp.size.results_1 <- rbind(fp.size.results_1, tmp)
  }
}


# Grid of conditions that we can use for post-processing ----
study3_conditions <- expand.grid(model = 1:3, 
                                 bounds = c("none", "pos.var"),
                                 ident = c("refmarker", "unitvar"),
                                 mlit = c(150, 1500),
                                 conv = c("F", "T"))

# Individual study conditions
study_con_model_1 <- study3_conditions %>%
  filter(model == 1)

# merging the data frame with the conditions
study_con_model_1$condition <- 1:16

res_pos_model_1 <- merge(fp.size.results_1, study_con_model_1)

res_pos_model_1$model <- factor(res_pos_model_1$model)

# saving the results
saveRDS(fp.size.results_1, file = here::here("Online-Supplement","Data-Files", "study3_poscor_model_1_true_dataspacesize.RDS"))

saveRDS(res_pos_model_1, file = here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_1_true_dataspacesize.RDS"))


# Model 2

# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits_2(model_2, size = size)))
  tmp <- (do.call(rbind, tmp))
  tmp <- cbind(tmp, rep(sizes, each = 16))
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl', 'size')
  
  if(i == 1) {
    fp.size.results_2 <- tmp
  } else {
    fp.size.results_2 <- rbind(fp.size.results_2, tmp)
  }
}


# Individual study conditions
study_con_model_2 <- study3_conditions %>%
  filter(model == 2)

# merging the data frame with the conditions
study_con_model_2$condition <- 1:16

res_pos_model_2 <- merge(fp.size.results_2, study_con_model_2)

res_pos_model_2$model <- factor(res_pos_model_2$model)

# saving the results
saveRDS(fp.size.results_2, file = here::here("Online-Supplement","Data-Files", "study3_poscor_model_2_true_dataspacesize.RDS"))

saveRDS(res_pos_model_2, file = here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_2_true_dataspacesize.RDS"))


# Model 3

# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits_3(model_3, size = size)))
  tmp <- (do.call(rbind, tmp))
  tmp <- cbind(tmp, rep(sizes, each = 16))
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr','logl', 'size')
  
  if(i == 1) {
    fp.size.results_3 <- tmp
  } else {
    fp.size.results_3 <- rbind(fp.size.results_3, tmp)
  }
}

# Individual study conditions
study_con_model_3 <- study3_conditions %>%
  filter(model == 3)

# merging the data frame with the conditions
study_con_model_3$condition <- 1:16

res_pos_model_3 <- merge(fp.size.results_3, study_con_model_3)

res_pos_model_3$model <- factor(res_pos_model_3$model)

# saving the results
saveRDS(fp.size.results_3, file = here::here("Online-Supplement","Data-Files", "study3_poscor_model_3_true_dataspacesize.RDS"))

saveRDS(res_pos_model_3, file = here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_3_true_dataspacesize.RDS"))

#-------------------------------------------------------------------------------
# All Correlations
#-------------------------------------------------------------------------------


# subsets
model_1 <- lapply(study3_all_fits, function(x) x %>% filter(model == 1))

model_2 <- lapply(study3_all_fits, function(x) x %>% filter(model == 2))

model_3 <- lapply(study3_all_fits, function(x) x %>% filter(model == 3))

# removing the empty data frames
model_1 <- model_1[sapply(model_1,nrow) > 0]

model_2 <- model_2[sapply(model_2,nrow) > 0]

model_3 <- model_3[sapply(model_3,nrow) > 0]


# Model 1
# function to get prop CFI >= .51
cfi.fp_1 <- function(fits, cutoff = .51) {
    
    return((sum(fits >= cutoff, na.rm = TRUE) / length(fits)))
}


# function to get prop RMSEA < .41
rmsea.fp_1 <- function(fits, cutoff = .41) {
  
  return((sum(fits < cutoff, na.rm = TRUE) / length(fits)))
}


# function to get prop SRMR < .17
srmr.fp_1 <- function(fits, cutoff = .17) {
  
  return((sum(fits < cutoff, na.rm = TRUE) / length(fits)))
}

# function to get the log likelihood with > -11,460---------------------------

logl.fp_1 <- function(fits, cutoff = -11460){
  
  return((sum(fits >= cutoff, na.rm = TRUE) / length(fits)))
  
}


# functions to sample X fit values and compute FP
samfits <- function(fits_all, size = 100) {
  
  rows <- sample(1:nrow(fits_all[[1]]), size, replace = FALSE)
  
  # container
  results <- matrix(NA, length(fits_all), 5)
  
  # renaming columns
  colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl')
  
  for(m in seq_along(fits_all)){
    
    # sub-setting  
    fits = fits_all[[m]][rows,]
    
    # CFI
    cfi<- cfi.fp_1(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp_1(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp_1(fits[,"srmr"])
    
    # LOGL
    logl <- logl.fp_1(fits[,"logl"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi, rmsea, srmr, logl)
    
    # saving results
    results[m,] <- res
    
  }
  
  # output
  return(results)
  
  
}
#-------------------------------------------------------------------------------
# Model 2
# function to get prop CFI >= .431
cfi.fp_2 <- function(fits, cutoff = .431) {
  
  return((sum(fits >= cutoff, na.rm = TRUE) / length(fits)))
}

# function to get prop RMSEA < .249
rmsea.fp_2 <- function(fits, cutoff = .249) {
  
  return((sum(fits < cutoff, na.rm = TRUE) / length(fits)))
}


# function to get prop SRMR < .09
srmr.fp_2 <- function(fits, cutoff = .084) {
  
  return((sum(fits < cutoff, na.rm = TRUE) / length(fits)))
}

# function to get the log likelihood with > -36,999------------------------------

logl.fp_2 <- function(fits, cutoff = -36999){
  
  return((sum(fits >= cutoff, na.rm = TRUE) / length(fits)))
  
}


# functions to sample X fit values and compute FP
samfits_2 <- function(fits_all, size = 100) {
  
  rows <- sample(1:nrow(fits_all[[1]]), size, replace = FALSE)
  
  # container
  results <- matrix(NA, length(fits_all), 5)
  
  # renaming columns
  colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl')
  
  for(m in seq_along(fits_all)){
    # sub-setting  
    fits = fits_all[[m]][rows,]
    
    # CFI
    cfi<- cfi.fp_2(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp_2(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp_2(fits[,"srmr"])
    
    # LOGL
    logl <- logl.fp_2(fits[,"logl"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi, rmsea, srmr, logl)
    
    # saving results
    results[m,] <- res
    
  }
  
  # output
  return(results)
}  
#-------------------------------------------------------------------------------

# Model 3
# function to get prop CFI >= .35
cfi.fp_3 <- function(fits, cutoff = .35) {
  
  return((sum(fits >= cutoff, na.rm = TRUE) / length(fits)))
}


# function to get prop RMSEA < .17
rmsea.fp_3 <- function(fits, cutoff = .17) {
  
  
  return((sum(fits < cutoff, na.rm = TRUE) / length(fits)))
}


# function to get prop SRMR < .065
srmr.fp_3 <- function(fits, cutoff = .065) {
  
  
  return((sum(fits < cutoff, na.rm = TRUE) / length(fits)))
}

# function to get the log likelihood with > -62,493------------------------------

logl.fp_3 <- function(fits, cutoff = -62493){
  
  return((sum(fits >= cutoff,  na.rm = TRUE) / length(fits)))
  
}


# functions to sample X fit values and compute FP
samfits_3 <- function(fits_all, size = 100) {
  
  rows <- sample(1:nrow(fits_all[[1]]), size, replace = FALSE)
  
  # container
  results <- matrix(NA, length(fits_all), 5)
  
  # renaming columns
  colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl')
  
  for(m in seq_along(fits_all)){
    # sub-setting  
    fits = fits_all[[m]][rows,]
    
    # CFI
    cfi<- cfi.fp_3(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp_3(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp_3(fits[,"srmr"])
    
    # LOGL
    logl <- logl.fp_3(fits[,"logl"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi, rmsea, srmr, logl)
    
    # saving results
    results[m,] <- res
    
  }
  
  # output
  return(results)
}  


#-------------------------------------------------------------------------------
## Get FP for different space sizes ----

# Model 1
# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits(model_1, size = size)))
  tmp <- (do.call(rbind, tmp))
  tmp <- cbind(tmp, rep(sizes, each = 16))
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl', 'size')
  
  if(i == 1) {
    fp.size.results_1 <- tmp
  } else {
    fp.size.results_1 <- rbind(fp.size.results_1, tmp)
  }
}


# Grid of conditions that we can use for post-processing ----
study3_conditions <- expand.grid(model = 1:3, 
                                 bounds = c("none", "pos.var"),
                                 ident = c("refmarker", "unitvar"),
                                 mlit = c(150, 1500),
                                 conv = c("F", "T"))

# Individual study conditions
study_con_model_1 <- study3_conditions %>%
  filter(model == 1)

# merging the data frame with the conditions
study_con_model_1$condition <- 1:16

res_all_model_1 <- merge(fp.size.results_1, study_con_model_1)

res_all_model_1$model <- factor(res_all_model_1$model)

# saving the results
saveRDS(fp.size.results_1, file = here::here("Online-Supplement","Data-Files", "study3_allcor_model_1_true_dataspacesize.RDS"))

saveRDS(res_all_model_1, file = here::here("Online-Supplement","Data-Files", "results_allcor_study_3_model_1_true_dataspacesize.RDS"))



# Model 2

# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits_2(model_2, size = size)))
  tmp <- (do.call(rbind, tmp))
  tmp <- cbind(tmp, rep(sizes, each = 16))
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl', 'size')
  
  if(i == 1) {
    fp.size.results_2 <- tmp
  } else {
    fp.size.results_2 <- rbind(fp.size.results_2, tmp)
  }
}


# Individual study conditions
study_con_model_2 <- study3_conditions %>%
  filter(model == 2)

# merging the data frame with the conditions
study_con_model_2$condition <- 1:16

res_all_model_2 <- merge(fp.size.results_2, study_con_model_2)

res_all_model_2$model <- factor(res_all_model_2$model)

# saving the results
saveRDS(fp.size.results_2, file = here::here("Online-Supplement","Data-Files", "study3_allcor_model_2_true_dataspacesize.RDS"))

saveRDS(res_all_model_2, file = here::here("Online-Supplement","Data-Files", "results_allcor_study_3_model_2_true_dataspacesize.RDS"))


# Model 3

# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits_3(model_3, size = size)))
  tmp <- (do.call(rbind, tmp))
  tmp <- cbind(tmp, rep(sizes, each = 16))
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr','logl', 'size')
  
  if(i == 1) {
    fp.size.results_3 <- tmp
  } else {
    fp.size.results_3 <- rbind(fp.size.results_3, tmp)
  }
}

# Individual study conditions
study_con_model_3 <- study3_conditions %>%
  filter(model == 3)

# merging the data frame with the conditions
study_con_model_3$condition <- 1:16

res_all_model_3 <- merge(fp.size.results_3, study_con_model_3)

res_all_model_3$model <- factor(res_all_model_3$model)

# saving the results
saveRDS(fp.size.results_3, file = here::here("Online-Supplement","Data-Files", "study3_allcor_model_3_true_dataspacesize.RDS"))

saveRDS(res_all_model_3, file = here::here("Online-Supplement","Data-Files", "results_allcor_study_3_model_3_true_dataspacesize.RDS"))

# clean work space
rm(list = ls())

#-------------------------------------------------------------------------------
# Analysis:
# loading the datasets
res_pos_model_1 <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_1_true_dataspacesize.RDS"))
res_all_model_1 <- readRDS(here::here("Online-Supplement","Data-Files","results_allcor_study_3_model_1_true_dataspacesize.RDS"))

res_pos_model_2 <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_2_true_dataspacesize.RDS"))
res_all_model_2 <- readRDS(here::here("Online-Supplement","Data-Files","results_allcor_study_3_model_2_true_dataspacesize.RDS"))

res_pos_model_3 <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_3_true_dataspacesize.RDS"))
res_all_model_3 <- readRDS(here::here("Online-Supplement","Data-Files","results_allcor_study_3_model_3_true_dataspacesize.RDS"))
#-------------------------------------------------------------------------------


# new variable in each data frame-----------------------------------------------
# Model 1
results_full_1 <- bind_rows(
  
  res_pos_model_1 %>% mutate(cor = 'pos') %>% mutate(cor = factor(cor)),
  
  res_all_model_1 %>% mutate(cor = 'all') %>% mutate(cor = factor(cor))
)

# Model 2
results_full_2 <- bind_rows(
  
  res_pos_model_2 %>% mutate(cor = 'pos') %>% mutate(cor = factor(cor)),
  
  res_all_model_2 %>% mutate(cor = 'all') %>% mutate(cor = factor(cor))
)

# Model 3
results_full_3 <- bind_rows(
  
  res_pos_model_3 %>% mutate(cor = 'pos') %>% mutate(cor = factor(cor)),
  
  res_all_model_3 %>% mutate(cor = 'all') %>% mutate(cor = factor(cor))
)

# checking for any missing data so we can confirm that our merging strategy worked
apply(results_full_1, 2, function(x) 100* mean(is.na(x)))

apply(results_full_2, 2, function(x) 100* mean(is.na(x)))

apply(results_full_3, 2, function(x) 100* mean(is.na(x)))

# Fitting propensity------------------------------------------------------------
results_full_1 %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(size == 1e4 & cor == "all") %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size, index, bounds,ident, mlit, conv) %>%
  summarize(m = mean(value),
            se = sd(value), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>%
  pivot_wider(names_from = c("cor", "model","conv", "mlit"), values_from = result) -> table_111

results_full_2 %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(size == 1e4 & cor == "all") %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size, index, bounds,ident, mlit, conv) %>%
  summarize(m = mean(value),
            se = sd(value), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>%
  pivot_wider(names_from = c("cor", "model","conv", "mlit"), values_from = result) -> table_222

results_full_3 %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(size == 1e4 & cor == "all") %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size, index, bounds,ident, mlit, conv) %>%
  summarize(m = mean(value),
            se = sd(value), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>%
  pivot_wider(names_from = c("cor", "model","conv", "mlit"), values_from = result) -> table_333
#-------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
results_full_1 %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  filter(size <= 10000) %>%
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c('Model 1C', 'Model 2C', 'Model 3C')),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr", "logl"), labels = c("CFI", "RMSEA", "SRMR", "LL")), 
         
         bounds = factor(bounds, levels = c('none', 'pos.var'), labels = c('Unrestricted', 'Positive')),
         
         ident = factor(ident, levels = c('refmarker', 'unitvar'), labels = c('Reference', 'Unit Variance')),
         
         mlit = factor(mlit, levels = c('150', '1500')),
         
         conv = factor(conv, levels = c('F', 'T'), labels = c('Default', 'Forced')),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == 'LL' & cor == 'Positive Subspace') %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  .width = 0.95, linewidth = 0.4, alpha = 0.9) +
  facet_grid(rows = vars(mlit),
             cols = vars(ident, conv),
             scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (proportion  > cutoff)",
       fill = 'Model', fill_ramp = 'Confidence Level', linetype = 'Model') +
  scale_fill_manual(name = 'Model',
                    values = c('Model 1C' = "#276419",
                               'Model 2C' = "#8e0152", 
                               'Model 3C' = "#b2182b")) +
  scale_linetype_manual(
    name = "Model",
    values = c("Model 1C" = "solid", "Model 2C" = "dashed", "Model 3C" = "dotted")
  ) +
  scale_x_continuous(breaks = seq(0, 1e4, 2e3),
                     labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                     sec.axis = sec_axis(~ . , name = "Positive Subspace", breaks = NULL, labels = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                     sec.axis = sec_axis(~ . , name = "Maximum Iterations", breaks = NULL, labels = NULL)) -> plot_1

# saving the graph
ggsave(plot = plot_1, path = here::here("Online-Supplement","Plots", "Study Three"), 
filename = 'Fit-positive-LL-dataspace-study-3-model-1-listwise.png', width = 10, height = 8,
units = 'in')


#---
results_full_1 %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  filter(size <= 10000) %>%
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c('Model 1C', 'Model 2C', 'Model 3C')),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr", "logl"), labels = c("CFI", "RMSEA", "SRMR", "LL")), 
         
         bounds = factor(bounds, levels = c('none', 'pos.var'), labels = c('Unrestricted', 'Positive')),
         
         ident = factor(ident, levels = c('refmarker', 'unitvar'), labels = c('Reference', 'Unit Variance')),
         
         mlit = factor(mlit, levels = c('150', '1500')),
         
         conv = factor(conv, levels = c('F', 'T'), labels = c('Default', 'Forced')),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == 'LL' & cor == 'Complete Space') %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  .width = 0.95, linewidth = 0.4, alpha = 0.9) +
  facet_grid(rows = vars(mlit),
             cols = vars(ident, conv),
             scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (proportion  > cutoff)",
       fill = 'Model', fill_ramp = 'Confidence Level', linetype = 'Model') +
  scale_fill_manual(name = 'Model',
                    values = c('Model 1C' = "#276419",
                               'Model 2C' = "#8e0152", 
                               'Model 3C' = "#b2182b")) +
  scale_linetype_manual(
    name = "Model",
    values = c("Model 1C" = "solid", "Model 2C" = "dashed", "Model 3C" = "dotted")
  ) +
  scale_x_continuous(breaks = seq(0, 1e4, 2e3),
                     labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                     sec.axis = sec_axis(~ . , name = "Complete Space", breaks = NULL, labels = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                     sec.axis = sec_axis(~ . , name = "Maximum Iterations", breaks = NULL, labels = NULL)) -> plot_11

# saving the graph
ggsave(plot = plot_11, path = here::here("Online-Supplement","Plots", "Study Three"), 
filename = 'Fit-complete-LL-dataspace-study-3-model-1-listwise.png', width = 10, height = 8,
units = 'in')

#-----------------------------------------------------------------
results_full_2 %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  filter(size <= 10000) %>%
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c('Model 1C', 'Model 2C', 'Model 3C')),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr", "logl"), labels = c("CFI", "RMSEA", "SRMR", "LL")), 
         
         bounds = factor(bounds, levels = c('none', 'pos.var'), labels = c('Unrestricted', 'Positive')),
         
         ident = factor(ident, levels = c('refmarker', 'unitvar'), labels = c('Reference', 'Unit Variance')),
         
         mlit = factor(mlit, levels = c('150', '1500')),
         
         conv = factor(conv, levels = c('F', 'T'), labels = c('Default', 'Forced')),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == 'LL' & cor == 'Positive Subspace') %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  .width = 0.95, linewidth = 0.4, alpha = 0.9) +
  facet_grid(rows = vars(mlit),
             cols = vars(ident, conv),
             scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (proportion  > cutoff)",
       fill = 'Model', fill_ramp = 'Confidence Level', linetype = 'Model') +
  scale_fill_manual(name = 'Model',
                    values = c('Model 1C' = "#276419",
                               'Model 2C' = "#8e0152", 
                               'Model 3C' = "#b2182b")) +
  scale_linetype_manual(
    name = "Model",
    values = c("Model 1C" = "solid", "Model 2C" = "dashed", "Model 3C" = "dotted")
  ) +
  scale_x_continuous(breaks = seq(0, 1e4, 2e3),
                     labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                     sec.axis = sec_axis(~ . , name = "Positive Subspace", breaks = NULL, labels = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                     sec.axis = sec_axis(~ . , name = "Maximum Iterations", breaks = NULL, labels = NULL)) -> plot_2

# saving the graph
ggsave(plot = plot_2, path = here::here("Online-Supplement","Plots", "Study Three"),
filename = 'Fit-positive-LL-dataspace-study-3-model-2-listwise.png', width = 10, height = 8,
units = 'in')


#---
results_full_2 %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  filter(size <= 10000) %>%
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c('Model 1C', 'Model 2C', 'Model 3C')),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr", "logl"), labels = c("CFI", "RMSEA", "SRMR", "LL")), 
         
         bounds = factor(bounds, levels = c('none', 'pos.var'), labels = c('Unrestricted', 'Positive')),
         
         ident = factor(ident, levels = c('refmarker', 'unitvar'), labels = c('Reference', 'Unit Variance')),
         
         mlit = factor(mlit, levels = c('150', '1500')),
         
         conv = factor(conv, levels = c('F', 'T'), labels = c('Default', 'Forced')),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == 'LL' & cor == 'Complete Space') %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  .width = 0.95, linewidth = 0.4, alpha = 0.9) +
  facet_grid(rows = vars(mlit),
             cols = vars(ident, conv),
             scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (proportion  > cutoff)",
       fill = 'Model', fill_ramp = 'Confidence Level', linetype = 'Model') +
  scale_fill_manual(name = 'Model',
                    values = c('Model 1C' = "#276419",
                               'Model 2C' = "#8e0152", 
                               'Model 3C' = "#b2182b")) +
  scale_linetype_manual(
    name = "Model",
    values = c("Model 1C" = "solid", "Model 2C" = "dashed", "Model 3C" = "dotted")
  ) +
  scale_x_continuous(breaks = seq(0, 1e4, 2e3),
                     labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                     sec.axis = sec_axis(~ . , name = "Complete Space", breaks = NULL, labels = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                     sec.axis = sec_axis(~ . , name = "Maximum Iterations", breaks = NULL, labels = NULL)) -> plot_22

# saving the graph
ggsave(plot = plot_22, path = here::here("Online-Supplement","Plots", "Study Three"), 
filename = 'Fit-complete-LL-dataspace-study-3-model-2-listwise.png', width = 10, height = 8,
units = 'in')

#--------------------------------------------------------------------------------
results_full_3 %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  filter(size <= 10000) %>%
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c('Model 1C', 'Model 2C', 'Model 3C')),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr", "logl"), labels = c("CFI", "RMSEA", "SRMR", "LL")), 
         
         bounds = factor(bounds, levels = c('none', 'pos.var'), labels = c('Unrestricted', 'Positive')),
         
         ident = factor(ident, levels = c('refmarker', 'unitvar'), labels = c('Reference', 'Unit Variance')),
         
         mlit = factor(mlit, levels = c('150', '1500')),
         
         conv = factor(conv, levels = c('F', 'T'), labels = c('Default', 'Forced')),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == 'LL' & cor == 'Positive Subspace') %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  .width = 0.95, linewidth = 0.4, alpha = 0.9) +
  facet_grid(rows = vars(mlit),
             cols = vars(ident, conv),
             scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (proportion  > cutoff)",
       fill = 'Model', fill_ramp = 'Confidence Level', linetype = 'Model') +
  scale_fill_manual(name = 'Model',
                    values = c('Model 1C' = "#276419",
                               'Model 2C' = "#8e0152", 
                               'Model 3C' = "#b2182b")) +
  scale_linetype_manual(
    name = "Model",
    values = c("Model 1C" = "solid", "Model 2C" = "dashed", "Model 3C" = "dotted")
  ) +
  scale_x_continuous(breaks = seq(0, 1e4, 2e3),
                     labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                     sec.axis = sec_axis(~ . , name = "Positive Subspace", breaks = NULL, labels = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                     sec.axis = sec_axis(~ ., name = "Maximum Iterations", breaks = NULL, labels = NULL)) -> plot_3

# saving the graph
ggsave(plot = plot_3, path = here::here("Online-Supplement","Plots", "Study Three"),
filename = 'Fit-positive-LL-dataspace-study-3-model-3-listwise.png', width = 10, height = 8,
units = 'in')


#---
results_full_3 %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  filter(size <= 10000) %>%
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c('Model 1C', 'Model 2C', 'Model 3C')),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr", "logl"), labels = c("CFI", "RMSEA", "SRMR", "LL")), 
         
         bounds = factor(bounds, levels = c('none', 'pos.var'), labels = c('Unrestricted', 'Positive')),
         
         ident = factor(ident, levels = c('refmarker', 'unitvar'), labels = c('Reference', 'Unit Variance')),
         
         mlit = factor(mlit, levels = c('150', '1500')),
         
         conv = factor(conv, levels = c('F', 'T'), labels = c('Default', 'Forced')),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == 'LL' & cor == 'Complete Space') %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  .width = 0.95, linewidth = 0.4, alpha = 0.9) +
  facet_grid(rows = vars(mlit),
             cols = vars(ident, conv),
             scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (proportion  > cutoff)",
       fill = 'Model', fill_ramp = 'Confidence Level', linetype = 'Model') +
  scale_fill_manual(name = 'Model',
                    values = c('Model 1C' = "#276419",
                               'Model 2C' = "#8e0152", 
                               'Model 3C' = "#b2182b")) +
  scale_linetype_manual(
    name = "Model",
    values = c("Model 1C" = "solid", "Model 2C" = "dashed", "Model 3C" = "dotted")
  ) +
  scale_x_continuous(breaks = seq(0, 1e4, 2e3),
                     labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                     sec.axis = sec_axis(~ . , name = "Complete Space", breaks = NULL, labels = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                     sec.axis = sec_axis(~ . , name = "Maximum Iterations", breaks = NULL, labels = NULL)) -> plot_33

# saving the graph
ggsave(plot = plot_33, path = here::here("Online-Supplement","Plots", "Study Three"), 
filename = 'Fit-all-LL-dataspace-study-3-model-3-listwise.png', width = 10, height = 8,
units = 'in')
#-------------------------------------------------------------------------------

## Total fitting propensity for each fit measure
results_full_1 %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(size == 1e4 & cor == "all") %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size, index, bounds,ident, mlit, conv) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop")%>% 
  pivot_wider(names_from = c("cor", "model","conv", "mlit"), values_from = total) -> total_111

results_full_2 %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(size == 1e4 & cor == "all") %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size, index, bounds,ident, mlit, conv) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop")%>% 
  pivot_wider(names_from = c("cor", "model","conv", "mlit"), values_from = total) -> total_222

results_full_3 %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(size == 1e4 & cor == "all") %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size, index, bounds,ident, mlit, conv) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop")%>% 
  pivot_wider(names_from = c("cor", "model","conv", "mlit"), values_from = total) -> total_333
#-------------------------------------------------------------------------------
# Saving the results above
saveRDS(table_111, file = here::here("Online-Supplement", "Data-Files", "table_111.RDS"))

saveRDS(table_222, file = here::here("Online-Supplement", "Data-Files", "table_222.RDS"))

saveRDS(table_333, file = here::here("Online-Supplement", "Data-Files", "table_333.RDS"))

saveRDS(total_111, file = here::here("Online-Supplement","Data-Files", "total_111.RDS"))

saveRDS(total_222, file =  here::here("Online-Supplement","Data-Files", "total_222.RDS"))

saveRDS(total_333, file =  here::here("Online-Supplement","Data-Files", "total_333.RDS"))
#-------------------------------------------------------------------------------




