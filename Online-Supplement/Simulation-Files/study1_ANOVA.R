#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 09/25/2025
#-------------------------------------------------------------------------------
# ANOVA Analysis for Study 1
#-------------------------------------------------------------------------------

# packages
pkgs <- c('ockhamSEM', 'parallel', 'doParallel', 'lavaan', 'ggthemes', 'tidyverse', 'ggdist', 'effectsize')

# check if this packages have been installed
rlang::check_required(pkgs)

# loading them
lapply(pkgs, library, character.only = TRUE)


# # create a local package library (run once)
# dir.create("C:/Rlibs", showWarnings = FALSE)
# 
# # put the new folder first in the library search path for this session
# .libPaths("C:/Rlibs")
# 
# # now install into the local folder
# install.packages("xfun", lib = .libPaths()[1], dependencies = TRUE)


## Load full data space(s) -----------------------------------------------------
study1_pos <- readRDS("results/study1_poscor.RDS")

study1_all <- readRDS("results/study1_allcor.RDS")

## Extract fit information (list with 1 element for each model)
study1_pos_fits <- study1_pos$fit_list

study1_all_fits <- study1_all$fit_list

study1_fits <- c(study1_pos_fits, study1_all_fits)

## Combine outputs in one data frame ----
for(i in 1:4) {
  tmp <- data.frame(study1_fits[[i]])
  
  tmp$model <- if_else(i %in% c(1,3), "A", "B")
  
  tmp$cor <- if_else(i < 3, "pos", "all")
  
  if(i == 1) {
    study1_total <- tmp
  } else {
    study1_total <- bind_rows(study1_total, tmp)
  }
}

# checking for any missing data so we can confirm that our merging strategy worked
apply(study1_total, 2, function(x) 100* mean(is.na(x)))

study1_total <- study1_total %>%
  mutate(
    model = factor(model, levels = c('A', 'B')),
    
    cor = factor(cor, levels = c('pos', 'all'))
  )

## ANOVA ----
aov_srmr <- aov(srmr ~ model + cor, data = study1_total)

# results
summary(aov_srmr)

# effect sizes
effectsize(aov_srmr)

aov_cfi <- aov(cfi ~ model + cor, data = study1_total)

# results
summary(aov_cfi)

# effect sizes
effectsize(aov_cfi)

aov_rmsea <- aov(rmsea ~ model + cor, data = study1_total)

# results
summary(aov_rmsea)

# effect sizes
effectsize(aov_rmsea)


aov_logl <- aov(logl ~ model + cor, data = study1_total)

# results
summary(aov_logl)

# effect sizes
effectsize(aov_logl)


# Sensitivity analysis----------------------------------------------------------
# using Blimp

# re-coding the the factors as 0 and 1
study1_total_coded <- study1_total %>%
  mutate(
    
    model = ifelse(model == 'A', 0, 1),
    
    cor = ifelse(cor == 'all', 0, 1)
  )

# package
library(rblimp)

# fitting a Bayesian imputation model
srmr_blimp <- rblimp(
  
  # data
  data = study1_total_coded,
  
  # factors
  ordinal = 'cor model',
  
  # model
  model = 'srmr ~ cor  model',
  
  # MCMC
  seed = 3242,
  
  # iterations
  iter = 5000,
  
  # chains
  chains = 6,
  
  # number of imputes
  nimps = 5
  
)


# loading the mitml package