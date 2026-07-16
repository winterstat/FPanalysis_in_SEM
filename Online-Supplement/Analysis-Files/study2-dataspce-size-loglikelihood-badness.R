#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 09/09/2025
#-------------------------------------------------------------------------------
# Relationship between the dataspace and the fitting propensity of the models
# in study 2
#-------------------------------------------------------------------------------

## Load packages ---------------------------------------------------------------
pkgs <- c('tidyverse', 'ggdist', 'ggthemes', 'effectsize', 'dplyr', 'ggplot2')

# check if packages are available
#rlang::check_installed(pkgs)

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


# Load full data
load(here::here("Online-Supplement","Data-Files", "output_poscor_2_all_conditions.RData"))

# positive correlations
study2_pos_fits <- output_conditions

load(here::here("Online-Supplement","Data-Files", "output_all_2_conditions.RData"))

study2_all_fits <- output_conditions

# removing the data
rm(output_conditions)

# Select range of data space sizes ----
sizes <- c(50, 100, 500, 1000, 5000, 10000, 50000, 100000, 150000, 175000, 200000, 500000, 750000)

# Set up functions to compute FP (proportion of fit index values above/below cutoff) ----

# We replace the non-convergence with some badness of fit value

#-------------------------------------------------------------------------------
# Positive correlations
#-------------------------------------------------------------------------------

# function to get prop CFI >= .45-----------------------------------------------
cfi.fp <- function(fits, cutoff = .45) {
  
  fits[is.na(fits)] <- 0
  
  return((sum(fits >= cutoff) / length(fits)))
}

# function to get prop RMSEA < .27----------------------------------------------
rmsea.fp <- function(fits, cutoff = .27) {
  
  fits[is.na(fits)] <- 1
  
  return((sum(fits < cutoff) / length(fits)))
}

# function to get prop SRMR < .09-----------------------------------------------
srmr.fp <- function(fits, cutoff = .09) {
  
  fits[is.na(fits)] <- 1
  
  return((sum(fits < cutoff) / length(fits)))
}

# function to get the log likelihood with > -27103------------------------------

logl.fp <- function(fits, cutoff = -27103){
  
  fits[is.na(fits)] <- -35000
  
  return((sum(fits >= cutoff) / length(fits)))
  
}

# functions to sample X fit values and compute FP-------------------------------
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
    cfi<- cfi.fp(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp(fits[,"srmr"])
    
    # LOGL
    logl <- logl.fp(fits[,"logl"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi, rmsea, srmr, logl)
    
    # saving results
    results[m,] <- res
    
  }
  
  # output
  return(results)
  
  
}

# Get FP for different space sizes --------------------------------------------

# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits(study2_pos_fits, size = size)))
  
  tmp <- (do.call(rbind, tmp))
  
  tmp <- cbind(tmp, rep(sizes, each = 48))
  
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl', 'size')
  
  if(i == 1) {
    fp.size.results <- tmp
  } else {
    fp.size.results <- rbind(fp.size.results, tmp)
  }
}

#-------------------------------------------------------------------------------
# Grid of conditions that we can use for post-processing ----
study2_conditions <- expand.grid(bounds = c("none", "pos.var"),
                                 ident = c("refmarker", "unitvar"),
                                 mlit = c('150', '1500'),
                                 conv = c("F", "T"),
                                 model = 1:3)

# merging the data frame with the conditions
study2_conditions$condition <- 1:48

res_pos <- merge(fp.size.results, study2_conditions)

res_pos$model <- factor(res_pos$model)

# saving the results
saveRDS(fp.size.results, file = here::here("Online-Supplement","Data-Files", "study2_poscor_badness_dataspacesize.RDS"))

saveRDS(res_pos, file = here::here("Online-Supplement","Data-Files", "results_poscor_badness_study_2_dataspacesize.RDS"))
#---------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# All correlations
#-------------------------------------------------------------------------------

# function to get prop CFI >= .20-----------------------------------------------
cfi.fp_1 <- function(fits, cutoff = .20) {
  
  fits[is.na(fits)] <- 0
  
  return((sum(fits >= cutoff) / length(fits)))
}

# function to get prop RMSEA < .27----------------------------------------------
rmsea.fp_1 <- function(fits, cutoff = .27) {
  
  fits[is.na(fits)] <- 1
  
  return((sum(fits < cutoff) / length(fits)))
}

# function to get prop SRMR < .17-----------------------------------------------
srmr.fp_1 <- function(fits, cutoff = .17) {
  
  fits[is.na(fits)] <- 1
  
  return((sum(fits < cutoff) / length(fits)))
}

# function to get the log likelihood with > -32755------------------------------

logl.fp_1 <- function(fits, cutoff = -32755){
  
  fits[is.na(fits)] <- -35000
  
  return((sum(fits >= cutoff) / length(fits)))
  
}

# functions to sample X fit values and compute FP-------------------------------
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

# Get FP for different space sizes --------------------------------------------

# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits(study2_all_fits, size = size)))
  
  tmp <- (do.call(rbind, tmp))
  
  tmp <- cbind(tmp, rep(sizes, each = 48))
  
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr', 'logl', 'size')
  
  if(i == 1) {
    fp.size.results_1 <- tmp
  } else {
    fp.size.results_1 <- rbind(fp.size.results_1, tmp)
  }
}

#-------------------------------------------------------------------------------

# Grid of conditions that we can use for post-processing ----
study2_conditions <- expand.grid(bounds = c("none", "pos.var"),
                                 ident = c("refmarker", "unitvar"),
                                 mlit = c('150', '1500'),
                                 conv = c("F", "T"),
                                 model = 1:3)

# merging the data frame with the conditions
study2_conditions$condition <- 1:48

res_all <- merge(fp.size.results_1, study2_conditions)

res_all$model <- factor(res_all$model)


# saving the results
saveRDS(fp.size.results_1, file = here::here("Online-Supplement","Data-Files", "study2_all_badness_dataspacesize.RDS"))

saveRDS(res_all, file = here::here("Online-Supplement","Data-Files", "results_all_badness_study_2_dataspacesize.RDS"))
#---------------------------------------------------------------------------------

# deleting everything
rm(list = ls())

# reading and merging the datasets----------------------------------------------
res_pos <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_badness_study_2_dataspacesize.RDS"))

res_all <- readRDS(here::here("Online-Supplement","Data-Files", "results_all_badness_study_2_dataspacesize.RDS"))

# new variable in each data frame-----------------------------------------------
results_full <- bind_rows(
  
  res_pos %>% mutate(cor = 'pos') %>% mutate(cor = factor(cor)),
  
  res_all %>% mutate(cor = 'all') %>% mutate(cor = factor(cor))
)

# checking for any missing data so we can confirm that our merging strategy worked
apply(results_full, 2, function(x) 100* mean(is.na(x)))

# Fitting propensity------------------------------------------------------------
results_full %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(size == 1e4 & cor == "all") %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size, index, bounds,ident, mlit, conv) %>%
  summarize(m = mean(value),
            se = sd(value), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>%
  pivot_wider(names_from = c("cor", "model", "conv", "mlit"), values_from = result) -> table_1

# size index bounds  ident     all_1_F_150 all_1_T_150   all_1_F_1500  all_1_T_1500  all_2_F_150   all_2_T_150   all_2_F_1500  all_2_T_1500  all_3_F_150   all_3_T_150   all_3_F_1500  all_3_T_1500 
# <dbl> <chr> <fct>   <fct>     <chr>       <chr>         <chr>         <chr>         <chr>         <chr>         <chr>         <chr>         <chr>         <chr>         <chr>         <chr>        
# 1 10000 logl  none    refmarker 0.001 (0)   0.215 (0.004) 0.005 (0.001) 0.26 (0.005)  0.021 (0.002) 0.701 (0.004) 0.057 (0.002) 0.766 (0.004) 0.013 (0.001) 0.996 (0.001) 0.082 (0.003) 0.999 (0)    
# 2 10000 logl  none    unitvar   0.001 (0)   0.257 (0.005) 0.001 (0)     0.302 (0.005) 0.013 (0.001) 0.535 (0.005) 0.019 (0.001) 0.569 (0.004) 0.008 (0.001) 0.994 (0.001) 0.058 (0.002) 0.996 (0.001)
# 3 10000 logl  pos.var refmarker 0.001 (0)   0.001 (0)     0.001 (0)     0.001 (0)     0.041 (0.002) 0.324 (0.005) 0.084 (0.003) 0.366 (0.005) 0.07 (0.002)  0.9 (0.003)   0.336 (0.005) 0.97 (0.002) 
# 4 10000 logl  pos.var unitvar   0.001 (0)   0.001 (0)     0.001 (0)     0.001 (0)     0.025 (0.002) 0.178 (0.004) 0.029 (0.002) 0.186 (0.004) 0.193 (0.004) 0.92 (0.002)  0.291 (0.004) 0.921 (0.003)


#-------------------------------------------------------------------------------
## Total fitting propensity for each fit measure
results_full %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(size == 1e4 & cor == "all") %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size, index, bounds,ident, mlit, conv) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop")%>% 
  pivot_wider(names_from = c("cor", "model", "conv", "mlit"), values_from = total) -> total_1

# # A tibble: 4 × 16
# size index bounds  ident     all_1_F_150 all_1_T_150 all_1_F_1500 all_1_T_1500 all_2_F_150 all_2_T_150 all_2_F_1500 all_2_T_1500 all_3_F_150 all_3_T_150 all_3_F_1500 all_3_T_1500
# <dbl> <chr> <fct>   <fct>           <dbl>       <dbl>        <dbl>        <dbl>       <dbl>       <dbl>        <dbl>        <dbl>       <dbl>       <dbl>        <dbl>        <dbl>
# 1 10000 logl  none    refmarker     0.00103     0.215        0.00485      0.260        0.0209       0.701       0.0575        0.766     0.0126        0.996       0.0817        0.999
# 2 10000 logl  none    unitvar       0.00123     0.257        0.0013       0.302        0.0134       0.535       0.0192        0.569     0.00812       0.994       0.0581        0.996
# 3 10000 logl  pos.var refmarker     0.00111     0.00121      0.00112      0.00121      0.0405       0.324       0.0837        0.366     0.0701        0.900       0.336         0.970
# 4 10000 logl  pos.var unitvar       0.00124     0.00124      0.00124      0.00124      0.0246       0.178       0.0291        0.186     0.193         0.920       0.291         0.921
#-------------------------------------------------------------------------------


# figures------------------------------------------------------------------------

# setting the theme
my_theme <- theme_grey() + theme(panel.background = element_blank(),
        legend.position = 'bottom',
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), strip.text = element_text(face = 'bold', size = 6.0),
        strip.background = element_rect(colour = 'grey79'), axis.line = element_line(colour = 'grey79', linewidth = .35),
        axis.text = element_text(size = 6.0)) 

theme_set(my_theme)

# figures-----------------------------------------------------------------------
results_full %>%
  filter(size <= 1e4) %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c('Model 2A', 'Model 2B', 'Model 2C')),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr", "logl"), labels = c("CFI", "RMSEA", "SRMR", "LL")), 
         
         bounds = factor(bounds, levels = c('none', 'pos.var'), labels = c('Unrestricted', 'Positive')),
         
         ident = factor(ident, levels = c('refmarker', 'unitvar'), labels = c('Reference', 'Unit Variance')),
         
         mlit = factor(mlit, levels = c('150', '1500')),
         
         conv = factor(conv, levels = c("F", "T"), labels = c("Default", "Forced")),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == 'LL' & cor == 'Complete Space') %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  .width = 0.95, linewidth = 0.4, alpha = 0.9) +
  facet_grid(rows = vars(mlit),
             cols = vars(ident, conv), scales = "free")+
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (LL > -32755)",
       fill = 'Model', fill_ramp = 'Confidence Level', linetype = 'Model') +
  scale_fill_manual(name = 'Model',
                    values = c('Model 2A' = "#276419",
                               'Model 2B' = "#8e0152", 
                               'Model 2C' = "#b2182b")) +
  scale_linetype_manual(
    name = "Model",
    values = c("Model 2A" = "solid", "Model 2B" = "dashed", 'Model 2C' = 'dotted')
  ) +
  scale_x_continuous(breaks = seq(0, 1e4, 2e3),
                     labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                     sec.axis = sec_axis(~ . , name = "Complete Space", breaks = NULL, labels = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                     sec.axis = sec_axis(~ . , name = "Maximum Iterations", breaks = NULL, labels = NULL)) -> plot_1

# saving the graph
ggsave(plot = plot_1, path = here::here("Online-Supplement","Plots", "Study Two"), 
       filename = 'Fit-logl-complete-space-study-2.png', width = 10, height = 8,
       units = 'in')

#----
results_full %>%
  filter(size <= 1e4) %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c('Model 2A', 'Model 2B', 'Model 2C')),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr", "logl"), labels = c("CFI", "RMSEA", "SRMR", "LL")), 
         
         bounds = factor(bounds, levels = c('none', 'pos.var'), labels = c('Unrestricted', 'Positive')),
         
         ident = factor(ident, levels = c('refmarker', 'unitvar'), labels = c('Reference', 'Unit Variance')),
         
         mlit = factor(mlit, levels = c('150', '1500')),
         
         conv = factor(conv, levels = c("F", "T"), labels = c("Default", "Forced")),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == 'LL' & cor == 'Positive Subspace') %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  .width = 0.95, linewidth = 0.4, alpha = 0.9) +
  facet_grid(rows = vars(mlit),
             cols = vars(ident, conv), scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (LL > -32755)",
       fill = 'Model', fill_ramp = 'Confidence Level', linetype = 'Model') +
  scale_fill_manual(name = 'Model',
                    values = c('Model 2A' = "#276419",
                               'Model 2B' = "#8e0152", 
                               'Model 3C' = "#b2182b")) +
  scale_linetype_manual(
    name = "Model",
    values = c("Model 2A" = "solid", "Model 2B" = "dashed", 'Model 2C' = 'dotted')
  ) +
  scale_x_continuous(breaks = seq(0, 1e4, 2e3),
                     labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                     sec.axis = sec_axis(~ . , name = "Positive Subspace", breaks = NULL, labels = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                     sec.axis = sec_axis(~ . , name = "Maximum Iterations", breaks = NULL, labels = NULL)) -> plot_2

# saving the graph
ggsave(plot = plot_2, path = here::here("Online-Supplement","Plots", "Study Two"), 
       filename = 'Fit-logl-positive-subspace-study-2.png', width = 10, height = 8,
       units = 'in')
