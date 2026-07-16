#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 09/09/2025
#-------------------------------------------------------------------------------
# Relationship between the dataspace size and the fitting propensity of the models
# in study 2
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# List wise deletion
#--------------------------------------------------------------------------------

# Load packages ---------------------------------------------------------------
pkgs <- c('tidyverse', 'ggdist', 'ggthemes', 'effectsize', 'dplyr', 'ggplot2')

# check if packages are available
#rlang::check_installed(pkgs)

# loading packages
lapply(pkgs, library, character.only = TRUE)


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
#-------------------------------------------------------------------------------
# Positive correlations
#-------------------------------------------------------------------------------

# function to get prop CFI >= .45-----------------------------------------------
cfi.fp <- function(fits, cutoff = .45){
  return(sum(fits >= cutoff, na.rm = T)/sum(!is.na(fits)))
}

# function to get prop RMSEA < .27----------------------------------------------
rmsea.fp <- function(fits, cutoff = .27) {
  
  return(sum(fits < cutoff,na.rm = T)/sum(!is.na(fits)))
}

# function to get prop SRMR < .09-----------------------------------------------
srmr.fp <- function(fits, cutoff = .09) {
  
  return(sum(fits < cutoff, na.rm = T)/sum(!is.na(fits)))
}

# function to get the log likelihood with > -27103------------------------------
logl.fp <- function(fits, cutoff = -27103){
  
  return(sum(fits >= cutoff, na.rm = T)/sum(!is.na(fits)))
  
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

res_pos_true <- merge(fp.size.results, study2_conditions)

res_pos_true$model <- factor(res_pos_true$model)

# saving the results
saveRDS(fp.size.results, file = here::here("Online-Supplement","Data-Files", "study2_poscor_true_dataspacesize.RDS"))

saveRDS(res_pos_true, file = here::here("Online-Supplement","Data-Files", "results_poscor_true_study_2_dataspacesize.RDS"))
#---------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# All correlations
#-------------------------------------------------------------------------------

# function to get prop CFI >= .20-----------------------------------------------
cfi.fp_1 <- function(fits, cutoff = .20) {
  
  return(sum(fits >= cutoff, na.rm = T)/sum(!is.na(fits)))
}

# function to get prop RMSEA < .27----------------------------------------------
rmsea.fp_1 <- function(fits, cutoff = .27) {
  
  return(sum(fits < cutoff, na.rm = T)/sum(!is.na(fits)))
}

# function to get prop SRMR < .17-----------------------------------------------
srmr.fp_1 <- function(fits, cutoff = .17) {
  
  return(sum(fits < cutoff, na.rm = T)/sum(!is.na(fits)))
}

# function to get the log likelihood with > -32755------------------------------

logl.fp_1 <- function(fits, cutoff = -32755){
  
  return(sum(fits >= cutoff, na.rm = T)/sum(!is.na(fits)))
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

res_all_true <- merge(fp.size.results_1, study2_conditions)

res_all_true$model <- factor(res_all_true$model)


# saving the results
saveRDS(fp.size.results_1, file = here::here("Online-Supplement","Data-Files", "study2_all_true_dataspacesize.RDS"))

saveRDS(res_all_true, file = here::here("Online-Supplement","Data-Files", "results_all_true_study_2_dataspacesize.RDS"))
#---------------------------------------------------------------------------------

# deleting everything
rm(list = ls())

# reading and merging the datasets----------------------------------------------
res_pos_true <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_true_study_2_dataspacesize.RDS"))

res_all_true <- readRDS(here::here("Online-Supplement","Data-Files", "results_all_true_study_2_dataspacesize.RDS"))

# new variable in each data frame-----------------------------------------------
results_full <- bind_rows(
  
  res_pos_true %>% mutate(cor = 'pos') %>% mutate(cor = factor(cor)),
  
  res_all_true %>% mutate(cor = 'all') %>% mutate(cor = factor(cor))
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
            se = sd(value), .groups = "drop") %>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")")) %>%
  select(-m, -se) %>%
  pivot_wider(names_from = c("cor", "model", "conv", "mlit"), values_from = result) -> table_11

#-------------------------------------------------------------------------------
## Total fitting propensity for each fit measure
results_full %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(size == 1e4 & cor == "all") %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size, index, bounds,ident, mlit, conv) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop")%>% 
  pivot_wider(names_from = c("cor", "model", "conv", "mlit"), values_from = total) -> total_11


#-------------------------------------------------------------------------------


# figures------------------------------------------------------------------------

# setting the theme
my_theme <- theme_grey() + theme(panel.background = element_blank(),
                                 legend.position = 'bottom',
                                 axis.title = element_text(face = "bold"),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(), 
                                 strip.text = element_text(face = 'bold', size = 6.0),
                                 strip.background = element_rect(colour = 'grey79'), 
                                 axis.line = element_line(colour = 'grey79', linewidth = .35),
                                 axis.text = element_text(size = 6.0)) 
theme_set(my_theme)

# figures-----------------------------------------------------------------------
results_full %>%
  drop_na() %>%
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
       filename = 'Fit-logl-complete-space-study-2-listwise.png', width = 10, height = 8,
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
                               'Model 2C' = "#b2182b")) +
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
       filename = 'Fit-logl-positive-subspace-study-2-listwise.png', width = 10, height = 8,
       units = 'in')
#-------------------------------------------

