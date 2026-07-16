#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 09/09/2025
#-------------------------------------------------------------------------------

# loading packages
pkgs <- c('tidyverse', 'gt', 'ggdist', 'ggthemes')

lapply(pkgs, library, character.only = T)

## Load full data
load(here::here("Online-Supplement","Data-Files" ,"output_poscor_3_all_conditions.RData"))

# positive correlations
study3_pos_fits <- output_conditions

load(here::here("Online-Supplement","Data-Files", "output_all_3_conditions.RData"))

study3_all_fits <- output_conditions

## Select range of data space sizes ----
sizes <- c(50, 100, 500, 1000, 5000, 10000, 50000, 100000, 150000, 175000, 200000, 500000, 750000)

## Set up functions to compute FP (proportion of fit index values above/below cutoff) ----

# Positive correlation

# subsets
model1_subsets <- lapply(study3_pos_fits, function(df) subset(df, model == 1))

# removing the empty data frames
model1_subsets <- model1_subsets[sapply(model1_subsets,nrow) > 0]

model2_subsets <- lapply(study3_pos_fits, function(df) subset(df, model == 2))

# removing the empty data frames
model2_subsets <- model2_subsets[sapply(model2_subsets,nrow) > 0]


model3_subsets <- lapply(study3_pos_fits, function(df) subset(df, model == 3))

# removing the empty data frames
model3_subsets <- model3_subsets[sapply(model3_subsets,nrow) > 0]


#-----------------------------------------------------------------------
# Fitting propensity for each respective model (Positive correlations)
#-----------------------------------------------------------------------
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


# functions to sample X fit values and compute FP
samfits <- function(fits_all, size = 100) {
  
  rows <- sample(1:nrow(fits_all[[1]]), size, replace = FALSE)
  
  # container
  results <- matrix(NA, length(fits_all), 4)
  
  # renaming columns
  colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr')
  
  for(m in seq_along(fits_all)){
    # sub-setting  
    fits = fits_all[[m]][rows,]
    
    # CFI
    cfi<- cfi.fp_1(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp_1(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp_1(fits[,"srmr"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi,rmsea,srmr)
    
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

# functions to sample X fit values and compute FP
samfits_2 <- function(fits_all, size = 100) {
  
  rows <- sample(1:nrow(fits_all[[1]]), size, replace = FALSE)
  
  # container
  results <- matrix(NA, length(fits_all), 4)
  
  # renaming columns
  colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr')
  
  for(m in seq_along(fits_all)){
    # sub-setting  
    fits = fits_all[[m]][rows,]
    
    # CFI
    cfi<- cfi.fp_2(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp_2(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp_2(fits[,"srmr"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi,rmsea,srmr)
    
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

# functions to sample X fit values and compute FP
samfits_3 <- function(fits_all, size = 100) {
  
  rows <- sample(1:nrow(fits_all[[1]]), size, replace = FALSE)
  
  # container
  results <- matrix(NA, length(fits_all), 4)
  
  # renaming columns
  colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr')
  
  for(m in seq_along(fits_all)){
    # sub-setting  
    fits = fits_all[[m]][rows,]
    
    # CFI
    cfi<- cfi.fp_3(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp_3(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp_3(fits[,"srmr"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi,rmsea,srmr)
    
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
  tmp <- (lapply(sizes, function(size) samfits(model1_subsets, size = size)))
  tmp <- (do.call(rbind, tmp))
  tmp <- cbind(tmp, rep(sizes, each = 16))
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr', 'size')
  
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

res_pos_model1 <- merge(fp.size.results_1, study_con_model_1)

res_pos_model1$model <- factor(res_pos_model1$model)

# saving the results
saveRDS(fp.size.results_1, file = here::here("Online-Supplement","Data-Files", "study3_poscor_model_1_dataspacesize.RDS"))
saveRDS(res_pos_model1, file = here::here("Online-Supplement", "Data-Files", "results_poscor_study_3_model_1_dataspacesize.RDS"))

# Model 2

# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits_2(model2_subsets, size = size)))
  tmp <- (do.call(rbind, tmp))
  tmp <- cbind(tmp, rep(sizes, each = 16))
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr', 'size')
  
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

res_pos_model2 <- merge(fp.size.results_2, study_con_model_2)

res_pos_model2$model <- factor(res_pos_model2$model)

# saving the results
saveRDS(fp.size.results_2, file = here::here("Online-Supplement","Data-Files", "study3_poscor_model_2_dataspacesize.RDS"))
saveRDS(res_pos_model2, file = here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_2_dataspacesize.RDS"))

# Model 3

# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits_3(model3_subsets, size = size)))
  tmp <- (do.call(rbind, tmp))
  tmp <- cbind(tmp, rep(sizes, each = 16))
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr', 'size')
  
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

res_pos_model3 <- merge(fp.size.results_3, study_con_model_3)

res_pos_model3$model <- factor(res_pos_model3$model)

# saving the results
saveRDS(fp.size.results_3, file = here::here("Online-Supplement","Data-Files", "study3_poscor_model_3_dataspacesize.RDS"))
saveRDS(res_pos_model3, file = here::here("Online-Supplement","Data-Files" ,"results_poscor_study_3_model_3_dataspacesize.RDS"))

#-------------------------------------------------------------------------------
# Visualizing the data spaces
#-------------------------------------------------------------------------------
my_theme <- theme_grey() + theme(panel.background = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 6.0),
        strip.background = element_rect(colour = 'grey79'), axis.line = element_line(colour = 'grey79', linewidth = .35),
        axis.text = element_text(size = 6.0)) 

set_theme(my_theme)

# Model 1:

# loading the data
res_pos_model1 <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_1_dataspacesize.RDS"))

plot_1 <- res_pos_model1 %>%
  pivot_longer(cols = cfi:srmr, names_to = "index") %>%
  filter(size <= 10000) %>%
  mutate(model = factor(model,levels = c(1, 2, 3),
    labels = c("Model 3A", "Model 3B", "Model 3C")),
    index = factor(index, levels = c("cfi", "rmsea", "srmr"),labels = c("CFI", "RMSEA", "SRMR")),
    bounds = factor(bounds,levels = c("none", "pos.var"),labels = c("Unrestricted", "Positive")),
    ident = factor(ident, levels = c("refmarker", "unitvar"), labels = c("Reference", "Unit Variance")),
    mlit = factor(mlit, levels = c("150", "1500")),
    conv = factor(conv, levels = c("F", "T"), labels = c("Default", "Forced"))) %>%
  ggplot(aes( x = size, y = value, group = model,
      fill = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)), .width = 0.95,
    linewidth = 0.4, alpha = 0.9) +
  facet_grid( rows = vars(index, mlit), cols = vars(conv, ident),
    scales = "free") +
  labs(x = "Data Space Size", y = "Fitting Propensity (proportion < or > cutoff)",
    fill_ramp = "Confidence Level") +
  scale_fill_manual( name = "Model",
    values = c(
      "Model 3A" = "#276419",
      "Model 3B" = "#8e0152",
      "Model 3C" = "#b2182b")) + scale_x_continuous( breaks = seq(0, 1e4, 2e3),
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),sec.axis = sec_axis(~ .,
      name = "Positive Subspace", breaks = NULL, labels = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)), sec.axis = sec_axis( ~ ., name = "Fit Index",
      breaks = NULL,
      labels = NULL))

# Save it
ggsave(plot = plot_1, filename = 'poscor_3_model_1_dataspace-badness.png',path =  here::here("Plots", "Study Three"), width = 16, height = 12, units = 'in')


# Model 2:

# loading the data
res_pos_model2 <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_2_dataspacesize.RDS"))

plot_2 <-  res_pos_model2 %>%
  pivot_longer(cols = cfi:srmr, names_to = "index") %>%
  filter(size <= 10000) %>%
  mutate(model = factor(model,levels = c(1, 2, 3),
        labels = c("Model 3A", "Model 3B", "Model 3C")),
         index = factor(index, levels = c("cfi", "rmsea", "srmr"),labels = c("CFI", "RMSEA", "SRMR")),
         bounds = factor(bounds,levels = c("none", "pos.var"),labels = c("Unrestricted", "Positive")),
         ident = factor(ident, levels = c("refmarker", "unitvar"), labels = c("Reference", "Unit Variance")),
         mlit = factor(mlit, levels = c("150", "1500")), conv = factor(conv, levels = c("F", "T"),
        labels = c("Default", "Forced"))) %>%
  ggplot(aes( x = size, y = value, group = model, fill = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)), .width = 0.95,
  linewidth = 0.4, alpha = 0.9) + facet_grid( rows = vars(index, mlit), cols = vars(conv, ident), scales = "free") +
  labs(x = "Data Space Size", y = "Fitting Propensity (proportion < or > cutoff)", fill_ramp = "Confidence Level") +
  scale_fill_manual(name = "Model",
                     values = c(
                       "Model 3A" = "#276419",
                       "Model 3B" = "#8e0152",
                       "Model 3C" = "#b2182b")) + 
  scale_x_continuous( breaks = seq(0, 1e4, 2e3),
  labels = scales::label_number(scale_cut = scales::cut_short_scale()),sec.axis = sec_axis(~ .,
  name = "Positive Subspace", breaks = NULL, labels = NULL)) + scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)), sec.axis = sec_axis( ~ ., name = "Fit Index",
 breaks = NULL, labels = NULL))

ggsave(plot = plot_2, path =  here::here("Plots", "Study Three"), filename = 'poscor_3_model_2_dataspace-badness.png', width = 16, height = 12,
       units = 'in')


# Model 3:

# loading the data
res_pos_model3 <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_3_dataspacesize.RDS"))

plot_3 <- res_pos_model3 %>%
  pivot_longer(cols = cfi:srmr, names_to = "index") %>%
  filter(size <= 10000) %>%
  mutate(model = factor(model,levels = c(1, 2, 3), labels = c("Model 3A", "Model 3B", "Model 3C")),
         index = factor(index, levels = c("cfi", "rmsea", "srmr"),labels = c("CFI", "RMSEA", "SRMR")),
         bounds = factor(bounds,levels = c("none", "pos.var"),labels = c("Unrestricted", "Positive")),
         ident = factor(ident, levels = c("refmarker", "unitvar"), labels = c("Reference", "Unit Variance")),
         mlit = factor(mlit, levels = c("150", "1500")), conv = factor(conv, levels = c("F", "T"),labels = c("Default", "Forced"))) %>%
  ggplot(aes(x = size, y = value, group = model, fill = model)) + 
  stat_lineribbon(aes(fill_ramp = after_stat(level)), .width = 0.95,
  linewidth = 0.4, alpha = 0.9) +
  facet_grid( rows = vars(index, mlit), cols = vars(conv, ident),  scales = "free") +
  labs(x = "Data Space Size", y = "Fitting Propensity (proportion < or > cutoff)",
       fill_ramp = "Confidence Level") + scale_fill_manual( name = "Model",
                     values = c(
                       "Model 3A" = "#276419",
                       "Model 3B" = "#8e0152",
                       "Model 3C" = "#b2182b")) + scale_x_continuous( breaks = seq(0, 1e4, 2e3),
labels = scales::label_number(scale_cut = scales::cut_short_scale()),sec.axis = sec_axis(~ .,
name = "Positive Subspace", breaks = NULL, labels = NULL)) + scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
sec.axis = sec_axis( ~ ., name = "Fit Index", breaks = NULL, labels = NULL))

ggsave(plot = plot_3, path =  here::here("Plots", "Study Three"), filename = 'poscor_3_model_3_dataspace-badness.png', width = 16, height = 12,
       units = 'in')


# All correlation

# subsets
model1_subsets <- lapply(study3_all_fits, function(df) subset(df, model == 1))

# removing the empty data frames
model1_subsets <- model1_subsets[sapply(model1_subsets,nrow) > 0]

model2_subsets <- lapply(study3_all_fits, function(df) subset(df, model == 2))

# removing the empty data frames
model2_subsets <- model2_subsets[sapply(model2_subsets,nrow) > 0]

model3_subsets <- lapply(study3_all_fits, function(df) subset(df, model == 3))

# removing the empty data frames
model3_subsets <- model3_subsets[sapply(model3_subsets,nrow) > 0]


# CFI plots for model 1
for(i in seq_along(model1_subsets)) {
  df <- model1_subsets[[i]]
  p <- ggplot(df, aes(x = cfi)) +
    geom_histogram(bins = 30, alpha = 0.6, color = "black", na.rm = T) +
    labs(title =  "CFI Distribution for Model 1",
         x = "CFI",
         y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # saving the plots
  ggsave(plot = gridExtra::grid.arrange(p, ncol = 1), path = here::here("Plots", "Study Three"), filename = paste0('CFI distribution_', i, '_model11.png'), width = 16, height = 12,
         units = 'in')
}




# CFI plots for model 2
for(i in seq_along(model2_subsets)) {
  df <- model2_subsets[[i]]
  p <- ggplot(df, aes(x = cfi)) +
    geom_histogram(bins = 30, alpha = 0.6, color = "black", na.rm = T) +
    labs(title =  "CFI Distribution for Model 2",
         x = "CFI",
         y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # saving the plots
  ggsave(plot = gridExtra::grid.arrange(p, ncol = 1), path =  here::here("Plots", "Study Three"), filename = paste0('CFI distribution_', i, '_model21.png'), width = 16, height = 12,
         units = 'in')
}


# CFI plots for model 3
for(i in seq_along(model3_subsets)) {
  df <- model3_subsets[[i]]
  p <- ggplot(df, aes(x = cfi)) +
    geom_histogram(bins = 30, alpha = 0.6, color = "black", na.rm = T) +
    labs(title =  "CFI Distribution for Model 3",
         x = "CFI",
         y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # saving the plots
  ggsave(plot = gridExtra::grid.arrange(p, ncol = 1), path =  here::here("Plots", "Study Three"), filename = paste0('CFI distribution_', i, '_model31.png'), width = 16, height = 12,
         units = 'in')
}


# SRMR plots for model 1
for(i in seq_along(model1_subsets)) {
  df <- model1_subsets[[i]]
  p <- ggplot(df, aes(x = srmr)) +
    geom_histogram(bins = 30, alpha = 0.6, color = "black", na.rm = T) +
    labs(title =  "SRMR Distribution for Model 1",
         x = "SRMR",
         y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # saving the plots
  ggsave(plot = gridExtra::grid.arrange(p, ncol = 1), path = here::here("Plots", "Study Three"), filename = paste0('SRMR distribution_', i, '_model11.png'), width = 16, height = 12,
         units = 'in')
}



# SRMR plots for model 2
for(i in seq_along(model2_subsets)) {
  df <- model2_subsets[[i]]
  p <- ggplot(df, aes(x = srmr)) +
    geom_histogram(bins = 30, alpha = 0.6, color = "black", na.rm = T) +
    labs(title =  "SRMR Distribution for Model 2",
         x = "SRMR",
         y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # saving the plots
  ggsave(plot = gridExtra::grid.arrange(p, ncol = 1), path = here::here("Plots", "Study Three"), filename = paste0('SRMR distribution_', i, '_model21.png'), width = 16, height = 12,
         units = 'in')
}


# SRMR plots for model 3
for(i in seq_along(model3_subsets)) {
  df <- model3_subsets[[i]]
  p <- ggplot(df, aes(x = srmr)) +
    geom_histogram(bins = 30, alpha = 0.6, color = "black", na.rm = T) +
    labs(title =  "SRMR Distribution for Model 3",
         x = "SRMR",
         y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # saving the plots
  ggsave(plot = gridExtra::grid.arrange(p, ncol = 1), path =  here::here("Plots", "Study Three"), filename = paste0('SRMR distribution_', i, '_model31.png'), width = 16, height = 12,
         units = 'in')
}



# RMSEA plots for model 1
for(i in seq_along(model1_subsets)) {
  df <- model1_subsets[[i]]
  p <- ggplot(df, aes(x = rmsea)) +
    geom_histogram(bins = 30, alpha = 0.6, color = "black", na.rm = T) +
    labs(title =  "RMSEA Distribution for Model 1",
         x = "RMSEA",
         y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # saving the plots
  ggsave(plot = gridExtra::grid.arrange(p, ncol = 1), path =  here::here("Plots", "Study Three"), filename = paste0('RMSEA distribution_', i, '_model11.png'), width = 16, height = 12,
         units = 'in')
}



# RMSEA plots for model 2
for(i in seq_along(model2_subsets)) {
  df <- model2_subsets[[i]]
  p <- ggplot(df, aes(x = rmsea)) +
    geom_histogram(bins = 30, alpha = 0.6, color = "black", na.rm = T) +
    labs(title =  "RMSEA Distribution for Model 2",
         x = "RMSEA",
         y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # saving the plots
  ggsave(plot = gridExtra::grid.arrange(p, ncol = 1), path =  here::here("Plots", "Study Three"), filename = paste0('RMSEA distribution_', i, '_model21.png'), width = 16, height = 12,
         units = 'in')
}


## RMSEA plots for model 3
for(i in seq_along(model3_subsets)) {
  df <- model3_subsets[[i]]
  p <- ggplot(df, aes(x = rmsea)) +
    geom_histogram(bins = 30, alpha = 0.6, color = "black", na.rm = T) +
    labs(title =  "RMSEA Distribution for Model 3",
         x = "RMSEA",
         y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # saving the plots
  ggsave(plot = gridExtra::grid.arrange(p, ncol = 1), path =  here::here("Plots", "Study Three"), filename = paste0('RMSEA distribution_', i, '_model31.png'), width = 16, height = 12,
         units = 'in')
}


## Based on the different histograms, it was evident that we assign each model different cut-offs for this simulation case.
# Hence, to assess the fitting propensities we will use different thresholds for each model.


# Model 1: CFI >= [.55, .65], RMSEA <= [.40, .55], SRMR <= [.15. .17]

# Model 2: CFI >= [.40, .45], RMSEA <= [.23, .25], SRMR <= [.08, .09]

# Model 3: CFI >= [.42, .43], RMSEA <= [.17, .20], SRMR <= [.071, .073]

#-------------------------------------------------------------------------------
# Fitting propensity for each respective model (Positive correlations)

# Model 1
# function to get prop CFI >= .60
cfi.fp_1 <- function(fits, cutoff = .60) {
  fits[is.na(fits)] <- 0
  return((sum(fits >= cutoff) / length(fits)))
}


# function to get prop RMSEA < .45
rmsea.fp_1 <- function(fits, cutoff = .45) {
  fits[is.na(fits)] <- 1
  return((sum(fits < cutoff) / length(fits)))
}


# function to get prop SRMR < .15
srmr.fp_1 <- function(fits, cutoff = .15) {
  fits[is.na(fits)] <- 1
  return((sum(fits < cutoff) / length(fits)))
}


# functions to sample X fit values and compute FP
samfits <- function(fits_all, size = 100) {
  
  rows <- sample(1:nrow(fits_all[[1]]), size, replace = FALSE)
  
  # container
  results <- matrix(NA, length(fits_all), 4)
  
  # renaming columns
  colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr')
  
  for(m in seq_along(fits_all)){
    # sub-setting  
    fits = fits_all[[m]][rows,]
    
    # CFI
    cfi<- cfi.fp_1(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp_1(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp_1(fits[,"srmr"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi,rmsea,srmr)
    
    # saving results
    results[m,] <- res
    
  }
  
  # output
  return(results)
  
  
}
#-------------------------------------------------------------------------------
# Model 2
# function to get prop CFI >= .45
cfi.fp_2 <- function(fits, cutoff = .45) {
  fits[is.na(fits)] <- 0
  return((sum(fits >= cutoff) / length(fits)))
}


# function to get prop RMSEA < .23
rmsea.fp_2 <- function(fits, cutoff = .23) {
  fits[is.na(fits)] <- 1
  return((sum(fits < cutoff) / length(fits)))
}


# function to get prop SRMR < .15
srmr.fp_2 <- function(fits, cutoff = .15) {
  fits[is.na(fits)] <- 1
  return((sum(fits < cutoff) / length(fits)))
}

# functions to sample X fit values and compute FP
samfits_2 <- function(fits_all, size = 100) {
  
  rows <- sample(1:nrow(fits_all[[1]]), size, replace = FALSE)
  
  # container
  results <- matrix(NA, length(fits_all), 4)
  
  # renaming columns
  colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr')
  
  for(m in seq_along(fits_all)){
    # sub-setting  
    fits = fits_all[[m]][rows,]
    
    # CFI
    cfi<- cfi.fp_2(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp_2(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp_2(fits[,"srmr"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi,rmsea,srmr)
    
    # saving results
    results[m,] <- res
    
  }
  
  # output
  return(results)
}  
#-------------------------------------------------------------------------------

# Model 3
# function to get prop CFI >= .42
cfi.fp_3 <- function(fits, cutoff = .42) {
  fits[is.na(fits)] <- 0
  return((sum(fits >= cutoff) / length(fits)))
}


# function to get prop RMSEA < .17
rmsea.fp_3 <- function(fits, cutoff = .17) {
  fits[is.na(fits)] <- 1
  return((sum(fits < cutoff) / length(fits)))
}


# function to get prop SRMR < .071
srmr.fp_3 <- function(fits, cutoff = .073) {
  fits[is.na(fits)] <- 1
  return((sum(fits < cutoff) / length(fits)))
}

# functions to sample X fit values and compute FP
samfits_3 <- function(fits_all, size = 100) {
  
  rows <- sample(1:nrow(fits_all[[1]]), size, replace = FALSE)
  
  # container
  results <- matrix(NA, length(fits_all), 4)
  
  # renaming columns
  colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr')
  
  for(m in seq_along(fits_all)){
    # sub-setting  
    fits = fits_all[[m]][rows,]
    
    # CFI
    cfi<- cfi.fp_3(fits[,"cfi"])
    
    # RMSEA
    rmsea<- rmsea.fp_3(fits[,"rmsea"])
    
    # SRMR
    srmr <- srmr.fp_3(fits[,"srmr"])
    
    # data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
    res <- c(m, cfi,rmsea,srmr)
    
    # saving results
    results[m,] <- res
    
  }
  
  # output
  return(results)
}  

# seed
set.seed(3242)

#-------------------------------------------------------------------------------
## Get FP for different space sizes ----

# Model 1
# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits(model1_subsets, size = size)))
  tmp <- (do.call(rbind, tmp))
  tmp <- cbind(tmp, rep(sizes, each = 16))
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr', 'size')
  
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

res_all_model1 <- merge(fp.size.results_1, study_con_model_1)

res_all_model1$model <- factor(res_all_model1$model)

# saving the results
saveRDS(fp.size.results_1, file = here::here("Online-Supplement", "Data-Files", "study3_all_model_1_dataspacesize.RDS"))

saveRDS(res_all_model1, file = here::here("Online-Supplement","Data-Files", "results_allcor_study_3_model_1_dataspacesize.RDS"))

# Model 2

# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits_2(model2_subsets, size = size)))
  tmp <- (do.call(rbind, tmp))
  tmp <- cbind(tmp, rep(sizes, each = 16))
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr', 'size')
  
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

res_all_model2 <- merge(fp.size.results_2, study_con_model_2)

res_all_model2$model <- factor(res_all_model2$model)

# saving the results
saveRDS(fp.size.results_2, file = here::here("Online-Supplement","Data-Files", "study3_all_model_2_dataspacesize.RDS"))

saveRDS(res_all_model2, file = here::here("Online-Supplement","Data-Files", "results_allcor_study_3_model_2_dataspacesize.RDS"))


# Model 3

# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  tmp <- (lapply(sizes, function(size) samfits_3(model3_subsets, size = size)))
  tmp <- (do.call(rbind, tmp))
  tmp <- cbind(tmp, rep(sizes, each = 16))
  colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr', 'size')
  
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

res_all_model3 <- merge(fp.size.results_3, study_con_model_3)

res_all_model3$model <- factor(res_all_model3$model)

# saving the results
saveRDS(fp.size.results_3, file = here::here("Online-Supplement","Data-Files", "study3_all_model_3_dataspacesize.RDS"))

saveRDS(res_all_model3, file = here::here("Online-Supplement","Data-Files", "results_allcor_study_3_model_3_dataspacesize.RDS"))


#-------------------------------------------------------------------------------
# Visualizing the data spaces- use tables
#-------------------------------------------------------------------------------

# Model 1:

# loading the data
res_all_model1 <- readRDS(here::here("Online-Supplement","Data-Files", "results_allcor_study_3_model_1_dataspacesize.RDS"))

plot_1 <-res_all_model1 %>%
  pivot_longer(cols = cfi:srmr, names_to = "index") %>%
  filter(size <= 10000) %>%
  mutate(model = factor(model,levels = c(1, 2, 3), labels = c("Model 3A", "Model 3B", "Model 3C")),
         index = factor(index, levels = c("cfi", "rmsea", "srmr"),labels = c("CFI", "RMSEA", "SRMR")),
         bounds = factor(bounds,levels = c("none", "pos.var"),labels = c("Unrestricted", "Positive")),
         ident = factor(ident, levels = c("refmarker", "unitvar"), labels = c("Reference", "Unit Variance")),
         mlit = factor(mlit, levels = c("150", "1500")), conv = factor(conv, levels = c("F", "T"), labels = c("Default", "Forced"))) %>%
  ggplot(aes( x = size, y = value, group = model, fill = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)), .width = 0.95, linewidth = 0.4, alpha = 0.9) +
  facet_grid( rows = vars(index, mlit), cols = vars(conv, ident), scales = "free") +
  labs(x = "Data Space Size", y = "Fitting Propensity (proportion < or > cutoff)",
       fill_ramp = "Confidence Level") +
  scale_fill_manual( name = "Model",
                     values = c(
                       "Model 3A" = "#276419",
                       "Model 3B" = "#8e0152",
                       "Model 3C" = "#b2182b")) + scale_x_continuous( breaks = seq(0, 1e4, 2e3),
  labels = scales::label_number(scale_cut = scales::cut_short_scale()),sec.axis = sec_axis(~ .,
  name = "Complete Space", breaks = NULL, labels = NULL)) + scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)), 
  sec.axis = sec_axis( ~ ., name = "Fit Index", breaks = NULL, labels = NULL))


ggsave(plot = plot_1, path =  here::here("Plots", "Study Three"), filename = 'allcor_3_model_1_dataspace-badness.png', width = 16, height = 12,
       units = 'in')

# Model 2:

# loading the data
res_all_model2 <- readRDS(here::here("Online-Supplement","Data-Files", "results_allcor_study_3_model_2_dataspacesize.RDS"))

plot_2 <-res_all_model2 %>%
  pivot_longer(cols = cfi:srmr, names_to = "index") %>%
  filter(size <= 10000) %>%
  mutate(model = factor(model,levels = c(1, 2, 3),labels = c("Model 3A", "Model 3B", "Model 3C")),
         index = factor(index, levels = c("cfi", "rmsea", "srmr"),labels = c("CFI", "RMSEA", "SRMR")),
         bounds = factor(bounds,levels = c("none", "pos.var"),labels = c("Unrestricted", "Positive")),
         ident = factor(ident, levels = c("refmarker", "unitvar"), labels = c("Reference", "Unit Variance")),
         mlit = factor(mlit, levels = c("150", "1500")), 
         conv = factor(conv, levels = c("F", "T"), labels = c("Default", "Forced"))) %>%
  ggplot(aes( x = size, y = value, group = model, fill = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)), .width = 0.95, linewidth = 0.4, alpha = 0.9) +
  facet_grid( rows = vars(index, mlit), cols = vars(conv, ident), scales = "free") +
  labs(x = "Data Space Size", y = "Fitting Propensity (proportion < or > cutoff)",
  fill_ramp = "Confidence Level") +
  scale_fill_manual( name = "Model",
                     values = c(
                       "Model 3A" = "#276419",
                       "Model 3B" = "#8e0152",
                       "Model 3C" = "#b2182b")) + scale_x_continuous( breaks = seq(0, 1e4, 2e3),
  labels = scales::label_number(scale_cut = scales::cut_short_scale()),sec.axis = sec_axis(~ .,
  name = "Complete Space", breaks = NULL, labels = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)), 
  sec.axis = sec_axis( ~ ., name = "Fit Index", breaks = NULL, labels = NULL))

ggsave(plot = plot_2, path = here::here("Plots", "Study Three"), filename = 'allcor_3_model_2_dataspace-badness.png', width = 16, height = 12,
       units = 'in')


# Model 3

# loading the data
res_all_model3 <- readRDS(here::here("Online-Supplement","Data-Files","results_allcor_study_3_model_3_dataspacesize.RDS"))

plot_3 <-res_all_model3 %>%
  pivot_longer(cols = cfi:srmr, names_to = "index") %>%
  filter(size <= 10000) %>%
  mutate(model = factor(model,levels = c(1, 2, 3),labels = c("Model 3A", "Model 3B", "Model 3C")),
         index = factor(index, levels = c("cfi", "rmsea", "srmr"),labels = c("CFI", "RMSEA", "SRMR")),
         bounds = factor(bounds,levels = c("none", "pos.var"),labels = c("Unrestricted", "Positive")),
         ident = factor(ident, levels = c("refmarker", "unitvar"), labels = c("Reference", "Unit Variance")),
         mlit = factor(mlit, levels = c("150", "1500")), 
         conv = factor(conv, levels = c("F", "T"), labels = c("Default", "Forced"))) %>%
  ggplot(aes( x = size, y = value, group = model, fill = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)), .width = 0.95, linewidth = 0.4, alpha = 0.9) +
  facet_grid( rows = vars(index, mlit), cols = vars(conv, ident), scales = "free") +
  labs(x = "Data Space Size", y = "Fitting Propensity (proportion < or > cutoff)",
       fill_ramp = "Confidence Level") +
  scale_fill_manual( name = "Model",
                     values = c(
                       "Model 3A" = "#276419",
                       "Model 3B" = "#8e0152",
                       "Model 3C" = "#b2182b")) + scale_x_continuous( breaks = seq(0, 1e4, 2e3),
  labels = scales::label_number(scale_cut = scales::cut_short_scale()),sec.axis = sec_axis(~ .,
  name = "Complete Space", breaks = NULL, labels = NULL)) +scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)), 
  sec.axis = sec_axis( ~ ., name = "Fit Index", breaks = NULL, labels = NULL))

ggsave(plot = plot_3, path = here::here("Plots", "Study Three"), filename = 'allcor_3_model_3_dataspace-badness.png', width = 16, height = 12,
       units = 'in')