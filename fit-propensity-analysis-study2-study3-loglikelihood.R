#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 09/09/2025
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Data pre-processing for studies 2 and 3
#-------------------------------------------------------------------------------

## Load packages ---------------------------------------------------------------
pkgs <- c('tidyverse', 'ggdist', 'ggthemes', 'effectsize')

# check if packages are available
rlang::check_installed(pkgs)

# loading packages
lapply(pkgs, library, character.only = TRUE)

# Load full data ---------------------------------------------------------------
load(here::here("Online-Supplement","Data-Files", "output_poscor_2_all_conditions.RData"))

# positive correlations
study2_pos_fits <- output_conditions

# deleting the data
rm(output_conditions)

load(here::here("Online-Supplement","Data-Files", "output_all_2_conditions.RData"))

## All correlations
study2_all_fits <- output_conditions

# deleting the data
rm(output_conditions)

## Combine all output in one data frame -----------------------------------------
for(i in 1:48) {
  tmp <- study2_all_fits[[i]]
  tmp$condition <- i
  
  if(i == 1) {
    study2_all_total <- tmp
  } else {
    study2_all_total <- bind_rows(study2_all_total, tmp)
  }
}

for(i in 1:48) {
  tmp <- study2_pos_fits[[i]]
  tmp$condition <- i
  
  if(i == 1) {
    study2_pos_total <- tmp
  } else {
    study2_pos_total <- bind_rows(study2_pos_total, tmp)
  }
}

rm(tmp, study2_all_fits, study2_pos_fits)


# Study 2 conditions -----------------------------------------------------------
study2_conditions <- expand.grid(bounds = c("none", "pos.var"),
                                 ident = c("refmarker", "unitvar"),
                                 mlit = c('150', '1500'),
                                 conv = c("F", "T"),
                                 model = 1:3)

## merging the data frame with the conditions
study2_conditions$condition <- 1:48

## Overall average and FP across 1 million matrices ----------------------------

## FOR CUTOFFS STUDY 3 ---------------------------------------------------------
study2_pos_total %>% summarize(m = mean(logl, na.rm = T))
study2_all_total  %>% summarize(m = mean(logl, na.rm = T))

# the overall means in both data spaces given that we used list wise deletion to deal with none convergence
# Positive: -27013.04 and All: -32754.76

# saving the results
save(study2_pos_total, file = here::here("Online-Supplement","Data-Files", "study2_pos_total.RData"))
save(study2_all_total, file = here::here("Online-Supplement","Data-Files", "study2_all_total.RData"))

#-------------------------------------------------------------------------------
# loading the datasets
load(here::here("Online-Supplement","Data-Files", "study2_pos_total.RData"))
load(here::here("Online-Supplement","Data-Files", "study2_all_total.RData"))
#-------------------------------------------------------------------------------

## Average fit value
study2_all_avg <- 
  study2_all_total %>% group_by(model, conv, mlit, ident, bounds) %>%
  summarize(m_cfi = mean(cfi, na.rm = T),
            m_srmr = mean(srmr, na.rm = T),
            m_rmsea = mean(rmsea, na.rm = T),
            m_logl = mean(logl, na.rm = T),
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()

# saving the results
save(study2_all_avg, file = here::here("Online-Supplement","Data-Files", "study2_all_averages.RData"))

study2_pos_avg <- 
  study2_pos_total %>% group_by(model, conv, mlit, ident, bounds) %>%
  summarize(m_cfi = mean(cfi, na.rm = T),
            m_srmr = mean(srmr, na.rm = T),
            m_rmsea = mean(rmsea, na.rm = T),
            m_logl = mean(logl, na.rm = T),
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()

# saving the results
save(study2_pos_avg, file = here::here("Online-Supplement","Data-Files", "study2_pos_averages.RData"))
#-------------------------------------------------------------------------------

# Fitting propensity analysis in study 2----------------------------------------

#-------------------------------------------------------------------------------
# here we have replaced NAs with assumed badness of fit
study2_all_logl_fp <- study2_all_total %>% 
  #drop_na() %>%
  mutate(logl = if_else(is.na(logl), -40000, logl)) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -32755) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()

# saving the results
save(study2_all_logl_fp, file = here::here("Online-Supplement","Data-Files", "study2_all_logl_fp_badness.RData"))

study2_pos_logl_fp <- study2_pos_total %>% 
  #drop_na() %>%
  mutate(logl = if_else(is.na(logl), -40000, logl)) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -27103) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()

# saving the results
save(study2_pos_logl_fp, file = here::here("Online-Supplement","Data-Files","study2_pos_logl_fp_badness.RData"))
#-------------------------------------------------------------------------------
# without replacing the NAs
#-------------------------------------------------------------------------------
study2_all_logl_fp_1 <- study2_all_total %>% 
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = mean(logl > -32755, na.rm = TRUE), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()

# saving the results
save(study2_all_logl_fp_1, file = here::here("Online-Supplement","Data-Files", "study2_all_logl_fp_true.RData"))

study2_pos_logl_fp_1 <- study2_pos_total %>% 
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = mean(logl > -27103, na.rm = TRUE), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()

# saving the results
save(study2_pos_logl_fp_1, file = here::here("Online-Supplement","Data-Files", "study2_pos_logl_fp_true.RData"))
#-------------------------------------------------------------------------------

# removing the files in the environment
rm(list = ls())

## Same for Study 3 ------------------------------------------------------------

# loading the data set
load(here::here("Online-Supplement","Data-Files","output_all_3_conditions.RData"))

study3_all_fits <- output_conditions

# loading the data set
load(here::here("Online-Supplement","Data-Files", "output_poscor_3_all_conditions.RData"))

study3_pos_fits <- output_conditions

rm(output_conditions)

## Combine all output in one dataframe -----------------------------------------
for(i in 1:48) {
  tmp <- study3_all_fits[[i]]
  tmp$condition <- i
  
  if(i == 1) {
    study3_all_total <- tmp
  } else {
    study3_all_total <- bind_rows(study3_all_total, tmp)
  }
}

for(i in 1:48) {
  tmp <- study3_pos_fits[[i]]
  tmp$condition <- i
  
  if(i == 1) {
    study3_pos_total <- tmp
  } else {
    study3_pos_total <- bind_rows(study3_pos_total, tmp)
  }
}

rm(tmp, study3_all_fits, study3_pos_fits)

# Study 3 conditions -----------------------------------------------------------
study3_conditions <- expand.grid(model = 1:3,
                                 bounds = c("none", "pos.var"),
                                 ident = c("refmarker", "unitvar"),
                                 mlit = c('150', '1500'),
                                 conv = c("F", "T"))

study3_conditions$condition <- 1:48



## FOR CUTOFFS STUDY 3 ---------------------------------------------------------
study3_pos_total %>% group_by(model) %>% summarize(m = mean(logl, na.rm = T))

# model       m
# <fct>   <dbl>
#   1 1     -10306.
# 2 2     -30479.
# 3 3     -50076.

study3_all_total %>% group_by(model) %>% summarize(m = mean(logl, na.rm = T))
# 
# model       m
# <fct>   <dbl>
#   1 1     -11460.
# 2 2     -36999.
# 3 3     -62493.

# saving the results------------------------------------------------------------
save(study3_pos_total, file = here::here("Online-Supplement","Data-Files", "study3_pos_total.RData"))

save(study3_all_total, file = here::here("Online-Supplement","Data-Files", "study3_all_total.RData"))

#-------------------------------------------------------------------------------
# loading the datasets
load(file = here::here("Online-Supplement","Data-Files", "study3_pos_total.RData"))

load(file = here::here("Online-Supplement","Data-Files", "study3_all_total.RData"))
#-------------------------------------------------------------------------------
## Average fit value
study3_all_avg <- 
  study3_all_total %>% group_by(model, conv, mlit, ident, bounds) %>%
  summarize(m_cfi = mean(cfi, na.rm = T),
            m_srmr = mean(srmr, na.rm = T),
            m_rmsea = mean(rmsea, na.rm = T),
            m_logl = mean(logl, na.rm = T),
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()

# saving the results
save(study3_all_avg, file = here::here("Online-Supplement","Data-Files","study3_all_averages.RData"))

study3_pos_avg <- 
  study3_pos_total %>% group_by(model, conv, mlit, ident, bounds) %>%
  summarize(m_cfi = mean(cfi, na.rm = T),
            m_srmr = mean(srmr, na.rm = T),
            m_rmsea = mean(rmsea, na.rm = T),
            m_logl = mean(logl, na.rm = T),
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()

# saving the results
save(study3_pos_avg, file = here::here("Online-Supplement","Data-Files", "study3_pos_averages.RData"))
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# here we have replaced NAs with assumed badness of fit
study3_all_logl_fp_model_1 <- study3_all_total %>%
  filter(model == 1) %>%
  mutate(logl = if_else(is.na(logl), -20000, logl)) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -11460) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()
  
study3_all_logl_fp_model_2 <- study3_all_total %>% 
    filter(model == 2) %>%
  mutate(logl = if_else(is.na(logl), -40000, logl)) %>%
    group_by(model, conv, mlit, ident, bounds) %>%
    summarize(fp_logl = (sum(logl > -36999) / n()), # CUTOFF
              condition_num = min(condition)) %>%
    mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
    ungroup()  

study3_all_logl_fp_model_3 <- study3_all_total %>% 
  filter(model == 3) %>%
  mutate(logl = if_else(is.na(logl), -70000, logl)) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -62493) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()  

# saving the results
save(study3_all_logl_fp_model_1, file = here::here("Online-Supplement","Data-Files","study3_all_logl_fp_mode1_1_badness.RData"))

save(study3_all_logl_fp_model_2, file = here::here("Online-Supplement","Data-Files", "study3_all_logl_fp_mode1_2_badness.RData"))

save(study3_all_logl_fp_model_3, file = here::here("Online-Supplement","Data-Files","study3_all_logl_fp_mode1_3_badness.RData"))
#-----------------------------------------------------------------------------------------

study3_pos_logl_fp_model_1 <- study3_pos_total %>% 
  filter(model == 1) %>%
  mutate(logl = if_else(is.na(logl), -15000, logl)) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -10306) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()

study3_pos_logl_fp_model_2 <- study3_pos_total %>% 
  filter(model == 2) %>%
  mutate(logl = if_else(is.na(logl), -40000, logl)) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -30479) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()  

study3_pos_logl_fp_model_3 <- study3_pos_total %>% 
  filter(model == 3) %>%
  mutate(logl = if_else(is.na(logl), -60000, logl)) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -50076) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()  

# saving results
save(study3_pos_logl_fp_model_1, file = here::here("Online-Supplement","Data-Files", "study3_pos_logl_fp_mode1_1_badness.RData"))

save(study3_pos_logl_fp_model_2, file = here::here("Online-Supplement","Data-Files", "study3_pos_logl_fp_mode1_2_badness.RData"))

save(study3_pos_logl_fp_model_3, file = here::here("Online-Supplement","Data-Files", "study3_pos_logl_fp_mode1_3_badness.RData"))
#-------------------------------------------------------------------------------
# default(list-wise)
#-------------------------------------------------------------------------------
study3_all_logl_fp_model_11 <- study3_all_total %>% 
  filter(model == 1) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = mean(logl > -11460, na.rm = TRUE), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()

study3_all_logl_fp_model_22 <- study3_all_total %>% 
  filter(model == 2) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = mean(logl > -36999, na.rm = TRUE), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()  

study3_all_logl_fp_model_33 <- study3_all_total %>% 
  drop_na() %>%
  filter(model == 3) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = mean(logl > -62493, na.rm = TRUE), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()  

# saving the results
save(study3_all_logl_fp_model_11, file = here::here("Online-Supplement","Data-Files", "study3_all_logl_fp_mode1_1_true.RData"))

save(study3_all_logl_fp_model_22, file = here::here("Online-Supplement","Data-Files", "study3_all_logl_fp_mode1_2_true.RData"))

save(study3_all_logl_fp_model_33, file = here::here("Online-Supplement","Data-Files", "study3_all_logl_fp_mode1_3_true.RData"))
#-----------------------------------------------------------------------------------------
# not replacing the NAs
#-------------------------------------------------------------------------------
study3_pos_logl_fp_model_11 <- study3_pos_total %>% 
  filter(model == 1) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = mean(logl > -10306, na.rm = TRUE), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()

study3_pos_logl_fp_model_22 <- study3_pos_total %>% 
  filter(model == 2) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = mean(logl > -30479, na.rm = TRUE), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()  

study3_pos_logl_fp_model_33 <- study3_pos_total %>% 
  filter(model == 3) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = mean(logl > -50076, na.rm = TRUE), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()  

# saving results
save(study3_pos_logl_fp_model_11, file = here::here("Online-Supplement","Data-Files","study3_pos_logl_fp_mode1_1_true.RData"))
save(study3_pos_logl_fp_model_22, file = here::here("Online-Supplement","Data-Files", "study3_pos_logl_fp_mode1_2_true.RData"))
save(study3_pos_logl_fp_model_33, file = here::here("Online-Supplement","Data-Files", "study3_pos_logl_fp_mode1_3_true.RData"))

#-------------------------------------------------------------------------------
# removing everything in the environment
rm(list = ls())
#-------------------------------------------------------------------------------