#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 09/30/2025
#-------------------------------------------------------------------------------

# loading packages
pkgs <- c('tidyverse', 'gt', 'ggdist', 'ggthemes')

lapply(pkgs, library, character.only = T)

# Load full data ----
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

## Combine all output in one dataframe ----
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


# Study 2 conditions ----
study2_conditions <- expand.grid(bounds = c("none", "pos.var"),
                                 ident = c("refmarker", "unitvar"),
                                 mlit = c('150', '1500'),
                                 conv = c("F", "T"),
                                 model = 1:3)

## merging the data frame with the conditions
study2_conditions$condition <- 1:48

## Overall average and FP across 1 million matrices ----

## FOR CUTOFFS STUDY 3 ----
study2_pos_total %>% summarize(m = mean(logl, na.rm = T))
study2_all_total  %>% summarize(m = mean(logl, na.rm = T))

# the overall means in both data spaces given that we used listwise deletion to deal with none convergence
# Positive: -27013.04 and All: -32754.76

# saving the results
save(study2_pos_total, file = here::here("Online-Supplement","Data-Files", "study2_pos_total.RData"))

save(study2_all_total, file = here::here("Online-Supplement","Data-Files","study2_all_total.RData"))

#-------------------------------------------------------------------------------
# # loading the datasets
# load(file = "results/codes-data-files/study3_pos_total.RData")
# 
# load(file = "results/codes-data-files/study3_all_total.RData")
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
save(study2_all_avg, file = "Data-Files/study2_all_averages.RData")

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
# ## FP
# 
# ## Method 1: Non-converged is replaced with "bad fit"
# study2_all_fp1 <- 
#   study2_all_total %>% 
#   mutate(cfi = if_else(is.na(cfi), 0, cfi),
#          rmsea = if_else(is.na(rmsea), 1, rmsea),
#          srmr = if_else(is.na(srmr), 1, srmr)) %>%
#   group_by(model, conv, mlit, ident, bounds) %>%
#   summarize(fp_cfi = (sum(cfi > .16) / n()),
#             fp_srmr = (sum(srmr < .17) / n()),
#             fp_rmsea = (sum(rmsea < .27) / n()),
#             condition_num = min(condition)) %>%
#   mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
#   ungroup()
# 
# # Method 2: FP based only on converged within condition (not across models)
# study2_all_fp2 <- 
#   study2_all_total %>% 
#   drop_na() %>%
#   group_by(model, conv, mlit, ident, bounds) %>%
#   summarize(fp_cfi = (sum(cfi > .16) / n()),
#             fp_srmr = (sum(srmr < .17) / n()),
#             fp_rmsea = (sum(rmsea < .27) / n()),
#             condition_num = min(condition)) %>%
#   mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
#   ungroup()
# 
# study2_pos_fp <- 
#   study2_pos_total %>% 
#   mutate(cfi = if_else(is.na(cfi), 0, cfi),
#          rmsea = if_else(is.na(rmsea), 1, rmsea),
#          srmr = if_else(is.na(srmr), 1, srmr)) %>%
#   group_by(model, conv, mlit, ident, bounds) %>%
#   summarize(fp_cfi = (sum(cfi > .45) / n()),
#             fp_srmr = (sum(srmr < .09) / n()),
#             fp_rmsea = (sum(rmsea < .27) / n()),
#             condition_num = min(condition)) %>%
#   mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
#   ungroup()
# 
# ## Figures for average fit values
# levels_cfi <- study2_all_avg %>% arrange(m_cfi) %>%  select(condition) %>% pull()
# levels_srmr <- study2_all_avg %>% arrange(m_srmr) %>%  select(condition) %>% pull()
# levels_rmsea <- study2_all_avg %>% arrange(m_rmsea) %>%  select(condition) %>% pull()
# levels_logl <- study2_all_avg %>% arrange(m_logl) %>%  select(condition) %>% pull()
# 
# study2_all_avg %>% 
#   mutate(condition = factor(condition, levels = levels_cfi)) %>%
#   ggplot(aes(y = condition, x = m_cfi, color = as.factor(model))) +
#   geom_point() + 
#   facet_wrap(vars(model), nrow = 3, ncol = 1, scales = "free_y") +
#   theme(axis.text.x = element_text(angle = 90))
# 
# study2_all_avg %>% 
#   mutate(condition = factor(condition, levels = levels_srmr)) %>%
#   ggplot(aes(y = condition, x = m_srmr, color = as.factor(model))) +
#   geom_point() + geom_line() +
#   #facet_wrap(vars(model), nrow = 3, ncol = 1, scales = "free_y") +
#   theme(axis.text.x = element_text(angle = 90))
# 
# study2_all_avg %>% 
#   mutate(condition = factor(condition, levels = levels_rmsea)) %>%
#   ggplot(aes(y = condition, x = m_rmsea, color = as.factor(model))) +
#   geom_point() + geom_line() +
#   facet_wrap(vars(model), nrow = 3, ncol = 1, scales = "free_y") +
#   theme(axis.text.x = element_text(angle = 90))
# 
# study2_all_avg %>% 
#   mutate(condition = factor(condition, levels = levels_logl)) %>%
#   ggplot(aes(y = condition, x = m_logl, color = as.factor(model))) +
#   geom_point() +
#   #facet_wrap(vars(model), nrow = 3, ncol = 1, scales = "free_y") +
#   theme(axis.text.x = element_text(angle = 90))
# 
# study2_pos_avg %>% 
#   mutate(condition = factor(condition, levels = levels_logl)) %>%
#   ggplot(aes(y = condition, x = m_logl, color = as.factor(model))) +
#   geom_point() +
#   #facet_wrap(vars(model), nrow = 3, ncol = 1, scales = "free_y") +
#   theme(axis.text.x = element_text(angle = 90))
## For the paper: 
## Guidance in paper: Logl is best
## Include all fit indices in study 1, then conclude that logl is only consistent one so continue with just that one
## logl also makes intuitive sense because it is directly related to normalized maximum likelihood (its sum across all
## data is the denominator)
## Then in discussion you can suggest future directions that look to better understand why certain fit indices respond to model
## specification issues in a different or inconsistent way
## Restructure paper so its intro > study 1 > study 1 methods > study 1 results > study 1 mini discussion > study 2 methods
## Final discussion includes best practices 
# 
# ## Figures for FP - Positive Correlations
# levels_cfi <- study2_pos_fp %>% arrange(fp_cfi) %>%  select(condition) %>% pull()
# levels_srmr <- study2_pos_fp %>% arrange(fp_srmr) %>%  select(condition) %>% pull()
# levels_rmsea <- study2_pos_fp %>% arrange(fp_rmsea) %>%  select(condition) %>% pull()
# 
# study2_pos_fp %>% 
#   mutate(condition = factor(condition, levels = levels_cfi)) %>%
#   ggplot(aes(y = condition, x = fp_cfi, color = as.factor(model))) +
#   geom_point() + geom_line() +
#   theme(axis.text.x = element_text(angle = 90))
# 
# study2_pos_fp %>% 
#   mutate(condition = factor(condition, levels = levels_srmr)) %>%
#   ggplot(aes(y = condition, x = fp_srmr, color = as.factor(model))) +
#   geom_point() + geom_line() +
#   theme(axis.text.x = element_text(angle = 90))
# 
# study2_pos_fp %>% 
#   mutate(condition = factor(condition, levels = levels_rmsea)) %>%
#   ggplot(aes(y = condition, x = fp_rmsea, color = as.factor(model))) +
#   geom_point() + geom_line() +
#   theme(axis.text.x = element_text(angle = 90))
# 
# ## Figures for FP - All Correlations - Method 1
# levels_cfi <- study2_all_fp1 %>% group_by(conv, mlit, ident, bounds) %>% arrange(fp_cfi, .by_group = T) %>% ungroup() %>%  select(condition) %>% pull()
# 
# study2_all_fp1 %>% 
#   mutate(condition = factor(condition, levels = levels_cfi)) %>%
#   ggplot(aes(y = condition, x = fp_cfi, color = as.factor(model))) +
#   geom_point() + geom_line() +
#   theme(axis.text.x = element_text(angle = 90))
# 
# levels_srmr <- study2_all_fp1 %>% arrange(fp_srmr) %>%  select(condition) %>% pull()
# 
# study2_all_fp1 %>% 
#   mutate(condition = factor(condition, levels = levels_srmr)) %>%
#   ggplot(aes(y = condition, x = fp_srmr, color = as.factor(model))) +
#   geom_point() + geom_line() +
#   theme(axis.text.x = element_text(angle = 90))
# 
# levels_rmsea <- study2_all_fp1 %>% arrange(fp_rmsea) %>%  select(condition) %>% pull()
# 
# study2_all_fp1 %>% 
#   mutate(condition = factor(condition, levels = levels_rmsea)) %>%
#   ggplot(aes(y = condition, x = fp_rmsea, color = as.factor(model))) +
#   geom_point() + geom_line() +
#   theme(axis.text.x = element_text(angle = 90))
# 
# ## Figures for FP - All Correlations - Method 2
# levels_cfi <- study2_all_fp2 %>% group_by(conv, mlit, ident, bounds) %>% arrange(fp_cfi, .by_group = T) %>% ungroup() %>%  select(condition) %>% pull()


# study2_all_fp2 %>% 
#   mutate(condition = factor(condition, levels = levels_cfi)) %>%
#   ggplot(aes(y = condition, x = fp_cfi, color = as.factor(model))) +
#   geom_point() + geom_line() +
#   theme(axis.text.x = element_text(angle = 90))
# 
# levels_srmr <- study2_all_fp2 %>% group_by(conv, mlit, ident, bounds) %>% arrange(fp_srmr, .by_group = T) %>% ungroup() %>%  select(condition) %>% pull()
# 
# 
# study2_all_fp2 %>% 
#   mutate(condition = factor(condition, levels = levels_srmr)) %>%
#   ggplot(aes(y = condition, x = fp_srmr, color = as.factor(model))) +
#   geom_point() + geom_line() +
#   theme(axis.text.x = element_text(angle = 90))
# 
# levels_rmsea <- study2_all_fp2 %>% arrange(fp_rmsea) %>%  select(condition) %>% pull()
# 
# study2_all_fp2 %>% 
#   mutate(condition = factor(condition, levels = levels_rmsea)) %>%
#   ggplot(aes(y = condition, x = fp_rmsea, color = as.factor(model))) +
#   geom_point() + geom_line() +
#   theme(axis.text.x = element_text(angle = 90))

## Maybe try FP based on logl?
 
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
save(study2_all_logl_fp, file = here::here("Online-Supplement","Data-Files" ,"study2_all_logl_fp_badness.RData"))


study2_pos_logl_fp <- study2_pos_total %>% 
  #drop_na() %>%
  mutate(logl = if_else(is.na(logl), -40000, logl)) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -27103) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()

# saving the results
save(study2_pos_logl_fp, file = here::here("Online-Supplement","Data-Files", "study2_pos_logl_fp_badness.RData"))
#-------------------------------------------------------------------------------
# without replacing the NAs
study2_all_logl_fp_1 <- study2_all_total %>% 
  drop_na() %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -32755) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()

# saving the results
save(study2_all_logl_fp_1, file = here::here("Online-Supplement","Data-Files", "study2_all_logl_fp_true.RData"))


study2_pos_logl_fp_1 <- study2_pos_total %>% 
  drop_na() %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -27103) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()

# saving the results
save(study2_pos_logl_fp_1, file = here::here("Online-Supplement","Data-Files", "study2_pos_logl_fp_true.RData"))
#-------------------------------------------------------------------------------

# levels_logl_fp <- study2_all_logl_fp %>% group_by(conv, mlit, ident, bounds) %>% arrange(fp_logl, .by_group = TRUE) %>%  ungroup() %>% select(condition) %>% pull()
# 
# study2_all_logl_fp %>% 
#   mutate(condition = factor(condition, levels = levels_logl_fp)) %>%
#   ggplot(aes(y = condition, x = fp_logl, color = as.factor(model))) +
#   geom_point() + geom_line() +
#   theme(axis.text.x = element_text(angle = 90))
# study2_all_logl_fp %>% arrange(fp_logl)
# 
# levels_logl_fp <- study2_pos_logl_fp %>% group_by(conv, mlit, ident, bounds) %>% arrange(fp_logl, .by_group = TRUE) %>%  ungroup() %>% select(condition) %>% pull()
# 
# study2_pos_logl_fp %>% 
#   mutate(condition = factor(condition, levels = levels_logl_fp)) %>%
#   ggplot(aes(y = condition, x = fp_logl, color = as.factor(model))) +
#   geom_point() + geom_line() +
#   theme(axis.text.x = element_text(angle = 90))

## Report point estimates first to show that pattern of FP is same regardless of whether you force convergence or not
## Then show how stable estimates are across data space sizes with forced convergence vs not. Look at how many 
## 

# removing the files in the environment
rm(list = ls())

## Same for Study 3 ------------------------------------------------------------

# loading the data set
load(here::here("Online-Supplement","Data-Files","output_all_3_conditions.RData"))

study3_all_fits <- output_conditions

# loading the data set
load(here::here("Online-Supplement","Data-Files","output_poscor_3_all_conditions.RData"))

study3_pos_fits <- output_conditions

rm(output_conditions)

## Combine all output in one dataframe ----
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

# Study 3 conditions ----
study3_conditions <- expand.grid(model = 1:3,
                                 bounds = c("none", "pos.var"),
                                 ident = c("refmarker", "unitvar"),
                                 mlit = c('150', '1500'),
                                 conv = c("F", "T"))
study3_conditions$condition <- 1:48

## FOR CUTOFFS STUDY 3 ----
study3_pos_total %>% group_by(model) %>% summarize(m = mean(logl, na.rm = T))

# model       m
# <fct>   <dbl>
# 1 1     -10306.
# 2 2     -30479.
# 3 3     -50076.

study3_all_total %>% group_by(model) %>% summarize(m = mean(logl, na.rm = T))
# 
# model       m
# <fct>   <dbl>
# 1 1     -11460.
# 2 2     -36999.
# 3 3     -62493.

# saving the results
save(study3_pos_total, file = here::here("Online-Supplement","Data-Files" ,"study3_pos_total.RData"))

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
save(study3_all_avg, file = here::here("Online-Supplement","Data-Files", "study3_all_averages.RData"))

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

## Maybe try FP based on logl?

#-------------------------------------------------------------------------------
# here we have replaced NAs with assumed badness of fit
study3_all_logl_fp <- study3_all_total %>% 
  #drop_na() %>%
  mutate(logl = if_else(is.na(logl), -40000, logl)) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -32755) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()

# saving the results
save(study3_all_logl_fp, file = here::here("Online-Supplement","Data-Files", "study3_all_logl_fp_badness.RData"))

study3_pos_logl_fp <- study3_pos_total %>% 
  #drop_na() %>%
  mutate(logl = if_else(is.na(logl), -40000, logl)) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -27103) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()

# saving the results
save(study3_pos_logl_fp, file = here::here("Online-Supplement","Data-Files", "study3_pos_logl_fp_badness.RData"))
#-------------------------------------------------------------------------------
# without replacing the NAs
study3_all_logl_fp_1 <- study3_all_total %>% 
  drop_na() %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -32755) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()

# saving the results
save(study3_all_logl_fp_1, file = here::here("Online-Supplement","Data-Files", "study3_all_logl_fp_true.RData"))


study3_pos_logl_fp_1 <- study3_pos_total %>% 
  drop_na() %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -27103) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
  ungroup()

# saving the results
save(study3_pos_logl_fp_1, file = here::here("Online-Supplement","Data-Files", "study3_pos_logl_fp_true.RData"))
#-------------------------------------------------------------------------------


# study3_pos_total %>%
#   mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
#   ggplot(aes(y = condition, x = logl)) +
#   geom_boxplot() +
#   facet_wrap(vars(model), nrow = 3, ncol = 1, scales = "free_y")
# 
# study3_all_total %>%
#   mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
#   ggplot(aes(y = condition, x = logl)) +
#   geom_boxplot() +
#   facet_wrap(vars(model), nrow = 3, ncol = 1, scales = "free_y")
# 
# study3_pos_avg <- 
#   study3_pos_total %>% group_by(model, conv, mlit, ident, bounds) %>%
#   summarize(m_cfi = mean(cfi, na.rm = T),
#             m_srmr = mean(srmr, na.rm = T),
#             m_rmsea = mean(rmsea, na.rm = T),
#             m_logl = mean(logl, na.rm = T),
#             condition_num = min(condition)) %>%
#   mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
#   ungroup()
# 
# study3_pos_avg2 <- study3_pos_avg %>%
#   select(m_cfi:condition_num) %>%
#   mutate(condition = condition_num) %>%
#   bind_cols(study3_conditions) %>%
#   mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
#   ungroup()
# 
# 
# study3_pos_fp <- 
#   study3_pos_total %>% 
#   mutate(cfi = if_else(is.na(cfi), 0, cfi),
#          rmsea = if_else(is.na(rmsea), 1, rmsea),
#          srmr = if_else(is.na(srmr), 1, srmr)) %>%
#   group_by(model, conv, mlit, ident, bounds) %>%
#   summarize(fp_cfi = (sum(cfi > .45) / n()),
#             fp_srmr = (sum(srmr < .09) / n()),
#             fp_rmsea = (sum(rmsea < .27) / n()),
#             condition_num = min(condition)) %>%
#   mutate(condition = paste(model, conv, mlit, ident, bounds, sep = "_")) %>%
#   ungroup()
# 
# ## Figures for average fit values
# levels_cfi <- study3_pos_avg2 %>% arrange(m_cfi) %>%  select(condition) %>% pull()
# levels_srmr <- study3_pos_avg %>% arrange(m_srmr) %>%  select(condition) %>% pull()
# levels_rmsea <- study3_pos_avg %>% arrange(m_rmsea) %>%  select(condition) %>% pull()
# levels_logl <- study3_pos_avg %>% arrange(m_logl) %>%  select(condition) %>% pull()
# 
# study3_pos_avg2 %>% 
#   mutate(condition = factor(condition, levels = levels_cfi)) %>%
#   ggplot(aes(y = condition, x = m_cfi, color = as.factor(model))) +
#   geom_point() + geom_line() +
#   facet_wrap(vars(model), nrow = 3, ncol = 1, scales = "free") +
#   theme(axis.text.x = element_text(angle = 90))
# 
# study3_pos_avg %>% 
#   mutate(condition = factor(condition, levels = levels_srmr)) %>%
#   ggplot(aes(y = condition, x = m_srmr, color = as.factor(model))) +
#   geom_point() + geom_line() +
#   facet_wrap(vars(model), nrow = 3, ncol = 1, scales = "free") +
#   theme(axis.text.x = element_text(angle = 90))
# 
# study3_pos_avg %>% 
#   mutate(condition = factor(condition, levels = levels_rmsea)) %>%
#   ggplot(aes(y = condition, x = m_rmsea, color = as.factor(model))) +
#   geom_point() + geom_line() +
#   facet_wrap(vars(model), nrow = 3, ncol = 1, scales = "free") +
#   theme(axis.text.x = element_text(angle = 90))
# 
# study3_pos_avg %>% 
#   mutate(condition = factor(condition, levels = levels_logl)) %>%
#   ggplot(aes(y = condition, x = m_logl, color = as.factor(model))) +
#   geom_point() + geom_line() +
#   facet_wrap(vars(model), nrow = 3, ncol = 1, scales = "free") +
#   theme(axis.text.x = element_text(angle = 90))
# 
