#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 01/23/2026
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

# loading data
# Analysis:
# loading the datasets
res_pos_model_1_true <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_1_true_dataspacesize.RDS"))
res_all_model_1_true <- readRDS(here::here("Online-Supplement","Data-Files","results_allcor_study_3_model_1_true_dataspacesize.RDS"))

res_pos_model_2_true <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_2_true_dataspacesize.RDS"))
res_all_model_2_true <- readRDS(here::here("Online-Supplement","Data-Files","results_allcor_study_3_model_2_true_dataspacesize.RDS"))

res_pos_model_3_true <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_3_true_dataspacesize.RDS"))
res_all_model_3_true <- readRDS(here::here("Online-Supplement","Data-Files","results_allcor_study_3_model_3_true_dataspacesize.RDS"))
#-------------------------------------------------------------------------------


# new variable in each data frame-----------------------------------------------
# Model 1
results_full_1_true <- bind_rows(
  
  res_pos_model_1_true %>% mutate(cor = 'pos') %>% mutate(cor = factor(cor)),
  
  res_all_model_1_true %>% mutate(cor = 'all') %>% mutate(cor = factor(cor))
)

# Model 2
results_full_2_true <- bind_rows(
  
  res_pos_model_2_true %>% mutate(cor = 'pos') %>% mutate(cor = factor(cor)),
  
  res_all_model_2_true %>% mutate(cor = 'all') %>% mutate(cor = factor(cor))
)

# Model 3
results_full_3_true <- bind_rows(
  
  res_pos_model_3_true %>% mutate(cor = 'pos') %>% mutate(cor = factor(cor)),
  
  res_all_model_3_true %>% mutate(cor = 'all') %>% mutate(cor = factor(cor))
)


# loading the datasets
res_pos_model_1_bad <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_1_dataspacesize_badness.RDS"))
res_all_model_1_bad <- readRDS(here::here("Online-Supplement","Data-Files", "results_allcor_study_3_model_1_dataspacesize_badness.RDS"))

res_pos_model_2_bad <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_2_dataspacesize_badness.RDS"))
res_all_model_2_bad <- readRDS(here::here("Online-Supplement","Data-Files", "results_allcor_study_3_model_2_dataspacesize_badness.RDS"))

res_pos_model_3_bad <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_study_3_model_3_dataspacesize_badness.RDS"))
res_all_model_3_bad <- readRDS(here::here("Online-Supplement","Data-Files", "results_allcor_study_3_model_3_dataspacesize_badness.RDS"))
#-------------------------------------------------------------------------------


# new variable in each data frame-----------------------------------------------
# Model 1
results_full_1_bad <- bind_rows(
  
  res_pos_model_1_bad %>% mutate(cor = 'pos') %>% mutate(cor = factor(cor)),
  
  res_all_model_1_bad %>% mutate(cor = 'all') %>% mutate(cor = factor(cor))
)

# Model 2
results_full_2_bad <- bind_rows(
  
  res_pos_model_2_bad %>% mutate(cor = 'pos') %>% mutate(cor = factor(cor)),
  
  res_all_model_2_bad %>% mutate(cor = 'all') %>% mutate(cor = factor(cor))
)

# Model 3
results_full_3_bad <- bind_rows(
  
  res_pos_model_3_bad %>% mutate(cor = 'pos') %>% mutate(cor = factor(cor)),
  
  res_all_model_3_bad %>% mutate(cor = 'all') %>% mutate(cor = factor(cor))
)


# Fitting propensity------------------------------------------------------------

# Model 1
results_full_1_true %>%
  select(-cfi, -srmr, -rmsea) %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size) %>%
  summarize(m = mean(value),
            se = sd(value), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>% print(n = 26)

results_full_1_bad %>%
  select(-cfi, -srmr, -rmsea) %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size) %>%
  summarize(m = mean(value),
            se = sd(value), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>% print(n = 26)

# Total
results_full_1_true %>%
select(-cfi, -srmr, -rmsea) %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(cor) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop") %>% print(n = 26)


results_full_1_bad %>%
  select(-cfi, -srmr, -rmsea) %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(cor) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop") %>% print(n = 26)

#--------------------------------------------------------------------

# Model 2
results_full_2_true %>%
  select(-cfi, -srmr, -rmsea) %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size) %>%
  summarize(m = mean(value, na.rm = TRUE),
            se = sd(value, na.rm = TRUE), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>% print(n = 26)

results_full_2_bad %>%
  select(-cfi, -srmr, -rmsea) %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size) %>%
  summarize(m = mean(value),
            se = sd(value), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>% print(n = 26)

# Total
results_full_2_true %>%
  select(-cfi, -srmr, -rmsea) %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(cor) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop") %>% print(n = 26)


results_full_2_bad %>%
  select(-cfi, -srmr, -rmsea) %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(cor) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop") %>% print(n = 26)


#--------------------------------------------------------------------

# Model 3
results_full_3_true %>%
  select(-cfi, -srmr, -rmsea) %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size) %>%
  summarize(m = mean(value, na.rm = TRUE),
  se = sd(value, na.rm = TRUE), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>% print(n = 26)

results_full_3_bad %>%
  select(-cfi, -srmr, -rmsea) %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size) %>%
  summarize(m = mean(value),
  se = sd(value), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>% print(n = 26)

# Total
results_full_3_true %>%
  select(-cfi, -srmr, -rmsea) %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(cor) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop") %>% print(n = 26)


results_full_3_bad %>%
  select(-cfi, -srmr, -rmsea) %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(cor) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop") %>% print(n = 26)


