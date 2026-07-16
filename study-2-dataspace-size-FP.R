#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 01/23/2026
#-------------------------------------------------------------------------------
# Relationship between the data space and the fitting propensity of the models
# in study 2
#-------------------------------------------------------------------------------

## Load packages ---------------------------------------------------------------
pkgs <- c('tidyverse', 'ggdist', 'ggthemes', 'effectsize', 'dplyr', 'ggplot2')

# check if packages are available
# rlang::check_installed(pkgs)

# loading packages
lapply(pkgs, library, character.only = TRUE)

# loading data

# reading and merging the datasets----------------------------------------------
res_pos_bad <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_badness_study_2_dataspacesize.RDS"))

res_all_bad <- readRDS(here::here("Online-Supplement","Data-Files", "results_all_badness_study_2_dataspacesize.RDS"))

# new variable in each data frame-----------------------------------------------
results_full_bad <- bind_rows(
  
  res_pos_bad %>% mutate(cor = 'pos') %>% mutate(cor = factor(cor)),
  
  res_all_bad %>% mutate(cor = 'all') %>% mutate(cor = factor(cor))
)


# reading and merging the datasets----------------------------------------------
res_pos_true <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_true_study_2_dataspacesize.RDS"))

res_all_true <- readRDS(here::here("Online-Supplement","Data-Files", "results_all_true_study_2_dataspacesize.RDS"))

# new variable in each data frame-----------------------------------------------
results_full_true <- bind_rows(
  
  res_pos_true %>% mutate(cor = 'pos') %>% mutate(cor = factor(cor)),
  
  res_all_true %>% mutate(cor = 'all') %>% mutate(cor = factor(cor))
)

# Fitting propensity------------------------------------------------------------

# Model 1
results_full_true %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(model == '1') %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size) %>%
  summarize(m = mean(value),
            se = sd(value), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>% print( n = 26)

results_full_bad %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(model == '1') %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size) %>%
  summarize(m = mean(value),
            se = sd(value), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>% print(n = 26)


# Total
results_full_true %>%
  filter(model == '1') %>%
  select(-cfi, -srmr, -rmsea) %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(cor) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop") %>% print(n = 26)


results_full_bad %>%
  filter(model == '1') %>%
  select(-cfi, -srmr, -rmsea) %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(cor) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop") %>% print(n = 26)

#------------------------------------------------------------------------------

# Model 2

results_full_true %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(model == '2') %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size) %>%
  summarize(m = mean(value, na.rm = TRUE),
            se = sd(value, na.rm = TRUE), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>% print( n = 26)

results_full_bad %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(model == '2') %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size) %>%
  summarize(m = mean(value),
            se = sd(value), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>% print(n = 26)


# Total
results_full_true %>%
  filter(model == '2') %>%
  select(-cfi, -srmr, -rmsea) %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(cor) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop") %>% print(n = 26)


results_full_bad %>%
  filter(model == '2') %>%
  select(-cfi, -srmr, -rmsea) %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(cor) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop") %>% print(n = 26)
#-------------------------------------------------------------------------------

# Model 3

results_full_true %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(model == '3') %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size) %>%
  summarize(m = mean(value, na.rm = TRUE),
            se = sd(value, na.rm = TRUE), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>% print( n = 26)

results_full_bad %>%
  select(-cfi, -srmr, -rmsea) %>%
  filter(model == '3') %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(model, cor, size) %>%
  summarize(m = mean(value, na.rm = TRUE),
            se = sd(value, na.rm = TRUE), .groups = "drop")%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>% print(n = 26)


# Total
results_full_true %>%
  filter(model == '3') %>%
  select(-cfi, -srmr, -rmsea) %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(cor) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop") %>% print(n = 26)


results_full_bad %>%
  filter(model == '3') %>%
  select(-cfi, -srmr, -rmsea) %>%
  pivot_longer(cols = logl, names_to = "index") %>%
  group_by(cor) %>%
  summarize(total = mean(value, na.rm = T), .groups = "drop") %>% print(n = 26)
