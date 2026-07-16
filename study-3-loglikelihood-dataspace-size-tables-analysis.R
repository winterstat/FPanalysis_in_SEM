#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 09/25/2025
#-------------------------------------------------------------------------------
# Relationship between the data-space size and the fitting propensity of the models
# in study 3
#-------------------------------------------------------------------------------

## Load packages ---------------------------------------------------------------
pkgs <- c('tidyverse', 'ggdist', 'ggthemes', 'effectsize')

# check if packages are available
rlang::check_installed(pkgs)

# loading packages
lapply(pkgs, library, character.only = TRUE)

# Case I: No-listwise Analysis--------------------------------------------------------------

# Loading the tables for this case

# descriptive statistics
table_11 <- readRDS(here::here("Online-Supplement", "Data-Files", "table_11_badness.RDS"))

table_22 <- readRDS(here::here("Online-Supplement", "Data-Files", "table_22_badness.RDS"))

table_33 <- readRDS(here::here("Online-Supplement","Data-Files", "table_33_badness.RDS"))

# totals
total_11 <- readRDS(here::here("Online-Supplement","Data-Files", "total_11_badness.RDS"))

total_22 <- readRDS(here::here("Online-Supplement","Data-Files", "total_22_badness.RDS"))

total_33 <- readRDS(here::here("Online-Supplement","Data-Files", "total_33_badness.RDS"))
#-------------------------------------------------------------------------------

# Case II: listwise Analysis----------------------------------------------------

# Descriptive statistics
table_111 <- readRDS(here::here("Online-Supplement", "Data-Files", "table_111_true.RDS"))

table_222 <-  readRDS(here::here("Online-Supplement", "Data-Files", "table_222_true.RDS"))

table_333 <-  readRDS(here::here("Online-Supplement", "Data-Files", "table_333_true.RDS"))

# totals
total_111 <- readRDS(here::here("Online-Supplement", "Data-Files", "total_111_true.RDS"))

total_222 <- readRDS(here::here("Online-Supplement", "Data-Files", "total_111_true.RDS"))

total_333 <-  readRDS(here::here("Online-Supplement", "Data-Files", "total_333_true.RDS"))

#-------------------------------------------------------------------------------
# Model 1:

# no-listwise
table_11 %>%
  select(contains('1'), bounds, ident) %>%
  relocate(all_1_F_150, all_1_F_1500)

# all_1_F_150   all_1_F_1500  all_1_T_150   all_1_T_1500  bounds  ident    
# <chr>         <chr>         <chr>         <chr>         <fct>   <fct>    
# 1 0.131 (0.003) 0.211 (0.004) 0.632 (0.005) 0.655 (0.005) none    refmarker
# 2 0.062 (0.002) 0.094 (0.003) 0.621 (0.005) 0.633 (0.005) none    unitvar  
# 3 0.208 (0.003) 0.266 (0.004) 0.323 (0.004) 0.342 (0.005) pos.var refmarker
# 4 0.121 (0.003) 0.14 (0.003)  0.303 (0.005) 0.308 (0.005) pos.var unitvar  

# listwise
table_111 %>%
  select(contains('1'), bounds, ident) %>%
  relocate(all_1_F_150, all_1_F_1500)

# all_1_F_150   all_1_F_1500  all_1_T_150   all_1_T_1500  bounds  ident    
# <chr>         <chr>         <chr>         <chr>         <fct>   <fct>    
# 1 0.577 (0.01)  0.656 (0.009) 0.633 (0.005) 0.656 (0.005) none    refmarker
# 2 0.519 (0.012) 0.558 (0.01)  0.621 (0.005) 0.634 (0.005) none    unitvar  
# 3 0.27 (0.005)  0.315 (0.005) 0.324 (0.005) 0.343 (0.005) pos.var refmarker
# 4 0.189 (0.004) 0.209 (0.005) 0.303 (0.004) 0.309 (0.004) pos.var unitvar  

# Model 2:----------------------------------------------------------------------

# no-listwise
table_22 %>%
  select(contains('2'), bounds, ident) %>%
  relocate(all_2_F_150, all_2_F_1500)

# A tibble: 4 × 6
# all_2_F_150   all_2_F_1500  all_2_T_150   all_2_T_1500  bounds  ident    
# <chr>         <chr>         <chr>         <chr>         <fct>   <fct>    
# 1 0.013 (0.001) 0.04 (0.002)  0.694 (0.005) 0.765 (0.004) none    refmarker
# 2 0.008 (0.001) 0.012 (0.001) 0.501 (0.005) 0.537 (0.005) none    unitvar  
# 3 0.026 (0.001) 0.059 (0.002) 0.32 (0.005)  0.368 (0.005) pos.var refmarker
# 4 0.015 (0.001) 0.018 (0.002) 0.167 (0.004) 0.177 (0.004) pos.var unitvar  

# listwise
table_222 %>%
  select(contains('2'), bounds, ident) %>%
  relocate(all_2_F_150, all_2_F_1500)
# 
# all_2_F_150   all_2_F_1500  all_2_T_150   all_2_T_1500  bounds  ident    
# <chr>         <chr>         <chr>         <chr>         <fct>   <fct>    
# 1 0.197 (0.018) 0.343 (0.015) 0.695 (0.005) 0.766 (0.005) none    refmarker
# 2 0.07 (0.008)  0.089 (0.009) 0.501 (0.005) 0.537 (0.005) none    unitvar  
# 3 0.06 (0.004)  0.106 (0.004) 0.32 (0.005)  0.367 (0.005) pos.var refmarker
# 4 0.019 (0.002) 0.023 (0.002) 0.168 (0.004) 0.177 (0.004) pos.var unitvar

# Model 3:----------------------------------------------------------------------

# no-listwise
table_33 %>%
  select(contains('3'), bounds, ident) %>%
  relocate(all_3_F_150, all_3_F_1500)

# all_3_F_150 all_3_F_1500  all_3_T_150   all_3_T_1500  bounds  ident    
# <chr>       <chr>         <chr>         <chr>         <fct>   <fct>    
# 1 0.001 (0)   0.006 (0.001) 0.691 (0.005) 0.77 (0.005)  none    refmarker
# 2 0 (0)       0.001 (0)     0.413 (0.005) 0.463 (0.005) none    unitvar  
# 3 0.002 (0)   0.011 (0.001) 0.254 (0.005) 0.422 (0.005) pos.var refmarker
# 4 0.001 (0)   0.001 (0)     0.176 (0.004) 0.203 (0.004) pos.var unitvar  

# listwise
table_333 %>%
  select(contains('3'), bounds, ident) %>%
  relocate(all_3_F_150, all_3_F_1500)

# all_3_F_150   all_3_F_1500  all_3_T_150   all_3_T_1500  bounds  ident    
# <chr>         <chr>         <chr>         <chr>         <fct>   <fct>    
# 1 0.016 (0.006) 0.061 (0.008) 0.69 (0.004)  0.769 (0.004) none    refmarker
# 2 0.002 (0.001) 0.004 (0.001) 0.414 (0.005) 0.464 (0.005) none    unitvar  
# 3 0.006 (0.002) 0.025 (0.002) 0.254 (0.004) 0.422 (0.005) pos.var refmarker
# 4 0.001 (0)     0.002 (0)     0.176 (0.004) 0.203 (0.004) pos.var unitvar  
