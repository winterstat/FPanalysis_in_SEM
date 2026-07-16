#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 09/09/2025
#-------------------------------------------------------------------------------
# Fitting propensity analysis for study 2
#-------------------------------------------------------------------------------

## Load packages ---------------------------------------------------------------
pkgs <- c('tidyverse', 'ggdist', 'ggthemes', 'effectsize', 'gt')

# check if packages are available
rlang::check_installed(pkgs)

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


# loading the datasets---------------------------------------------------------
load(file = here::here("Online-Supplement","Data-Files", "study2_pos_total.RData"))

load(file = here::here("Online-Supplement","Data-Files", "study2_all_total.RData"))

# Fitting propensity analysis in study 2----------------------------------------
# here we have replaced NAs with assumed badness of fit
study2_all_logl_fp <- study2_all_total %>% 
  mutate(logl = if_else(is.na(logl), -35000, logl)) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -32755) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(fp_logl = fp_logl*100) 

# saving the results
save(study2_all_logl_fp, file = here::here("Online-Supplement","Data-Files", "study2_all_logl_fp_badness.RData"))

# creating a table--------------------------------------------------------------
study2_all_logl_fp %>% # removing the frequency column
  pivot_wider(names_from = c(conv, mlit),
              values_from = fp_logl, names_expand = F, id_cols = c(model, ident, bounds)) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  ) -> study2_all_logl_fp_wider

#------------------------------------------------------------------------------- 
# the actual names
actual_names <- colnames(study2_all_logl_fp_wider)

# new names
new_names <- actual_names %>%
  str_remove('(F|T)_') %>%
  str_to_title()

names(new_names) <- actual_names

# Table using gt----------------------------------------------------------------
library(gt)
#-------------------------------------------------------------------------------
study2_all_logl_fp_wider %>%
  gt(groupname_col = 'model', rowname_col = 'ident') %>% 
  cols_label(.list = new_names) %>%
  tab_spanner(
    label = md('T'),
    columns = contains('T')
  ) %>%
  tab_spanner(
    label = md('F'),
    columns = contains('F')
  ) %>%
  tab_spanner(
    label = md('Forced Convergence'),
    columns = contains('T') | contains('F')
  ) %>%
  cols_align(
    align = 'left',
    columns = c(where(is.factor), where(is.character))
  ) %>%
  cols_label(
    bounds = 'Maximum Iterations'  
  )%>%
  fmt_number(decimals = 3) -> study2_all_logl_fp_table

# saving the table
gt::gtsave(data = study2_all_logl_fp_table, 
filename = "study2_all_badness.docx", path=here::here("Online-Supplement","All-Correlations-Tables", "study-2"))
#-------------------------------------------------------------------------------
study2_pos_logl_fp <- study2_pos_total %>% 
  mutate(logl = if_else(is.na(logl), -35000, logl)) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -27103) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(fp_logl = fp_logl*100) 


# saving the results
save(study2_pos_logl_fp, file = here::here("Online-Supplement","Data-Files", "study2_pos_logl_fp_badness.RData"))
#-------------------------------------------------------------------------------

# creating a table--------------------------------------------------------------
study2_pos_logl_fp %>% # removing the frequency column
  pivot_wider(names_from = c(conv, mlit),
              values_from = fp_logl, names_expand = F, id_cols = c(model, ident, bounds)) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  ) -> study2_pos_logl_fp_wider

#------------------------------------------------------------------------------- 
# the actual names
actual_names <- colnames(study2_pos_logl_fp_wider)

# new names
new_names <- actual_names %>%
  str_remove('(F|T)_') %>%
  str_to_title()

names(new_names) <- actual_names
#-------------------------------------------------------------------------------
study2_pos_logl_fp_wider %>%
  gt(groupname_col = 'model', rowname_col = 'ident') %>% 
  cols_label(.list = new_names) %>%
  tab_spanner(
    label = md('T'),
    columns = contains('T')
  ) %>%
  tab_spanner(
    label = md('F'),
    columns = contains('F')
  ) %>%
  tab_spanner(
    label = md('Forced Convergence'),
    columns = contains('T') | contains('F')
  ) %>%
  cols_align(
    align = 'left',
    columns = c(where(is.factor), where(is.character))
  ) %>%
  cols_label(
    bounds = 'Maximum Iterations'  
  )%>%
  fmt_number(decimals = 3) -> study2_pos_logl_fp_table

# saving the table
gt::gtsave(data = study2_pos_logl_fp_table, 
filename = "study2_pos_badness.docx", path= here::here("Online-Supplement","Positive-Correlations-Tables", "study-2"))


#-------------------------------------------------------------------------------
# not replacing the NAs with badness of fit
#-------------------------------------------------------------------------------
study2_all_logl_fp_true <- study2_all_total %>% 
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = mean(logl > -32755, na.rm = TRUE), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(fp_logl = fp_logl*100) 

# saving the results
save(study2_all_logl_fp_true, file = here::here("Online-Supplement","Data-Files","study2_all_logl_fp_true.RData"))

# creating a table--------------------------------------------------------------
study2_all_logl_fp_true %>% # removing the frequency column
  pivot_wider(names_from = c(conv, mlit),
              values_from = fp_logl, names_expand = F, id_cols = c(model, ident, bounds)) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  ) -> study2_all_logl_fp_wider_true

#------------------------------------------------------------------------------- 
# the actual names
actual_names <- colnames(study2_all_logl_fp_wider_true)

# new names
new_names <- actual_names %>%
  str_remove('(F|T)_') %>%
  str_to_title()

names(new_names) <- actual_names

#-------------------------------------------------------------------------------
study2_all_logl_fp_wider_true %>%
  gt(groupname_col = 'model', rowname_col = 'ident') %>% 
  cols_label(.list = new_names) %>%
  tab_spanner(
    label = md('T'),
    columns = contains('T')
  ) %>%
  tab_spanner(
    label = md('F'),
    columns = contains('F')
  ) %>%
  tab_spanner(
    label = md('Forced Convergence'),
    columns = contains('T') | contains('F')
  ) %>%
  cols_align(
    align = 'left',
    columns = c(where(is.factor), where(is.character))
  ) %>%
  cols_label(
    bounds = 'Maximum Iterations'  
  )%>%
  fmt_number(decimals = 3) -> study2_all_logl_fp_table_true

# saving the table
gt::gtsave(data = study2_all_logl_fp_table_true, 
filename = "study2_all_true.docx", path= here::here("Online-Supplement","All-Correlations-Tables", "study-2"))


#-------------------------------------------------------------------------------
study2_pos_logl_fp_true <- study2_pos_total %>% 
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = mean(logl > -27103, na.rm = TRUE), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(fp_logl = fp_logl*100) 


# saving the results
save(study2_pos_logl_fp_true, file = here::here("Online-Supplement","Data-Files","study2_pos_logl_fp_true.RData"))
#-------------------------------------------------------------------------------

# creating a table--------------------------------------------------------------
study2_pos_logl_fp_true %>% 
  pivot_wider(names_from = c(conv, mlit),
              values_from = fp_logl, names_expand = F, id_cols = c(model, ident, bounds)) %>%
  relocate(model) %>% 
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  ) -> study2_pos_logl_fp_wider_true

#------------------------------------------------------------------------------- 
# the actual names
actual_names <- colnames(study2_pos_logl_fp_wider_true)

# new names
new_names <- actual_names %>%
  str_remove('(F|T)_') %>%
  str_to_title()

names(new_names) <- actual_names
#-------------------------------------------------------------------------------
study2_pos_logl_fp_wider %>%
  gt(groupname_col = 'model', rowname_col = 'ident') %>% 
  cols_label(.list = new_names) %>%
  tab_spanner(
    label = md('T'),
    columns = contains('T')
  ) %>%
  tab_spanner(
    label = md('F'),
    columns = contains('F')
  ) %>%
  tab_spanner(
    label = md('Forced Convergence'),
    columns = contains('T') | contains('F')
  ) %>%
  cols_align(
    align = 'left',
    columns = c(where(is.factor), where(is.character))
  ) %>%
  cols_label(
    bounds = 'Maximum Iterations'  
  )%>%
  fmt_number(decimals = 3) -> study2_pos_logl_fp_table_true

# saving the table
gt::gtsave(data = study2_pos_logl_fp_table_true, 
  filename = "study2_pos_true.docx", path= here::here("Online-Supplement","Positive-Correlations-Tables", "study-2"))
