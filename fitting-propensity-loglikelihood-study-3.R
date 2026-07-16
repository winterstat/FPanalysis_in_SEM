#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 09/09/2025
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Fitting propensity analysis for study 3
#-------------------------------------------------------------------------------

## Load packages ---------------------------------------------------------------
pkgs <- c('tidyverse', 'ggdist', 'ggthemes', 'effectsize')

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
load(file = here::here("Online-Supplement","Data-Files", "study3_pos_total.RData"))

load(file = here::here("Online-Supplement","Data-Files", "study3_all_total.RData"))
#-------------------------------------------------------------------------------
# here we have replaced NAs with assumed badness of fit
#-------------------------------------------------------------------------------

study3_all_logl_fp_model_1 <- study3_all_total %>% 
  filter(model == 1) %>%
  mutate(logl = if_else(is.na(logl), -20000, logl)) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -11460) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(fp_logl = fp_logl*100) 

study3_all_logl_fp_model_2 <- study3_all_total %>% 
  filter(model == 2) %>%
  mutate(logl = if_else(is.na(logl), -40000, logl)) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -36999) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(fp_logl = fp_logl*100) 


study3_all_logl_fp_model_3 <- study3_all_total %>% 
  filter(model == 3) %>%
  mutate(logl = if_else(is.na(logl), -70000, logl)) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -62493) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(fp_logl = fp_logl*100) 


# saving the results
save(study3_all_logl_fp_model_1, file = here::here("Online-Supplement","Data-Files","study3_all_logl_fp_mode1_1_badness.RData"))

save(study3_all_logl_fp_model_2, file = here::here("Online-Supplement","Data-Files", "study3_all_logl_fp_mode1_2_badness.RData"))

save(study3_all_logl_fp_model_3, file = here::here("Online-Supplement","Data-Files", "study3_all_logl_fp_mode1_3_badness.RData"))
#-----------------------------------------------------------------------------------------

study3_pos_logl_fp_model_1 <- study3_pos_total %>% 
  filter(model == 1) %>%
  mutate(logl = if_else(is.na(logl), -15000, logl)) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -10306) / n()), # CUTOFF
            condition_num = min(condition))  %>%
  mutate(fp_logl = fp_logl*100) 


study3_pos_logl_fp_model_2 <- study3_pos_total %>% 
  filter(model == 2) %>%
  mutate(logl = if_else(is.na(logl), -40000, logl)) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -30479) / n()), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(fp_logl = fp_logl*100) 

study3_pos_logl_fp_model_3 <- study3_pos_total %>% 
  filter(model == 3) %>%
  mutate(logl = if_else(is.na(logl), -60000, logl)) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = (sum(logl > -50076) / n()), # CUTOFF
            condition_num = min(condition))  %>%
  mutate(fp_logl = fp_logl*100) 


# saving results
save(study3_pos_logl_fp_model_1, file = here::here("Online-Supplement","Data-Files","study3_pos_logl_fp_mode1_1_badness.RData"))

save(study3_pos_logl_fp_model_2, file = here::here("Online-Supplement","Data-Files", "study3_pos_logl_fp_mode1_2_badness.RData"))

save(study3_pos_logl_fp_model_3, file = here::here("Online-Supplement","Data-Files", "study3_pos_logl_fp_mode1_3_badness.RData"))
#-------------------------------------------------------------------------------
# here we do not assign worst fit
#-------------------------------------------------------------------------------
study3_all_logl_fp_model_11_true <- study3_all_total %>% 
  filter(model == 1) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = mean(logl > -11460, na.rm = TRUE), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(fp_logl = fp_logl*100) 

study3_all_logl_fp_model_22_true <- study3_all_total %>% 
  filter(model == 2) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = mean(logl > -36999, na.rm = TRUE), # CUTOFF
            condition_num = min(condition))  %>%
  mutate(fp_logl = fp_logl*100) 


study3_all_logl_fp_model_33_true <- study3_all_total %>% 
  filter(model == 3) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = mean(logl > -62493, na.rm = TRUE), # CUTOFF
            condition_num = min(condition))  %>%
  mutate(fp_logl = fp_logl*100) 


# saving the results
save(study3_all_logl_fp_model_11_true, file = here::here("Online-Supplement","Data-files", "study3_all_logl_fp_mode1_1_true.RData"))
save(study3_all_logl_fp_model_22_true, file = here::here("Online-Supplement","Data-files", "study3_all_logl_fp_mode1_2_true.RData"))
save(study3_all_logl_fp_model_33_true, file = here::here("Online-Supplement","Data-files", "study3_all_logl_fp_mode1_3_true.RData"))
#-----------------------------------------------------------------------------------------

study3_pos_logl_fp_model_11_true <- study3_pos_total %>% 
  filter(model == 1) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = mean(logl > -10306, na.rm = TRUE), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(fp_logl = fp_logl*100) 


study3_pos_logl_fp_model_22_true<- study3_pos_total %>% 
  filter(model == 2) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = mean(logl > -30479, na.rm = TRUE), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(fp_logl = fp_logl*100) 


study3_pos_logl_fp_model_33_true <- study3_pos_total %>% 
  filter(model == 3) %>%
  group_by(model, conv, mlit, ident, bounds) %>%
  summarize(fp_logl = mean(logl > -50076, na.rm = TRUE), # CUTOFF
            condition_num = min(condition)) %>%
  mutate(fp_logl = fp_logl*100) 

# saving results
save(study3_pos_logl_fp_model_11_true, file = here::here("Online-Supplement","Data-files", "study3_pos_logl_fp_mode1_1_true.RData"))
save(study3_pos_logl_fp_model_22_true, file = here::here("Online-Supplement","Data-files", "study3_pos_logl_fp_mode1_2_true.RData"))
save(study3_pos_logl_fp_model_33_true, file = here::here("Online-Supplement","Data-files", "study3_pos_logl_fp_mode1_3_true.RData"))
#-------------------------------------------------------------------------------
# creating tables
# loading the gt package 
library(gt)

# Part 1: Badness of fit--------------------------------------------------------

# Model 1
study3_all_logl_fp_model_1 %>% # removing the frequency column
  pivot_wider(names_from = c(conv, mlit),
              values_from = fp_logl, names_expand = F, id_cols = c(model, ident, bounds)) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  ) -> study3_all_logl_fp_wider_model_1

#------------------------------------------------------------------------------- 
# the actual names
actual_names <- colnames(study3_all_logl_fp_wider_model_1)

# new names
new_names <- actual_names %>%
  str_remove('(F|T)_') %>%
  str_to_title()

names(new_names) <- actual_names

#-------------------------------------------------------------------------------
study3_all_logl_fp_wider_model_1 %>%
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
  fmt_number(decimals = 3) -> study3_all_logl_fp_table_model_1

# saving the table
gt::gtsave(data = study3_all_logl_fp_table_model_1, 
filename = 'study3_all_badness_model_1.docx', path = here::here("Online-Supplement","All-Correlations-Tables", "study-3"))
# Model 2------------------------------------------------------------------------
study3_all_logl_fp_model_2 %>% # removing the frequency column
  pivot_wider(names_from = c(conv, mlit),
              values_from = fp_logl, names_expand = F, id_cols = c(model, ident, bounds)) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  ) -> study3_all_logl_fp_wider_model_2

#------------------------------------------------------------------------------- 
# the actual names
actual_names <- colnames(study3_all_logl_fp_wider_model_2)

# new names
new_names <- actual_names %>%
  str_remove('(F|T)_') %>%
  str_to_title()

names(new_names) <- actual_names

#-------------------------------------------------------------------------------
study3_all_logl_fp_wider_model_2 %>%
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
  fmt_number(decimals = 3) -> study3_all_logl_fp_table_model_2

# saving the table
gt::gtsave(data = study3_all_logl_fp_table_model_2, 
    filename = 'study3_all_badness_model_2.docx', path=here::here("Online-Supplement","All-Correlations-Tables", "study-3"))


# Model 3------------------------------------------------------------------------
study3_all_logl_fp_model_3 %>% # removing the frequency column
  pivot_wider(names_from = c(conv, mlit),
              values_from = fp_logl, names_expand = F, id_cols = c(model, ident, bounds)) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  ) -> study3_all_logl_fp_wider_model_3

#------------------------------------------------------------------------------- 
# the actual names
actual_names <- colnames(study3_all_logl_fp_wider_model_3)

# new names
new_names <- actual_names %>%
  str_remove('(F|T)_') %>%
  str_to_title()

names(new_names) <- actual_names

#-------------------------------------------------------------------------------
study3_all_logl_fp_wider_model_3 %>%
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
  fmt_number(decimals = 3) -> study3_all_logl_fp_table_model_3

# saving the table
gt::gtsave(data = study3_all_logl_fp_table_model_3, 
           filename = 'study3_all_badness_model_3.docx', path = here::here("Online-Supplement", "All-Correlations-Tables", "study-3"))

# Positive----------------------------------------------------------------------

# Model 1
study3_pos_logl_fp_model_1 %>% # removing the frequency column
  pivot_wider(names_from = c(conv, mlit),
              values_from = fp_logl, names_expand = F, id_cols = c(model, ident, bounds)) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  ) -> study3_pos_logl_fp_wider_model_1

#------------------------------------------------------------------------------- 
# the actual names
actual_names <- colnames(study3_pos_logl_fp_wider_model_1)

# new names
new_names <- actual_names %>%
  str_remove('(F|T)_') %>%
  str_to_title()

names(new_names) <- actual_names

#-------------------------------------------------------------------------------
study3_pos_logl_fp_wider_model_1 %>%
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
  fmt_number(decimals = 3) -> study3_pos_logl_fp_table_model_1

# saving the table
gt::gtsave(data = study3_pos_logl_fp_table_model_1, 
  filename = 'study3_pos_badness_model_1.docx',path = here::here("Online-Supplement","Positive-Correlations-Tables", "study-3"))


# Model 2------------------------------------------------------------------------
study3_pos_logl_fp_model_2 %>% # removing the frequency column
  pivot_wider(names_from = c(conv, mlit),
              values_from = fp_logl, names_expand = F, id_cols = c(model, ident, bounds)) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  ) -> study3_pos_logl_fp_wider_model_2

#------------------------------------------------------------------------------- 
# the actual names
actual_names <- colnames(study3_pos_logl_fp_wider_model_2)

# new names
new_names <- actual_names %>%
  str_remove('(F|T)_') %>%
  str_to_title()

names(new_names) <- actual_names

#-------------------------------------------------------------------------------
study3_pos_logl_fp_wider_model_2 %>%
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
  fmt_number(decimals = 3) -> study3_pos_logl_fp_table_model_2

# saving the table
gt::gtsave(data = study3_pos_logl_fp_table_model_2, 
  filename = 'study3_pos_badness_model_2.docx', path = here::here("Online-Supplement","Positive-Correlations-Tables", "study-3"))

# Model 3------------------------------------------------------------------------
study3_pos_logl_fp_model_3 %>% # removing the frequency column
  pivot_wider(names_from = c(conv, mlit),
              values_from = fp_logl, names_expand = F, id_cols = c(model, ident, bounds)) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  ) -> study3_pos_logl_fp_wider_model_3

#------------------------------------------------------------------------------- 
# the actual names
actual_names <- colnames(study3_pos_logl_fp_wider_model_3)

# new names
new_names <- actual_names %>%
  str_remove('(F|T)_') %>%
  str_to_title()

names(new_names) <- actual_names

#-------------------------------------------------------------------------------
study3_pos_logl_fp_wider_model_3 %>%
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
  fmt_number(decimals = 3) -> study3_pos_logl_fp_table_model_3

# saving the table
gt::gtsave(data = study3_pos_logl_fp_table_model_3, 
filename = 'study3_pos_badness_model_3.docx', path = here::here("Online-Supplement","Positive-Correlations-Tables", "study-3"))

# Part 2------------------------------------------------------------------------

# Model 1
study3_all_logl_fp_model_11_true %>% # removing the frequency column
  pivot_wider(names_from = c(conv, mlit),
              values_from = fp_logl, names_expand = F, id_cols = c(model, ident, bounds)) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  ) -> study3_all_logl_fp_wider_model_11

#------------------------------------------------------------------------------- 
# the actual names
actual_names <- colnames(study3_all_logl_fp_wider_model_11)

# new names
new_names <- actual_names %>%
  str_remove('(F|T)_') %>%
  str_to_title()

names(new_names) <- actual_names

#-------------------------------------------------------------------------------
study3_all_logl_fp_wider_model_11 %>%
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
  fmt_number(decimals = 3) -> study3_all_logl_fp_table_model_11

# saving the table
gt::gtsave(data = study3_all_logl_fp_table_model_11, 
           filename = 'study3_all_true_model_1.docx', path = here::here("Online-Supplement", "All-Correlations-Tables", "study-3"))


# Model 2------------------------------------------------------------------------
study3_all_logl_fp_model_22_true %>% # removing the frequency column
  pivot_wider(names_from = c(conv, mlit),
              values_from = fp_logl, names_expand = F, id_cols = c(model, ident, bounds)) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  ) -> study3_all_logl_fp_wider_model_22

#------------------------------------------------------------------------------- 
# the actual names
actual_names <- colnames(study3_all_logl_fp_wider_model_22)

# new names
new_names <- actual_names %>%
  str_remove('(F|T)_') %>%
  str_to_title()

names(new_names) <- actual_names

#-------------------------------------------------------------------------------
study3_all_logl_fp_wider_model_22 %>%
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
  fmt_number(decimals = 3) -> study3_all_logl_fp_table_model_22

# saving the table
gt::gtsave(data = study3_all_logl_fp_table_model_22, 
   filename = 'study3_all_true_model_2.docx', path = here::here("Online-Supplement", "All-Correlations-Tables", "study-3"))


# Model 3------------------------------------------------------------------------
study3_all_logl_fp_model_33_true %>% # removing the frequency column
  pivot_wider(names_from = c(conv, mlit),
              values_from = fp_logl, names_expand = F, id_cols = c(model, ident, bounds)) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  ) -> study3_all_logl_fp_wider_model_33

#------------------------------------------------------------------------------- 
# the actual names
actual_names <- colnames(study3_all_logl_fp_wider_model_33)

# new names
new_names <- actual_names %>%
  str_remove('(F|T)_') %>%
  str_to_title()

names(new_names) <- actual_names

#-------------------------------------------------------------------------------
study3_all_logl_fp_wider_model_33 %>%
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
  fmt_number(decimals = 3) -> study3_all_logl_fp_table_model_33

# saving the table
gt::gtsave(data = study3_all_logl_fp_table_model_33, 
           filename = 'study3_all_true_model_3.docx', path = here::here("Online-Supplement", "All-Correlations-Tables", "study-3"))

# Positive----------------------------------------------------------------------
# Model 1
study3_pos_logl_fp_model_11_true %>% # removing the frequency column
  pivot_wider(names_from = c(conv, mlit),
              values_from = fp_logl, names_expand = F, id_cols = c(model, ident, bounds)) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  ) -> study3_pos_logl_fp_wider_model_11

#------------------------------------------------------------------------------- 
# the actual names
actual_names <- colnames(study3_pos_logl_fp_wider_model_11)

# new names
new_names <- actual_names %>%
  str_remove('(F|T)_') %>%
  str_to_title()

names(new_names) <- actual_names

#-------------------------------------------------------------------------------
study3_pos_logl_fp_wider_model_11 %>%
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
  fmt_number(decimals = 3) -> study3_pos_logl_fp_table_model_11

# saving the table
gt::gtsave(data = study3_pos_logl_fp_table_model_11, 
           filename = 'study3_pos_true_model_1.docx', path = here::here("Online-Supplement","Positive-Correlations-Tables", "study-3"))


# Model 2------------------------------------------------------------------------
study3_pos_logl_fp_model_22_true %>% # removing the frequency column
  pivot_wider(names_from = c(conv, mlit),
              values_from = fp_logl, names_expand = F, id_cols = c(model, ident, bounds)) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  ) -> study3_pos_logl_fp_wider_model_22

#------------------------------------------------------------------------------- 
# the actual names
actual_names <- colnames(study3_pos_logl_fp_wider_model_22)

# new names
new_names <- actual_names %>%
  str_remove('(F|T)_') %>%
  str_to_title()

names(new_names) <- actual_names

#-------------------------------------------------------------------------------
study3_pos_logl_fp_wider_model_22 %>%
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
  fmt_number(decimals = 3) -> study3_pos_logl_fp_table_model_22

# saving the table
gt::gtsave(data = study3_pos_logl_fp_table_model_22, 
           filename = 'study3_pos_true_model_2.docx', path = here::here("Online-Supplement", "Positive-Correlations-Tables", "study-3"))

# Model 3------------------------------------------------------------------------
study3_pos_logl_fp_model_33_true %>% # removing the frequency column
  pivot_wider(names_from = c(conv, mlit),
              values_from = fp_logl, names_expand = F, id_cols = c(model, ident, bounds)) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  ) -> study3_pos_logl_fp_wider_model_33

#------------------------------------------------------------------------------- 
# the actual names
actual_names <- colnames(study3_pos_logl_fp_wider_model_33)

# new names
new_names <- actual_names %>%
  str_remove('(F|T)_') %>%
  str_to_title()

names(new_names) <- actual_names

#-------------------------------------------------------------------------------
study3_pos_logl_fp_wider_model_33 %>%
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
  fmt_number(decimals = 3) -> study3_pos_logl_fp_table_model_33

# saving the table
gt::gtsave(data = study3_pos_logl_fp_table_model_33, 
           filename = 'study3_pos_true_model_3.docx', path = here::here("Online-Supplement", "Positive-Correlations-Tables", "study-3"))
