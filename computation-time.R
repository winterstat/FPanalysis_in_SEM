# Authors:Sonja Winter (P.h.D), Tive Khumalo (M.A),
# Date: 06.25.2025
# ------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Computational time computations
#-------------------------------------------------------------------------------

# packages
pkgs <- c('tidyverse', 'gt')
lapply(pkgs, library, character.only = T)

#-------------------------------------------------------------------------------
# Study 2:
#-------------------------------------------------------------------------------

# Positive Correlations
# list all output files
files <- list.files(here::here("Online-Supplement", "Data-Files", "output_pos_sim2_full"),full.names = TRUE)

# container to save the computational time in each conditions
time <- vector(mode = 'list', length = 48)

# for each output file, load file (temp), save model number, extract fit values
# once a model has been added to its output element once, we use rbind to add the
# later seed results

# extracting the time
for(i in seq_along(files)){
  
  temp <- readRDS(files[i])
  
  mod_num <- temp[[1]]
  
  if(is.null(time[[mod_num]])) {
    
    time[[mod_num]] <- as.data.frame(temp[[2]])

  }else{
    time[[mod_num]] <- rbind(time[[mod_num]], as.data.frame(temp[[2]]))
  }
  
  print(i)
  
}

for(i in seq_along(time)){
  
  # column names
  colnames(time[[i]]) <- 'time'
}

 
# Grid of conditions that we can use for post-processing ----
study2_conditions <- expand.grid(bounds = c("none", "pos.var"),
                                 ident = c("refmarker", "unitvar"),
                                 mlit = c('150', '1500'),
                                 conv = c("F", "T"),
                                 model = 1:3)
# models as factors
study2_conditions$model <- factor(study2_conditions$model) 

# number of conditions
num_con <- nrow(study2_conditions)

# number of data frames in the output
k <- length(time)

# check if the above is equal
if(k != num_con){
  stop('error k and num_con should be the same')
}

# container to save the results
time_conditions <- list()

for(i in 1:k){
  time_conditions[[i]] <- bind_cols(time[[i]], study2_conditions[i,])
  
}

# data as a data frame
time_data_1 <- do.call(rbind, time_conditions) 

# table 
time_data_1 %>%
  group_by(model, conv, bounds, ident, mlit) %>%
  summarize(
  # total time
    total = sum(time),
    .groups = 'drop'
  )  %>%
  pivot_wider(names_from = c(conv, mlit),
              values_from = total, names_expand = T) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  )

#-------------------------------------------------------------------------------
# All Correlations
#-------------------------------------------------------------------------------
# list all output files
files <- list.files(here::here("Online-Supplement","Data-Files", "output_all_2"),full.names = T)

# container to save the computational time in each conditions
time <- vector(mode = 'list', length = 48)

# for each output file, load file (temp), save model number, extract fit values
# once a model has been added to its output element once, we use rbind to add the
# later seed results

# extracting the time
for(i in seq_along(files)){
  
  temp <- readRDS(files[i])
  
  mod_num <- temp[[1]]
  
  if(is.null(time[[mod_num]])) {
    
    time[[mod_num]] <- as.data.frame(temp[[2]])
    
  }else{
    time[[mod_num]] <- rbind(time[[mod_num]], as.data.frame(temp[[2]]))
  }
  
  print(i)
  
  
}

for(i in seq_along(time)){

# column names
colnames(time[[i]]) <- 'time'
}

# Grid of conditions that we can use for post-processing ----
study2_conditions <- expand.grid(bounds = c("none", "pos.var"),
                                 ident = c("refmarker", "unitvar"),
                                 mlit = c('150', '1500'),
                                 conv = c("F", "T"),
                                 model = 1:3)
# models as factors
study2_conditions$model <- factor(study2_conditions$model) 

# number of conditions
num_con <- nrow(study2_conditions)

# number of data frames in the output
k <- length(time)

# check if the above is equal
if(k != num_con){
  stop('error k and num_con should be the same')
}

# container to save the results
time_conditions <- list()

for(i in 1:k){
  time_conditions[[i]] <- bind_cols(time[[i]], study2_conditions[i,])
  
}

# data as a data frame
time_data_2 <- do.call(rbind, time_conditions)

# table 
time_data_2 %>%
  group_by(model, conv, bounds, ident, mlit) %>%
  summarize(
    # total time
    total = sum(time),
    .groups = 'drop'
  )  %>%
  pivot_wider(names_from = c(conv, mlit),
              values_from = total, names_expand = T) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  )


#-------------------------------------------------------------------------------
# Study 3:
#-------------------------------------------------------------------------------
# Positive Correlations
#--------------------------------------------------------------------------------

# list all output files
files <- files <- list.files(here::here("Online-Supplement","Data-Files",  "output_pos_sim3_full"),full.names = T)
  
# container to save the computational time in each conditions
time <- vector(mode = 'list', length = 48)

# for each output file, load file (temp), save model number, extract fit values
# once a model has been added to its output element once, we use rbind to add the
# later seed results

# extracting the time
for(i in seq_along(files)){
  
  temp <- readRDS(files[i])
  
  mod_num <- temp[[1]]
  
  if(is.null(time[[mod_num]])) {
    
    time[[mod_num]] <- as.data.frame(temp[[2]])
    
  }else{
    time[[mod_num]] <- rbind(time[[mod_num]], as.data.frame(temp[[2]]))
  }
  
  
  print(i)
  
}

for(i in seq_along(time)){
  
  # column names
  colnames(time[[i]]) <- 'time'
}

# Grid of conditions that we can use for post-processing ----
study3_conditions <- expand.grid(model = 1:3, 
                                 bounds = c("none", "pos.var"),
                                 ident = c("refmarker", "unitvar"),
                                 mlit = c(150, 1500),
                                 conv = c("F", "T"))

# models as factors
study3_conditions$model <- factor(study3_conditions$model)

# number of conditions
num_con <- nrow(study3_conditions)

# number of data frames in the output
k <- length(time)

# check if the above is equal
if(k != num_con){
  stop('error k and num_con should be the same')
}

# container to save the results
time_conditions <- list()

for(i in 1:k){
  time_conditions[[i]] <- bind_cols(time[[i]], study3_conditions[i,])
  
}

# data as a data frame
time_data_3 <- do.call(rbind, time_conditions)


# table 
time_data_3 %>%
  group_by(model, conv, bounds, ident, mlit) %>%
  summarize(
    # total time
    total = sum(time),
    .groups = 'drop'
  )  %>%
  pivot_wider(names_from = c(conv, mlit),
              values_from = total, names_expand = T) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  )

#-------------------------------------------------------------------------------
# All Correlations
#-------------------------------------------------------------------------------

# list all output files
files <- list.files(here::here("Online-Supplement","Data-Files", "output_all_2"),full.names = T)

# container to save the computational time in each conditions
time <- vector(mode = 'list', length = 48)

# for each output file, load file (temp), save model number, extract fit values
# once a model has been added to its output element once, we use rbind to add the
# later seed results

# extracting the time
for(i in seq_along(files)){
  
  temp <- readRDS(files[i])
  
  mod_num <- temp[[1]]
  
  if(is.null(time[[mod_num]])) {
    
    time[[mod_num]] <- as.data.frame(temp[[2]])
    
  }else{
    time[[mod_num]] <- rbind(time[[mod_num]], as.data.frame(temp[[2]]))
  }
  
  print(i)
  
}

for(i in seq_along(time)){
  
  # column names
  colnames(time[[i]]) <- 'time'
}

# Grid of conditions that we can use for post-processing ----
study3_conditions <- expand.grid(model = 1:3, 
                                 bounds = c("none", "pos.var"),
                                 ident = c("refmarker", "unitvar"),
                                 mlit = c(150, 1500),
                                 conv = c("F", "T"))

# models as factors
study3_conditions$model <- factor(study3_conditions$model)

# number of conditions
num_con <- nrow(study3_conditions)

# number of data frames in the output
k <- length(time)

# check if the above is equal
if(k != num_con){
  stop('error k and num_con should be the same')
}

# container to save the results
time_conditions <- list()

for(i in 1:k){
  time_conditions[[i]] <- bind_cols(time[[i]], study3_conditions[i,])
  
}


# data as a data frame
time_data_4 <- do.call(rbind, time_conditions)

# table 
time_data_4 %>%
  group_by(model, conv, bounds, ident, mlit) %>%
  summarize(
    # total time
    total = sum(time),
    .groups = 'drop'
  )  %>%
  pivot_wider(names_from = c(conv, mlit),
              values_from = total, names_expand = T) %>%
  relocate(model) %>% # interchanging the first column with the model column
  mutate(
    ident = as.character(ident),
    model = as.character(model), 
    bounds = paste0('Bounds: ', bounds)
  )

