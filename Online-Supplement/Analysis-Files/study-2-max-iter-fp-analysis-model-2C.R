#------------------------------------------------
# Author: Tive B.P. Khumalo and Sonja D. Winter
# Date: 07/09/2026
#------------------------------------------------

#-------------------------------------------------------------------
# Fitting propensity: Relationship between maximum iterations and FP
#-------------------------------------------------------------------

# Condition(1): complete data space, 10,000 correlation matrices, forced convergence
# reference-indicator, no variance bounds and max.iter = (150, 250, 750,
# 1500, 2500, 4000, 6000, 10,000)

con_1 <- data.frame('scaling' = 'reference', 'bounds' = 'none')

# Condition(2): complete data space, 10,000 correlation matrices, forced convergence
# reference-indicator, positive variance and max.iter = (150, 250, 750,
# 1500, 2500, 4000, 6000, 10,000)

con_2 <- data.frame('scaling' = 'reference', 'bounds' = 'positive')

# Condition(3): complete data space, 10,000 correlation matrices, forced convergence
# unit variance, no variance bounds and max.iter = (150, 250, 750,
# 1500, 2500, 4000, 6000, 10,000)

con_3 <- data.frame('scaling' = 'unit', 'bounds' = 'none')

# Condition(4): complete data space, 10,000 correlation matrices, forced convergence
# unit variance, positive variance and max.iter = (150, 250, 750,
# 1500, 2500, 4000, 6000, 10,000)

con_4 <- data.frame('scaling' = 'unit', 'bounds' = 'positive')


# maximum number of iterations
max_iter <- c(150, 250, 750, 1500, 2500, 4000, 6000, 10000)

# loading the results file
model_2c_results_1 <- readRDS(here::here("Online-Supplement","Data-Files" ,"model_2C_fp_con_1.rds"))
model_2c_results_2 <- readRDS(here::here("Online-Supplement","Data-Files" ,"model_2C_fp_con_2.rds"))
model_2c_results_3 <- readRDS(here::here("Online-Supplement","Data-Files" ,"model_2C_fp_con_3.rds"))
model_2c_results_4 <- readRDS(here::here("Online-Supplement","Data-Files" ,"model_2C_fp_con_4.rds"))

# cut points
#ll_cufoff_1 <- -32000
ll_cufoff_1 <- -32755

# Fit propensity across the maximum number of iterations
fp_values_1 <- sapply(model_2c_results_1, function(x){
  
  # fit data
  fit <- x$fit_list
  
  # FP values
  res <- sapply(fit, function(x) mean(x[, "logl"] >= ll_cufoff_1)) |> as.vector()
  return(res)
})

fp_values_2 <- sapply(model_2c_results_2, function(x){
  
  # fit data
  fit <- x$fit_list
  
  # FP values
  res <- sapply(fit, function(x) mean(x[, "logl"] >= ll_cufoff_1)) |> as.vector()
  return(res)
})


fp_values_3 <- sapply(model_2c_results_3, function(x){
  
  # fit data
  fit <- x$fit_list
  
  # FP values
  res <- sapply(fit, function(x) mean(x[, "logl"] >= ll_cufoff_1)) |> as.vector()
  return(res)
})

fp_values_4 <- sapply(model_2c_results_4, function(x){
  
  # fit data
  fit <- x$fit_list
  
  # FP values
  res <- sapply(fit, function(x) mean(x[, "logl"] >= ll_cufoff_1)) |> as.vector()
  return(res)
})

# loading tidyverse
library(tidyverse)

# column names
col_names <- c('values', 'scaling', 'bounds')

# data frames
fp_1 <- cbind(fp_values_1, con_1)
colnames(fp_1) <- col_names

fp_2 <- cbind(fp_values_2, con_2)
colnames(fp_2) <- col_names

fp_3 <- cbind(fp_values_3, con_3)
colnames(fp_3) <- col_names

fp_4 <- cbind(fp_values_4, con_4)
colnames(fp_4) <- col_names

# merging the data frames
data_results <- bind_rows(
fp_1, fp_2,
fp_3, fp_4)

# adding a new column
data_results$max_iter <- max_iter

# setting the theme
my_theme <- theme_grey() + theme(panel.background = element_blank(),legend.position = 'bottom',
axis.title = element_text(face = "bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), strip.text = element_text(face = 'bold', size = 8.0),
strip.background = element_rect(colour = 'grey79'), axis.line = element_line(colour = 'grey79', linewidth = .35),
axis.text = element_text(size = 8.0)) 

theme_set(my_theme)


# plotting
data_results |>
  filter(scaling == 'reference', bounds == 'none') |>
  ggplot(aes(x = max_iter, y = values)) +
  geom_point() + geom_line(color = "#276419") + labs(x = "Maximum number of iterations", y = "LL > Cutoff") +
  scale_x_continuous(breaks = seq(150, 1.5e4, 500)) +
  geom_vline(xintercept = 1500, color = "grey67", linetype = 'dashed', linewidth = 1) +
  annotate("text", x = 2500, y = Inf, label = "At max.iter = 1500",vjust = 6.5, color = "grey40") -> plot_1


# saving the plots 
ggsave(plot = plot_1, path = here::here("Online-Supplement","Plots", "Study Two"), 
       filename = 'max-iterations-fp-study-2-condition_1-model-2C.png', width = 10, height = 8,
       units = 'in')

data_results |>
  filter(scaling == 'reference', bounds == 'positive') |>
  ggplot(aes(x = max_iter, y = values)) +
  geom_point() + geom_line(color = "#276419") + labs(x = "Maximum number of iterations", y = "Fitting Propensity") +
  scale_x_continuous(breaks = seq(150, 1.5e4, 500)) +
  geom_vline(xintercept = 1500, color = "grey67", linetype = 'dashed', linewidth = 1) +
  annotate("text", x = 2500, y = Inf, label = "At max.iter = 1500",vjust = 6.5, color = "grey40") -> plot_2


# saving the plots 
ggsave(plot = plot_2, path = here::here("Online-Supplement","Plots", "Study Two"), 
       filename = 'max-iterations-fp-study-2-condition_2-model-2C.png', width = 10, height = 8,
       units = 'in')

data_results |>
  filter(scaling == 'unit', bounds == 'none') |>
  ggplot(aes(x = max_iter, y = values)) +
  geom_point() + geom_line(color = "#276419") + labs(x = "Maximum number of iterations", y = "Fitting Propensity") +
  scale_x_continuous(breaks = seq(150, 1.5e4, 500)) +
  geom_vline(xintercept = 4000, color = "grey67", linetype = 'dashed', linewidth = 1) +
  annotate("text", x = 5000, y = Inf, label = "At max.iter = 4000",vjust = 6.5, color = "grey40") -> plot_3


# saving the plots 
ggsave(plot = plot_3, path = here::here("Online-Supplement","Plots", "Study Two"), 
       filename = 'max-iterations-fp-study-2-condition_3-model-2C.png', width = 10, height = 8,
       units = 'in')

data_results |>
  filter(scaling == 'unit', bounds == 'positive') |>
  ggplot(aes(x = max_iter, y = values)) +
  geom_point() + geom_line(color = "#276419") + labs(x = "Maximum number of iterations", y ="Fitting Propensity") +
  scale_x_continuous(breaks = seq(150, 1.5e4, 500)) +
  geom_vline(xintercept = 6000, color = "grey67", linetype = 'dashed', linewidth = 1) +
  annotate("text", x = 7000, y = Inf, label = "At max.iter = 6000",vjust = 6.5, color = "grey40") -> plot_4


# saving the plots 
ggsave(plot = plot_4, path = here::here("Online-Supplement","Plots", "Study Two"), 
       filename = 'max-iterations-fp-study-2-condition_4-model-2C.png', width = 10, height = 8,
       units = 'in')


