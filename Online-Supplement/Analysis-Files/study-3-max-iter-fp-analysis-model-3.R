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


# maximum number of iterations
max_iter <- c(150, 250, 750, 1500, 2500, 4000, 6000, 10000)

# loading the results file
model_3a_results_1 <- readRDS(here::here("Online-Supplement","Data-Files" ,"model_3A_fp_con_1.rds"))
model_3b_results_1 <- readRDS(here::here("Online-Supplement","Data-Files" ,"model_3B_fp_con_1.rds"))
model_3c_results_1 <- readRDS(here::here("Online-Supplement","Data-Files" ,"model_3C_fp_con_1.rds"))

# cutoff
ll_cufoff_1 <- -11460
ll_cufoff_2 <- -36999
ll_cufoff_3 <- -62493

# Fit propensity across the maximum number of iterations
fp_values_1 <- sapply(model_3a_results_1, function(x){
  
  # fit data
  fit <- x$fit_list
  
  # FP values
  res <- sapply(fit, function(x) mean(x[, "logl"] >= ll_cufoff_1)) |> as.vector()
  return(res)
})

fp_values_2 <- sapply(model_3b_results_1, function(x){
  
  # fit data
  fit <- x$fit_list
  
  # FP values
  res <- sapply(fit, function(x) mean(x[, "logl"] >= ll_cufoff_2)) |> as.vector()
  return(res)
})


fp_values_3 <- sapply(model_3c_results_1, function(x){
  
  # fit data
  fit <- x$fit_list
  
  # FP values
  res <- sapply(fit, function(x) mean(x[, "logl"] >= ll_cufoff_3)) |> as.vector()
  return(res)
})

# models
models <- factor(rep(c("A", "B", "C"), each = 8))

# data frame
results <- data.frame(values = c(fp_values_1, fp_values_2, fp_values_3),
                      models = models, max_iter = max_iter)

# loading tidyverse
library(tidyverse)

# setting the theme
my_theme <- theme_grey() + theme(panel.background = element_blank(),legend.position = 'bottom',
                                 axis.title = element_text(face = "bold"),panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(), strip.text = element_text(face = 'bold', size = 8.0),
                                 strip.background = element_rect(colour = 'grey79'), 
                                 axis.line = element_line(colour = 'grey79', linewidth = .35),
                                 axis.text = element_text(size = 9.0)) 

theme_set(my_theme)

# figure
results |>
  ggplot(aes(x = max_iter, y = values, color = models)) +
  geom_point(size = 2) + geom_line(linewidth = 1.3) +labs(x = "Maximum Number of Iterations", y = "Fitting Propensity",
  color = "Model") + scale_color_manual(name = "Model 3",
  values = c("A" = "#276419","B" = "#8e0152", "C" = "#b2182b"))  +
  scale_x_continuous(breaks = seq(150, 15000, 500)) + facet_grid(rows = vars(models), scale = 'free_y') +
  geom_vline(xintercept = 1500, color = "grey67",linetype = "dashed", linewidth = 1.3) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                     sec.axis = sec_axis(~ . , name = "Model 3", breaks = NULL, labels = NULL)) -> plot_1

# saving the plots 
ggsave(plot = plot_1, path = here::here("Online-Supplement","Plots", "Study Three"), 
       filename = 'max-iterations-fp-study-3-condition_1-model-3.png', width = 10, height = 8,
       units = 'in')
