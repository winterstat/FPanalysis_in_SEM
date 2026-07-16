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
model_2a_results_1 <- readRDS(here::here("Online-Supplement","Data-Files" ,"model_2A_fp_con_1.rds"))
model_2b_results_1 <- readRDS(here::here("Online-Supplement","Data-Files" ,"model_2B_fp_con_1.rds"))
model_2c_results_1 <- readRDS(here::here("Online-Supplement","Data-Files" ,"model_2C_fp_con_1.rds"))

# cutoff
ll_cufoff_1 <- -32755

# Fit propensity across the maximum number of iterations
fp_values_1 <- sapply(model_2a_results_1, function(x){
  
  # fit data
  fit <- x$fit_list
  
  # FP values
  res <- sapply(fit, function(x) mean(x[, "logl"] >= ll_cufoff_1)) |> as.vector()
  return(res)
})

fp_values_2 <- sapply(model_2b_results_1, function(x){
  
  # fit data
  fit <- x$fit_list
  
  # FP values
  res <- sapply(fit, function(x) mean(x[, "logl"] >= ll_cufoff_1)) |> as.vector()
  return(res)
})


fp_values_3 <- sapply(model_2c_results_1, function(x){
  
  # fit data
  fit <- x$fit_list
  
  # FP values
  res <- sapply(fit, function(x) mean(x[, "logl"] >= ll_cufoff_1)) |> as.vector()
  return(res)
})

# models
models <- factor(rep(c("2A", "2B", "2C"), each = 8), levels = c("2C", "2B", "2A"))

# data frame
results <- data.frame(values = c(fp_values_1, fp_values_2, fp_values_3),
          models = models, max_iter = max_iter)

# loading tidyverse
library(tidyverse); library(ggh4x)

# setting the theme
my_theme <- theme_grey() + theme(panel.background = element_blank(),legend.position = c(.95, .05),
                    legend.justification = c("right", "bottom"),
                    axis.title = element_text(face = "bold"),panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), strip.text = element_text(face = 'bold', size = 8.0),
                    strip.background = element_rect(colour = 'grey79'), 
                    axis.line = element_line(colour = 'grey79', linewidth = .35),
                    axis.text = element_text(size = 9.0)) 

theme_set(my_theme)

# figure
library(ggh4x)

results |>
  ggplot(aes(x = max_iter, y = values, color = models, linetype = models)) +
  geom_point(size = 2.5) + geom_line(linewidth = 1.35) +
  labs(x = "Maximum Number of Iterations", y = "Fitting Propensity",
  color = "Model") +
  scale_color_manual(name = "Model ",
  values = c("2A" = "#332288", "2B" = "#88CCEE", "2C" = "#CC6677")) + scale_linetype_manual(name = "Model ",
  values = c("2A" = "solid", "2B" = "dashed", "2C" = "dotted")) + scale_x_continuous(breaks = seq(150, 15000, 500)) +
  facet_grid(rows = vars(models), scales = "free_y") + facetted_pos_scales(y = list(scale_y_continuous(limits = c(.93, 1.00),
  expand = expansion(mult = c(0.02, 0.02))), scale_y_continuous(limits = c(.69, .77), expand = expansion(mult = c(0.02, 0.02))),
  scale_y_continuous(limits = c(.20, .27), expand = expansion(mult = c(0.02, 0.02))))) +
  geom_vline(xintercept = 1500, color = "grey67", linetype = "dashed",linewidth = 1)-> plot_1


# labels on the other side of the y-axis
gridExtra::grid.arrange(
  plot_1,
  right = grid::textGrob(
    "Model",
    rot = -90,
    gp = grid::gpar(fontsize = 12, fontface = "bold")
  )
) 

png(
  filename = here::here(
    "Online-Supplement", "Plots", "Study Two",
    "max-iterations-fp-study-2-condition_1-model-2.png"
  ),
  width = 10,
  height = 8,
  units = "in",
  res = 300
)

grid::grid.draw(plot_1)

dev.off()