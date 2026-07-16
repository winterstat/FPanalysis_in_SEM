#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 09/10/2025
#-------------------------------------------------------------------------------
# Understanding the effect of data space size on FP results(average)
#-------------------------------------------------------------------------------

# loading packages
pkgs <- c('tidyverse', 'ggdist', 'ggthemes', 'magrittr', 'dplyr', 'tidyr', 'ggplot2')

lapply(pkgs, library, character.only = T)

## Load full data space(s) ----
study1_pos <- readRDS(here::here("Online-Supplement","Data-Files","study1_poscor_1014.RDS"))

study1_all <- readRDS(here::here("Online-Supplement","Data-Files","study1_allcor_1014.RDS"))

## Extract fit information (list with 1 element for each model)
study1_pos_fits <- study1_pos$fit_list

study1_all_fits <- study1_all$fit_list

## Select range of data space sizes ----
sizes <- c(50, 100, 500, 1000, 5000, 10000, 50000, 100000, 150000, 175000, 200000, 500000, 750000)

## Set up functions to compute average fit index values ----
samfits <- function(fits_all, size = 100) {
  
  rows <- sample(nrow(fits_all), size, replace = FALSE)
  fits = fits_all[rows,]
  
  cfi<- mean(fits[,"cfi"], na.rm = T)
  
  rmsea<- mean(fits[,"rmsea"], na.rm = T)
  
  srmr <- mean(fits[,"srmr"], na.rm = T)
  
  ll <- mean(fits[,"logl"], na.rm = T)
  
  #data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
  
  res<- c(cfi,rmsea,srmr, ll)
  
  names(res)<- c("cfi","rmsea","srmr", "logl")
  res
}

## Get average for different space sizes ----

# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
set.seed(2306)

for(i in 1:100) {
  # Pos cor model 1
  tmp1 <- t(sapply(sizes, function(size) samfits(study1_pos_fits[[1]], size = size)) )
  tmp1 <- data.frame("iter" = i, "model" = 1, "cor" = 'pos', "size" = sizes, tmp1)
  
  # Pos cor model 2
  tmp2 <- data.frame(t(sapply(sizes, function(size) samfits(study1_pos_fits[[2]], size = size)) ))
  tmp2 <- data.frame("iter" = i, "model" = 2, "cor" = 'pos', "size" = sizes, tmp2)
  
  # All cor model 1
  tmp3 <- data.frame(t(sapply(sizes, function(size) samfits(study1_all_fits[[1]], size = size)) ))
  tmp3 <- data.frame("iter" = i, "model" = 1, "cor" = 'all', "size" = sizes, tmp3)
  
  # All cor model 2
  tmp4 <- data.frame(t(sapply(sizes, function(size) samfits(study1_all_fits[[2]], size = size)) ))
  tmp4 <- data.frame("iter" = i, "model" = 2, "cor" = 'all', "size" = sizes, tmp4)
  
  # Combine all cor and model
  tmp <- rbind(tmp1, tmp2, tmp3, tmp4)
  
  if(i == 1) {
    avg.size.results <- tmp
  } else {
    avg.size.results <- rbind(avg.size.results, tmp)
  }
}

saveRDS(avg.size.results, file = here::here("Online-Supplement","Data-Files" , "average_size_results_study_1.RDS"))

# loading the data set
avg.size.results <- readRDS(here::here("Online-Supplement","Data-Files" , "average_size_results_study_1.RDS"))

## Display FP results in a table ----
avg.size.results %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  group_by(model, cor, size, index) %>%
  summarize(m = mean(value),
            se = sd(value))%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>%
  pivot_wider(names_from = c("cor", "index"), values_from = result) %>%
  print(n = 26) -> table_res



# theme-------------------------------------------------------------------------
my_theme <- theme_grey() +
  theme(panel.background = element_blank(),
        legend.position = 'bottom',
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), strip.text = element_text(face = 'bold'),
        strip.background = element_rect(colour = 'grey79'), axis.line = element_line(colour = 'grey79', linewidth = .35))

theme_set(my_theme)

## Display FP results as a figure (limiting to max spacesize of 100,000) ----
plot_1 <- avg.size.results %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  filter(size < 100001) %>% # selecting fit values less than the value specified
  group_by(model, cor, size, index) %>%
  mutate(
    model = factor(model, levels = c(1, 2),
                   labels = c("Model 1A", "Model 1B")),
    index = factor(index,
                   levels = c("cfi", "rmsea", "srmr", "logl"),
                   labels = c("CFI", "RMSEA", "SRMR", "LL")),
    cor = factor(cor,
                 levels = c("pos", "all"),
                 labels = c("Positive Subspace", "Complete Dataspace"))
  ) %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) + stat_lineribbon(
  aes(fill_ramp = after_stat(level)), .width = 0.95, linewidth = 0.4,
  alpha = 0.9) + facet_grid(rows = vars(index), cols = vars(cor), scales = "free_y",
  switch = "y") + labs(x = "Dataspace size", y = "Fit index",
  fill = "Model", linetype = "Model", fill_ramp = 'Confidence Level') + scale_fill_manual(values = c("Model 1A" = "#276419",
  "Model 1B" = "#8e0152")) + scale_linetype_manual(values = c("Model 1A" = "solid","Model 1B" = "dashed")) + scale_x_continuous(breaks = seq(0, 1e6, 1e4),
  labels = scales::label_number(scale_cut = scales::cut_short_scale())) + scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) + geom_vline(xintercept = 1e4, color = "grey29")


# saving the graph
ggsave(plot = plot_1, path = here::here('Online-Supplement','Plots', 'Study One'),
filename = 'Average-Fit-Indices-dataspace-study-1.png', width = 10, height = 8, units = 'in')


#-------------------------------------------------------------------------------
## Total fitting propensity for each fit measure
avg.size.results %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  filter(size < 200001) %>%
  group_by(model, cor, index) %>%
  summarize(total = mean(value)) %>%  
  mutate(total = round(total,3)) %>%
  pivot_wider(names_from = c("cor", "model", "index"), values_from = total) -> total
#-------------------------------------------------------------------------------



# Commented out 
# ggplot(aes(x = size, y = value, fill = model)) +
# stat_lineribbon(aes(fill_ramp = after_stat(level)),
#                 linewidth = .1, alpha = 1, show.legend = c('fill_ramp' = F, 'fill' = T)) +
# scale_fill_discrete("Model") +
# scale_x_continuous(sec.axis = sec_axis(~ . , name = "Data Space", breaks = NULL, labels = NULL)) +
# scale_y_continuous(sec.axis = sec_axis(~ . , name = "Fit Index", breaks = NULL, labels = NULL)) +
# facet_grid(rows = vars(index),cols =vars(cor), scales = "free") +
# labs(x = "Spacesize", y = "Average Fit Index Value") +
# jtools::theme_apa() +
# theme(panel.grid = element_line(colour = 'grey', linewidth = 1.5),
#   legend.position = "bottom",
#   axis.title = element_text(face = "bold"))



