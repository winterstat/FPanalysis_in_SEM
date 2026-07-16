#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 09/10/2025
#-------------------------------------------------------------------------------
# Understanding the effect of data space size on FP results
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


## Select range of data space sizes --------------------------------------------
sizes <- c(50, 100, 500, 1000, 5000, 10000, 50000, 100000, 150000, 175000, 200000, 500000, 750000)


# selecting the cut-off for the log-likelihood----------------------------------
pos_log <- sapply(study1_pos_fits, function(x) x[,'logl']) %>%
  as.data.frame()

# labeling the columns
colnames(pos_log) <- c('1', '2')

all_log <- sapply(study1_all_fits, function(x) x[,'logl']) %>%
  as.data.frame()

# labeling the columns
colnames(all_log) <- c('1', '2')

# merge by a third variable called space
pos_log$space <- 'pos'; all_log$space <- 'all'

# merging the above
log_res <- bind_rows(pos_log, all_log)

# reshape
log_res_log <- log_res %>%
  pivot_longer(cols = 1:2, values_to = 'logl', names_to = 'model') %>%
  mutate(
    
    # model
    model = factor(model, levels = c('1', '2'), labels = c('Model 1A', 'Model 1B')),
    
    space = factor(space, levels = c('pos', 'all'), labels = c('Positive Subspace', 'Complete Space'))
    
  ) %>%
  rename(cor = space)
  
# looking at the means
log_res_log %>%
  group_by(cor, model) %>%
  summarize(
   # Mean
    avg = mean(logl, na.rm = TRUE),
    
    # standard deviation
    std = sd(logl, na.rm = TRUE), .groups = 'drop'
  )

#  Summary results
# cor                  model      avg   std
# <fct>                <fct>    <dbl> <dbl>
# 1 Positive Correlation Model 1 -3793.  394.
# 2 Positive Correlation Model 2 -3761.  417.
# 3 All Correlation      Model 1 -3755.  500.
# 4 All Correlation      Model 2 -3869.  381.
#   

# Based on the above table we decided to use -3900

log_res_log %>%
  drop_na() %>%
ggplot(aes(x = logl, colour = model)) +  stat_ecdf(geom = 'step') + 
  labs(x = 'Log-likelihood', y = 'Cumulative proportion', color = 'Factor Scaling') +
  facet_grid(rows = vars(cor))

## Set up functions to compute FP (proportion of fit index values above/below cutoff) ----

# function to get prop CFI >= .95
cfi.fp <- function(fits, cutoff = .95) {
  
  fits[is.na(fits)] <- 0
  
  return((sum(fits >= cutoff) / length(fits)))
}

# function to get prop RMSEA < .10
rmsea.fp <- function(fits, cutoff = .10) {
  
  fits[is.na(fits)] <- 1
  
  return((sum(fits < cutoff) / length(fits)))
}

# function to get prop SRMR < .08
srmr.fp <- function(fits, cutoff = .08) {
  
  fits[is.na(fits)] <- 1
  
  return((sum(fits < cutoff) / length(fits)))
}


# function to get prop LOGL
logl.fp <- function(fits, cutoff = -3900) {
  
  fits[is.na(fits)] <- -4000
  
  return((sum(fits >= cutoff) / length(fits)))
}



# functions to sample X fit values and compute FP
samfits <- function(fits_all, size = 100) {
  
  rows <- sample(nrow(fits_all), size, replace = FALSE)
  
  fits = fits_all[rows,]
  
  cfi<- cfi.fp(fits[,"cfi"])
  
  rmsea<- rmsea.fp(fits[,"rmsea"])
  
  srmr <- srmr.fp(fits[,"srmr"])
  
  logl <- logl.fp(fits[,"logl"])
  
  #data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
  res<- c(cfi, rmsea, srmr, logl)
  
  names(res)<- c("cfi", "rmsea", "srmr", "logl")
  
  res
}

## Get FP for different space sizes ----
set.seed(3242)

# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
  # Positive cor model 1
  tmp1 <- t(sapply(sizes, function(size) samfits(study1_pos_fits[[1]], size = size)))
  
  tmp1 <- data.frame("iter" = i, "model" = 1, "cor" = "pos", "size" = sizes, tmp1)
  
  # Positive cor model 2
  tmp2 <- data.frame(t(sapply(sizes, function(size) samfits(study1_pos_fits[[2]], size = size))))
                     
  tmp2 <- data.frame("iter" = i, "model" = 2, "cor" = "pos", "size" = sizes, tmp2)
  
  # All cor model 1
  tmp3 <- data.frame(t(sapply(sizes, function(size) samfits(study1_all_fits[[1]], size = size))))
  
  tmp3 <- data.frame("iter" = i, "model" = 1, "cor" = "all", "size" = sizes, tmp3)
  
  # All cor model 2
  tmp4 <- data.frame(t(sapply(sizes, function(size) samfits(study1_all_fits[[2]], size = size))))
  
  tmp4 <- data.frame("iter" = i, "model" = 2, "cor" = "all", "size" = sizes, tmp4)
  
  # Combine all cor and model
  tmp <- rbind(tmp1, tmp2, tmp3, tmp4)
  
  if(i == 1) {
    
    fp.size.results <- tmp
    
  } else {
    
    fp.size.results <- rbind(fp.size.results, tmp)
    
  }
}

# save data
saveRDS(fp.size.results, file = here::here("Online-Supplement","Data-Files","study1_dataspacesize.RDS"))

# clear environment
rm(list = ls())

# load the data
fp.size.results <- readRDS(here::here("Online-Supplement","Data-Files","study1_dataspacesize.RDS"))

## Display FP results in a table ----
fp.size.results %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  group_by(model, cor, size, index) %>%
  summarize(m = mean(value),
            se = sd(value))%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>%
  pivot_wider(names_from = c("cor", "model", "index"), values_from = result) -> table_res

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
plot_1 <- fp.size.results %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  filter(size < 100001) %>%
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
    fill = "Model", linetype = "Model", fill_ramp = 'Confidence Level') +  
   scale_fill_manual(values = c("Model 1A" = "#276419",
    "Model 1B" = "#8e0152")) + scale_linetype_manual(values = c("Model 1A" = "solid",
  "Model 1B" = "dashed")) + scale_x_continuous(breaks = seq(0, 1e6, 1e4),
    labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) + geom_vline(xintercept = 1e4, color = "grey29")


# saving the graph
ggsave(plot = plot_1, path = here::here('Online-Supplement', 'Plots', 'Study One'),
filename = 'Fit-Indices-dataspace-study-1.png', width = 10, height = 8, units = 'in')

#-------------------------------------------------------------------------------
## Total fitting propensity for each fit measure
fp.size.results %>%
  pivot_longer(cols = cfi:logl, names_to = "index") %>%
  group_by(cor, model, index) %>%
  summarize(total = mean(value, na.rm = T), .groups = 'drop')%>% 
  pivot_wider(names_from = c("cor", "model", "index"), values_from = total) 
#-------------------------------------------------------------------------------






