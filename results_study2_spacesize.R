#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 09/30/2025
#-------------------------------------------------------------------------------

# loading packages
pkgs <- c('tidyverse', 'gt', 'ggdist', 'ggthemes')
lapply(pkgs, library, character.only = T)

#Load full data
load(here::here("Online-Supplement","Data-Files","output_poscor_2_all_conditions.RData"))

# positive correlations
study2_pos_fits <- output_conditions

load(here::here("Online-Supplement","Data-Files", "output_all_2_conditions.RData"))
study2_all_fits <- output_conditions

## Select range of data space sizes ----
sizes <- c(50, 100, 500, 1000, 5000, 10000, 50000, 100000, 150000, 175000, 200000, 500000, 750000)
# 
## Set up functions to compute FP (proportion of fit index values above/below cutoff) ----
# 
#------------------------------------------
# Positive correlations
#------------------------------------------
# 
#  function to get prop CFI >= .45
cfi.fp <- function(fits, cutoff = .45) {

  fits[is.na(fits)] <- 0

  return((sum(fits >= cutoff) / length(fits)))
}

# function to get prop rmsea < .27
rmsea.fp <- function(fits, cutoff = .27) {

  fits[is.na(fits)] <- 1

  return((sum(fits < cutoff) / length(fits)))
}

# function to get prop srmr < .09
srmr.fp <- function(fits, cutoff = .09) {

  fits[is.na(fits)] <- 1

  return((sum(fits < cutoff) / length(fits)))
}

# functions to sample x fit values and compute fp
samfits <- function(fits_all, size = 100) {

rows <- sample(1:nrow(fits_all[[1]]), size, replace = false)

# container
results <- matrix(na, length(fits_all), 4)

# renaming columns
colnames(results) <- c('condition', 'cfi', 'rmsea', 'srmr')

for(m in seq_along(fits_all)){

# sub-setting
fits = fits_all[[m]][rows,]

# cfi
cfi<- cfi.fp(fits[,"cfi"])

# rmsea
rmsea<- rmsea.fp(fits[,"rmsea"])

# srmr
srmr <- srmr.fp(fits[,"srmr"])

# data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
res <- c(m, cfi,rmsea,srmr)

# saving results
results[m,] <- res

  }

# output
return(results)

}

## get fp for different space sizes ----

# using 100 random samples of varying sizes so we can get an estimate of uncertainty
for(i in 1:100) {
tmp <- (lapply(sizes, function(size) samfits(study2_pos_fits, size = size)))

tmp <- (do.call(rbind, tmp))

tmp <- cbind(tmp, rep(sizes, each = 48))

colnames(tmp) <- c('condition', 'cfi', 'rmsea', 'srmr', 'size')

if(i == 1) {
   fp.size.results <- tmp
  } else {
     fp.size.results <- rbind(fp.size.results, tmp)
   }
}


# grid of conditions that we can use for post-processing ----
study2_conditions <- expand.grid(bounds = c("none", "pos.var"),
                                 ident = c("refmarker", "unitvar"),
                                 mlit = c('150', '1500'),
                                 conv = c("f", "t"),
                                 model = 1:3)

# merging the data frame with the conditions
study2_conditions$condition <- 1:48

res_pos <- merge(fp.size.results, study2_conditions)

res_pos$model <- factor(res_pos$model)
# saving the results
saverds(fp.size.results, file = here::here("online-supplement","data-files", "study2_poscor_dataspacesize.rds"))
saverds(res_pos, file = here::here("online-supplement","data-files", "results_poscor_study_2_dataspacesize.rds"))

#------------------------------------------
# All correlations
#------------------------------------------

# function to get prop cfi >= .20
cfi.fp <- function(fits, cutoff = .20) {

  fits[is.na(fits)] <- 0

  return((sum(fits >= cutoff) / length(fits)))
}

# function to get prop rmsea < .27
rmsea.fp <- function(fits, cutoff = .27) {

  fits[is.na(fits)] <- 1

  return((sum(fits < cutoff) / length(fits)))
}

# function to get prop srmr < .17
srmr.fp <- function(fits, cutoff = .17) {

  fits[is.na(fits)] <- 1

  return((sum(fits < cutoff) / length(fits)))
}

#Get FP for different space sizes ----
set.seed(3242)
# Using 100 random samples of varying sizes so we can get an estimate of uncertainty
 for(i in 1:100) {
   tmp1 <- (lapply(sizes, function(size) samfits(study2_all_fits, size = size)))

   tmp1 <- (do.call(rbind, tmp1))

   tmp1 <- cbind(tmp1, rep(sizes, each = 48))

   colnames(tmp1) <- c('condition', 'cfi', 'rmsea', 'srmr', 'size')

   if(i == 1) {
     fp.size.results1 <- tmp1
   } else {
     fp.size.results1 <- rbind(fp.size.results1, tmp1)
   }
 }


# merging the data frame with the conditions
res_all <- merge(fp.size.results1, study2_conditions)
res_all$model <- factor(res_all$model)

#saving the results
saveRDS(fp.size.results1, file = here::here("Online-Supplement","Data-Files", "study2_all_dataspacesize.RDS"))
saveRDS(res_all, file = here::here("Online-Supplement","Data-Files", "results_all_study_2_dataspacesize.RDS"))

# merging the two data sets by a common variable called cor
res_pos <- readRDS(here::here("Online-Supplement","Data-Files", "results_poscor_study_2_dataspacesize.RDS"))
res_all <- readRDS(here::here("Online-Supplement","Data-Files", "results_all_study_2_dataspacesize.RDS"))

# new variable in each data frame
results_full <- bind_rows(
  
  res_pos %>% mutate(cor = 'pos') %>% mutate(cor = factor(cor)),
  
  res_all %>% mutate(cor = 'all') %>% mutate(cor = factor(cor))
)

## Display FP results in a table ----
results_full %>%
  pivot_longer(cols = cfi:srmr, names_to = "index") %>%
  group_by(model, cor, size, index) %>%
  summarize(m = mean(value, na.rm = T),
            se = sd(value, na.rm = T))%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>%
  pivot_wider(names_from = c("cor", "index"), values_from = result) %>%
  print(n = 26)
#--------------------------------------------------------------------------

# for this condition: Positive correlations
# Factor scaling: Reference/Unit, Variance Bounds: None

# plots

# CFI----------------------------------------------------------------------------
plot_1_positive <- results_full %>%
  pivot_longer(cols = cfi:srmr, names_to = "index") %>%
  drop_na() %>%
  filter(size <= 10000) %>%
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c("Model 1B", "Model 2B", "Model 3B")),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr"), labels = c("CFI", "RMSEA", "SRMR")), 
         
         bounds = factor(bounds, levels = c("none", "pos.var"), labels = c("Unrestricted", "Positive")),
         
         ident = factor(ident, levels = c("refmarker", "unitvar"), labels = c("Reference", "Unit Variance")),
         
         mlit = factor(mlit, levels = c("150", "1500")),
         
         conv = factor(conv, levels = c("F", "T"), labels = c("Default", "Forced")),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == "CFI" & bounds == "Unrestricted" & cor == "Positive Subspace") %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  .width = 0.95, linewidth = 0.4, alpha = 0.9) +
  facet_grid(rows = vars(mlit),
             cols = vars(conv, ident),
             scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (CFI > .20)",
       fill = "Model", fill_ramp = "Confidence Level", linetype = "Model") +
  theme_grey() +
  theme(panel.background = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 6.0),
        strip.background = element_rect(colour = 'grey79'), axis.line = element_line(colour = 'grey79', linewidth = .35),
        axis.text = element_text(size = 6.0)) +
  scale_fill_manual(name = "Model",
                    values = c("Model 1B" = "#276419",
                               "Model 2B" = "#8e0152", 
                               "Model 3B" = "#b2182b")) + 
  scale_linetype_manual(
    name = "Model",
    values = c("Model 1B" = "solid", "Model 2B" = "dashed", "Model 3B" = "dotted")
  ) +
  scale_x_continuous(breaks = seq(0, 1e4, 2e3),
                     labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                     sec.axis = sec_axis(~ . , name = "Positive Subspace", breaks = NULL, labels = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                     sec.axis = sec_axis(~ . , name = "Maximum Iterations", breaks = NULL, labels = NULL))

ggsave(plot = plot_1_positive, path = here::here("Online-Supplement","Plots", "Study Two"), 
       filename = 'cfi-positive-subspace-unrestricted-bounds-study-2.png', width = 10, height = 8,
       units = 'in')

plot_1_all <- results_full %>%
  pivot_longer(cols = cfi:srmr, names_to = "index") %>%
  drop_na() %>%
  filter(size <= 10000) %>%
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c("Model 1B", "Model 2B", "Model 3B")),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr"), labels = c("CFI", "RMSEA", "SRMR")), 
         
         bounds = factor(bounds, levels = c("none", "pos.var"), labels = c("Unrestricted", "Positive")),
         
         ident = factor(ident, levels = c("refmarker", "unitvar"), labels = c("Reference", "Unit Variance")),
         
         mlit = factor(mlit, levels = c("150", "1500")),
         
         conv = factor(conv, levels = c("F", "T"), labels = c("Default", "Forced")),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == "CFI" & bounds == "Unrestricted" & cor == "Complete Space") %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  .width = 0.95, linewidth = 0.4, alpha = 0.9) +
  facet_grid(rows = vars(mlit),
             cols = vars(conv, ident),
             scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (CFI > .20)",
       fill = "Model", fill_ramp = "Confidence Level", linetype = "Model") +
  theme_grey() +
  theme(panel.background = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 6.0),
        strip.background = element_rect(colour = 'grey79'), axis.line = element_line(colour = 'grey79', linewidth = .35),
        axis.text = element_text(size = 6.0)) +
  scale_fill_manual(name = "Model",
                    values = c("Model 1B" = "#276419",
                               "Model 2B" = "#8e0152", 
                               "Model 3B" = "#b2182b")) + 
  scale_linetype_manual(
    name = "Model",
    values = c("Model 1B" = "solid", "Model 2B" = "dashed", "Model 3B" = "dotted")
  ) +
  scale_x_continuous(breaks = seq(0, 1e4, 2e3),
                     labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                     sec.axis = sec_axis(~ . , name = "Complete Space", breaks = NULL, labels = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                     sec.axis = sec_axis(~ . , name = "Maximum Iterations", breaks = NULL, labels = NULL))

ggsave(plot = plot_1_all, path = here::here("Online-Supplement","Plots", "Study Two"), 
       filename = 'cfi-complete-space-unrestricted-bounds-study-2.png', width = 10, height = 8,
       units = 'in')

# SRMR----------------------------------------------------------------------------
plot_2_positive <- results_full %>%
  pivot_longer(cols = cfi:srmr, names_to = "index") %>%
  drop_na() %>%
  filter(size <= 10000) %>%
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c("Model 1B", "Model 2B", "Model 3B")),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr"), labels = c("CFI", "RMSEA", "SRMR")), 
         
         bounds = factor(bounds, levels = c("none", "pos.var"), labels = c("Unrestricted", "Positive")),
         
         ident = factor(ident, levels = c("refmarker", "unitvar"), labels = c("Reference", "Unit Variance")),
         
         mlit = factor(mlit, levels = c("150", "1500")),
         
         conv = factor(conv, levels = c("F", "T"), labels = c("Default", "Forced")),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == "SRMR" & bounds == "Unrestricted" & cor == "Positive Subspace") %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  .width = 0.95, linewidth = 0.4, alpha = 0.9) +
  facet_grid(rows = vars(mlit),
             cols = vars(conv, ident),
             scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (SRMR < .17)",
       fill = "Model", fill_ramp = "Confidence Level", linetype = "Model") +
  theme_grey() +
  theme(panel.background = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 6.0),
        strip.background = element_rect(colour = 'grey79'), axis.line = element_line(colour = 'grey79', linewidth = .35),
        axis.text = element_text(size = 6.0)) +
  scale_fill_manual(name = "Model",
                    values = c("Model 1B" = "#276419",
                               "Model 2B" = "#8e0152", 
                               "Model 3B" = "#b2182b")) + 
  scale_linetype_manual(
    name = "Model",
    values = c("Model 1B" = "solid", "Model 2B" = "dashed", "Model 3B" = "dotted")
  ) +
  scale_x_continuous(breaks = seq(0, 1e4, 2e3),
                     labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                     sec.axis = sec_axis(~ . , name = "Positive Subspace", breaks = NULL, labels = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                     sec.axis = sec_axis(~ . , name = "Maximum Iterations", breaks = NULL, labels = NULL))

ggsave(plot = plot_2_positive, path = here::here("Online-Supplement","Plots", "Study Two"), 
       filename = 'srmr-positive-subspace-unrestricted-bounds-study-2.png', width = 10, height = 8,
       units = 'in')

plot_2_all <- results_full %>%
  pivot_longer(cols = cfi:srmr, names_to = "index") %>%
  drop_na() %>%
  filter(size <= 10000) %>%
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c("Model 1B", "Model 2B", "Model 3B")),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr"), labels = c("CFI", "RMSEA", "SRMR")), 
         
         bounds = factor(bounds, levels = c("none", "pos.var"), labels = c("Unrestricted", "Positive")),
         
         ident = factor(ident, levels = c("refmarker", "unitvar"), labels = c("Reference", "Unit Variance")),
         
         mlit = factor(mlit, levels = c("150", "1500")),
         
         conv = factor(conv, levels = c("F", "T"), labels = c("Default", "Forced")),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == "SRMR" & bounds == "Unrestricted" & cor == "Complete Space") %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  .width = 0.95, linewidth = 0.4, alpha = 0.9) +
  facet_grid(rows = vars(mlit),
             cols = vars(conv, ident),
             scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (SRMR < .17)",
       fill = "Model", fill_ramp = "Confidence Level", linetype = "Model") +
  theme_grey() +
  theme(panel.background = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 6.0),
        strip.background = element_rect(colour = 'grey79'), axis.line = element_line(colour = 'grey79', linewidth = .35),
        axis.text = element_text(size = 6.0)) +
  scale_fill_manual(name = "Model",
                    values = c("Model 1B" = "#276419",
                               "Model 2B" = "#8e0152", 
                               "Model 3B" = "#b2182b")) + 
  scale_linetype_manual(
    name = "Model",
    values = c("Model 1B" = "solid", "Model 2B" = "dashed", "Model 3B" = "dotted")
  ) +
  scale_x_continuous(breaks = seq(0, 1e4, 2e3),
                     labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                     sec.axis = sec_axis(~ . , name = "Complete Space", breaks = NULL, labels = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                     sec.axis = sec_axis(~ . , name = "Maximum Iterations", breaks = NULL, labels = NULL))

ggsave(plot = plot_2_all, path = here::here("Online-Supplement","Plots", "Study Two"),
       filename = 'srmr-complete-space-unrestricted-bounds.png-study-2.png', width = 10, height = 8,
       units = 'in')

# RMSEA----------------------------------------------------------------------------
plot_3_positive <- results_full %>%
  pivot_longer(cols = cfi:srmr, names_to = "index") %>%
  drop_na() %>%
  filter(size <= 10000) %>%
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c("Model 1B", "Model 2B", "Model 3B")),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr"), labels = c("CFI", "RMSEA", "SRMR")), 
         
         bounds = factor(bounds, levels = c("none", "pos.var"), labels = c("Unrestricted", "Positive")),
         
         ident = factor(ident, levels = c("refmarker", "unitvar"), labels = c("Reference", "Unit Variance")),
         
         mlit = factor(mlit, levels = c("150", "1500")),
         
         conv = factor(conv, levels = c("F", "T"), labels = c("Default", "Forced")),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == "RMSEA" & bounds == "Unrestricted" & cor == "Positive Subspace") %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  .width = 0.95, linewidth = 0.4, alpha = 0.9) +
  facet_grid(rows = vars(mlit),
             cols = vars(conv, ident),
             scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (RMSEA < .27)",
       fill = "Model", fill_ramp = "Confidence Level", linetype = "Model") +
  theme_grey() +
  theme(panel.background = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 6.0),
        strip.background = element_rect(colour = 'grey79'), axis.line = element_line(colour = 'grey79', linewidth = .35),
        axis.text = element_text(size = 6.0)) +
  scale_fill_manual(name = "Model",
                    values = c("Model 1B" = "#276419",
                               "Model 2B" = "#8e0152", 
                               "Model 3B" = "#b2182b")) + 
  scale_linetype_manual(
    name = "Model",
    values = c("Model 1B" = "solid", "Model 2B" = "dashed", "Model 3B" = "dotted")
  ) +
  scale_x_continuous(breaks = seq(0, 1e4, 2e3),
                     labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                     sec.axis = sec_axis(~ . , name = "Positive Subspace", breaks = NULL, labels = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                     sec.axis = sec_axis(~ . , name = "Maximum Iterations", breaks = NULL, labels = NULL))

ggsave(plot = plot_3_positive, path = here::here("Online-Supplement","Plots","Study Two"),
       filename = 'rmsea-positive-subspace-unrestricted-bounds-study-2.png', width = 10, height = 8,
       units = 'in')

plot_3_all <- results_full %>%
  pivot_longer(cols = cfi:srmr, names_to = "index") %>%
  drop_na() %>%
  filter(size <= 10000) %>%
  group_by(model, cor, size, index, bounds, ident) %>%
  mutate(model = factor(model, levels = c(1, 2, 3),  labels = c("Model 1B", "Model 2B", "Model 3B")),
         
         index = factor(index, levels = c("cfi", "rmsea", "srmr"), labels = c("CFI", "RMSEA", "SRMR")), 
         
         bounds = factor(bounds, levels = c("none", "pos.var"), labels = c("Unrestricted", "Positive")),
         
         ident = factor(ident, levels = c("refmarker", "unitvar"), labels = c("Reference", "Unit Variance")),
         
         mlit = factor(mlit, levels = c("150", "1500")),
         
         conv = factor(conv, levels = c("F", "T"), labels = c("Default", "Forced")),
         
         cor = factor(cor, levels = c("pos", "all"), labels = c("Positive Subspace", "Complete Space"))) %>%
  filter(index == "RMSEA" & bounds == "Unrestricted" & cor == "Complete Space") %>%
  ggplot(aes(x = size, y = value, fill = model, linetype = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  .width = 0.95, linewidth = 0.4, alpha = 0.9) +
  facet_grid(rows = vars(mlit),
             cols = vars(conv, ident),
             scales = "free") +
  labs(x = "Data  spacesize", 
       y = "Fitting Propensity (RMSEA < .27)",
       fill = "Model", fill_ramp = "Confidence Level", linetype = "Model") +
  theme_grey() +
  theme(panel.background = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 6.0),
        strip.background = element_rect(colour = 'grey79'), axis.line = element_line(colour = 'grey79', linewidth = .35),
        axis.text = element_text(size = 6.0)) +
  scale_fill_manual(name = "Model",
                    values = c("Model 1B" = "#276419",
                               "Model 2B" = "#8e0152", 
                               "Model 3B" = "#b2182b")) + 
  scale_linetype_manual(
    name = "Model",
    values = c("Model 1B" = "solid", "Model 2B" = "dashed", "Model 3B" = "dotted")
  ) +
  scale_x_continuous(breaks = seq(0, 1e4, 2e3),
                     labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                     sec.axis = sec_axis(~ . , name = "Complete Space", breaks = NULL, labels = NULL)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                     sec.axis = sec_axis(~ . , name = "Maximum Iterations", breaks = NULL, labels = NULL))

ggsave(plot = plot_3_all, path = here::here("Online-Supplement","Plots","Study Two"),
       filename = 'rmsea-complete-space-unrestricted-bounds-study-2.png', width = 10, height = 8,
       units = 'in')