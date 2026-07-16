# Understanding the effect of data space size on FP results

## Load packages ----
library(ockhamSEM)
library(tidyverse)
library(ggdist)

## Load full data space(s) ----
study1_pos <- readRDS("C:/Users/tivek/OneDrive/Desktop/CoSME/OckhamSEM_simulation/results/study1_poscor.RDS")
study1_all <- readRDS("C:/Users/tivek/OneDrive/Desktop/CoSME/OckhamSEM_simulation/results/study1_allcor.RDS")

## Extract fit information (list with 1 element for each model)
study1_pos_fits <- study1_pos$fit_list
study1_all_fits <- study1_all$fit_list

## Select range of data space sizes ----
sizes <- c(50, 100, 500, 1000, 5000, 10000, 50000, 100000, 150000, 175000, 200000, 500000, 750000)

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

# functions to sample X fit values and compute FP
samfits <- function(fits_all, size = 100) {
  rows <- sample(nrow(fits_all), size, replace = FALSE)
  fits = fits_all[rows,]
  
  cfi<- cfi.fp(fits[,"cfi"])
  rmsea<- rmsea.fp(fits[,"rmsea"])
  srmr <- srmr.fp(fits[,"srmr"])
  
  #data.frame("cfi" = cfi, "rmsea" = rmsea, "srmr" = srmr)
  res<- c(cfi,rmsea,srmr)
  names(res)<- c("cfi","rmsea","srmr")
  res
}

## Get FP for different space sizes ----

# Using 100 random samples of varrying sizes so we can get an estimate of uncertainty
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
    fp.size.results <- tmp
  } else {
    fp.size.results <- rbind(fp.size.results, tmp)
  }
}

saveRDS(fp.size.results, file = "results/study1_dataspacesize.RDS")

## Display FP results in a table ----
fp.size.results %>%
  pivot_longer(cols = cfi:srmr, names_to = "index") %>%
  group_by(model, cor, size, index) %>%
  summarize(m = mean(value),
            se = sd(value))%>%
  mutate(result = paste0(round(m,3), " (", round(se, 3),")"))%>%
  select(-m, -se) %>%
  pivot_wider(names_from = c("cor", "model", "index"), values_from = result)

## Display FP results as a figure (limiting to max spacesize of 100,000) ----
fp.size.results %>%
  pivot_longer(cols = cfi:srmr, names_to = "index") %>%
  filter(size < 100001) %>%
  group_by(model, cor, size, index) %>%
  #median_qi(value, .width = c(.5, .8, .95)) %>%
  mutate(model = factor(model)) %>%
  ggplot(aes(x = size, y = value, fill = model)) +
  stat_lineribbon(aes(fill_ramp = after_stat(level)),
                  linewidth = .1, alpha = 1) +
  facet_grid(rows = vars(index),cols =vars(cor), scales = "free") +
  labs(x = "Spacesize", y = "Fitting Propensity (proportion < or > cutoff)") +
  theme_ggdist()

