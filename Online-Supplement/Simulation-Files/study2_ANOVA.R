#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 09/25/2025
#-------------------------------------------------------------------------------
# ANOVA Analysis for Study 2
#-------------------------------------------------------------------------------

# packages
pkgs <- c('ockhamSEM', 'parallel', 'doParallel', 'lavaan', 'ggthemes', 'tidyverse', 'ggdist', 'effectsize')

# check if this packages have been installed
rlang::check_required(pkgs)

# loading them
lapply(pkgs, library, character.only = TRUE)

# Load full data ----
load("results/positive_correlations_res/output_poscor_2_all_conditions.RData")

# positive correlations
study2_pos_fits <- output_conditions

load("results/all_correlations_res/output_all_2.RData")

## All correlations
study2_all_fits <- output

rm(output, output_conditions)

## Combine all output in one dataframe ----
for(i in 1:48) {
  tmp <- study2_all_fits[[i]]
  tmp$condition <- i
  
  if(i == 1) {
    study2_all_total <- tmp
  } else {
    study2_all_total <- bind_rows(study2_all_total, tmp)
  }
}

for(i in 1:48) {
  tmp <- study2_pos_fits[[i]]
  tmp$condition <- i
  
  if(i == 1) {
    study2_pos_total <- tmp
  } else {
    study2_pos_total <- bind_rows(study2_pos_total, tmp)
  }
}

rm(tmp, study2_all_fits, study2_pos_fits)

# Study 2 conditions ----
study2_conditions <- expand.grid(bounds = c("none", "pos.var"),
                                 ident = c("refmarker", "unitvar"),
                                 mlit = c('150', '1500'),
                                 conv = c("F", "T"),
                                 model = 1:3)

## merging the data frame with the conditions
study2_conditions$condition <- 1:48

study2_all_total <- full_join(study2_all_total, study2_conditions, by = "condition")

study2_all_total <- study2_all_total %>% group_by(condition) %>% mutate(id = row_number())



## ANOVA: all correlations ----
aov_srmr <- aov(srmr ~ bounds + ident + mlit + conv + factor(model), data = study2_all_total)
summary(aov_srmr)
effectsize(aov_srmr) # no big effects

aov_cfi <- aov(cfi ~ bounds + ident + mlit + conv + factor(model), data = study2_all_total)
summary(aov_cfi)
effectsize(aov_cfi) # bounds, conv, model

aov_rmsea <- aov(rmsea ~ bounds + ident + mlit + conv + factor(model), data = study2_all_total)
summary(aov_rmsea)
effectsize(aov_rmsea) # bounds, conv, model (much smaller effects than for cfi and logl)

aov_logl <- aov(logl ~ bounds + ident + mlit + conv + factor(model), data = study2_all_total)
summary(aov_logl)
effectsize(aov_logl) # bounds, conv, model

## ANOVA: positive correlations ----
aov_srmr <- aov(srmr ~ bounds + ident + mlit + conv + factor(model), data = study2_pos_total)
summary(aov_srmr)
effectsize(aov_srmr) # no big effects

aov_cfi <- aov(cfi ~ bounds + ident + mlit + conv + factor(model), data = study2_pos_total)
summary(aov_cfi)
effectsize(aov_cfi) # bounds, conv, model

aov_rmsea <- aov(rmsea ~ bounds + ident + mlit + conv + factor(model), data = study2_pos_total)
summary(aov_rmsea)
effectsize(aov_rmsea) # bounds, conv, model (much smaller effects than for cfi and logl)

aov_logl <- aov(logl ~ bounds + ident + mlit + conv + factor(model), data = study2_pos_total)
summary(aov_logl)
effectsize(aov_logl) # bounds, conv, model

## ANOVA: both data spaces combined ----
aov_srmr <- aov(srmr ~ bounds + ident + mlit + conv + factor(model), data = study2_pos_total)
summary(aov_srmr)
effectsize(aov_srmr) # no big effects

aov_cfi <- aov(cfi ~ bounds + ident + mlit + conv + factor(model), data = study2_pos_total)
summary(aov_cfi)
effectsize(aov_cfi) # bounds, conv, model

aov_rmsea <- aov(rmsea ~ bounds + ident + mlit + conv + factor(model), data = study2_pos_total)
summary(aov_rmsea)
effectsize(aov_rmsea) # bounds, conv, model (much smaller effects than for cfi and logl)

aov_logl <- aov(logl ~ bounds + ident + mlit + conv + factor(model), data = study2_pos_total)
summary(aov_logl)
effectsize(aov_logl) # bounds, conv, model
