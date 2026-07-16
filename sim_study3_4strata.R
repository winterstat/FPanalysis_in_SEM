#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 09/30/2025
#-------------------------------------------------------------------------------
# Study 3: Example of potential outcomes in the face of death (non-convergence)
#-------------------------------------------------------------------------------

# Load packages ----
pkgs <- c('tidyverse', 'gt', 'ggdist', 'ggthemes')

lapply(pkgs, library, character.only = T)

load(here::here("Online-Supplement", "Data-Files", "output_all_3_conditions.RData"))

# all correlations--------------------------------------------------------------

study3_all_fits <- output_conditions

# removing these objects from the environment
rm(output_conditions)

# Study 3 conditions ----
study3_conditions <- expand.grid(model = 1:3, 
                                 bounds = c("none", "pos.var"),
                                 ident = c("refmarker", "unitvar"),
                                 mlit = c(150, 1500),
                                 conv = c("F", "T"))

## merging the data frame with the conditions
study3_conditions$condition <- 1:48

# Example 1: Select conditions to compare ----
comparison <- study3_conditions %>% filter(bounds == "none", mlit == 150, model == 3)

## Extract relevant SRMR values-------------------------------------------------
for(i in 1:nrow(comparison)) {
  element <- comparison$condition[i]
  
  tmp <- data.frame(study3_all_fits[[element]][, "srmr"])
  names(tmp) <- paste0("srmr_", element)
  
  if(i == 1) {
    fourcond_srmr <- tmp
  } else {
    fourcond_srmr <- bind_cols(fourcond_srmr, tmp)
  }
}

## Assign stratum, combine SRMR across converged and forced converged iterations

# column names
colnames(fourcond_srmr) <- c("refmarker_F","unitvar_F", "refmarker_T", "unitvar_T")
fourcond_srmr <- fourcond_srmr %>%
  mutate(stratum = case_when(!is.na(refmarker_F) & !is.na(unitvar_F) ~ "Always-taker", 
                             is.na(refmarker_F) & is.na(unitvar_F) ~ "Never-taker",
                             !is.na(refmarker_F) & is.na(unitvar_F) ~ "Defier (Lref)",
                             is.na(refmarker_F) & !is.na(unitvar_F) ~ "Complier (Lunit)",),
         srmr_refmarker = if_else(is.na(refmarker_F), refmarker_T, refmarker_F),
         srmr_unitvar = if_else(is.na(unitvar_F), unitvar_T, unitvar_F)) 

## Compute stratum specific FP results
fourcond_srmr %>% group_by(stratum) %>%
  summarize(
    freq = n(),
    
    m_refmarker = mean(srmr_refmarker),
    
    m_unitvar = mean(srmr_unitvar),
    
    diff_m = m_unitvar - m_refmarker,
    
    fp_refmarker = sum(srmr_refmarker < .15) / n(),
    
    fp_unitvar = sum(srmr_unitvar < .15) / n(),
    
    diff_fp = fp_unitvar - fp_refmarker) -> srmr_all

#-------------------------------------------------------------------------------

## Extract relevant RMSEA values-------------------------------------------------
for(i in 1:nrow(comparison)) {
  element <- comparison$condition[i]
  
  tmp <- data.frame(study3_all_fits[[element]][, "rmsea"])
  names(tmp) <- paste0("rmsea_", element)
  
  if(i == 1) {
    fourcond_rmsea <- tmp
  } else {
    fourcond_rmsea <- bind_cols(fourcond_rmsea, tmp)
  }
}

## Assign stratum, combine RMSEA across converged and forced converged iterations

# column names
colnames(fourcond_rmsea) <- c("refmarker_F","unitvar_F", "refmarker_T", "unitvar_T")
fourcond_rmsea <- fourcond_rmsea %>%
  mutate(stratum = case_when(!is.na(refmarker_F) & !is.na(unitvar_F) ~ "Always-taker", 
                             is.na(refmarker_F) & is.na(unitvar_F) ~ "Never-taker",
                             !is.na(refmarker_F) & is.na(unitvar_F) ~ "Defier (Lref)",
                             is.na(refmarker_F) & !is.na(unitvar_F) ~ "Complier (Lunit)",),
         rmsea_refmarker = if_else(is.na(refmarker_F), refmarker_T, refmarker_F),
         rmsea_unitvar = if_else(is.na(unitvar_F), unitvar_T, unitvar_F)) 

## Compute stratum specific FP results
fourcond_rmsea %>% group_by(stratum) %>%
  summarize(
    freq = n(),
  
    m_refmarker = mean(rmsea_refmarker),
    
    m_unitvar = mean(rmsea_unitvar),
    
    diff_m = m_unitvar - m_refmarker,
    
    fp_refmarker = sum(rmsea_refmarker < .25) / n(),
    
    fp_unitvar = sum(rmsea_unitvar < .25) / n(),
    
    diff_fp = fp_unitvar - fp_refmarker)  -> rmsea_all


#-------------------------------------------------------------------------------

## Extract relevant CFI values-------------------------------------------------
for(i in 1:nrow(comparison)) {
  element <- comparison$condition[i]
  
  tmp <- data.frame(study3_all_fits[[element]][, "cfi"])
  names(tmp) <- paste0("cfi_", element)
  
  if(i == 1) {
    fourcond_cfi <- tmp
  } else {
    fourcond_cfi <- bind_cols(fourcond_cfi, tmp)
  }
}

## Assign stratum, combine CFI across converged and forced converged iterations
colnames(fourcond_cfi) <- c("refmarker_F","unitvar_F", "refmarker_T", "unitvar_T")
fourcond_cfi <-fourcond_cfi %>%
  mutate(stratum = case_when(!is.na(refmarker_F) & !is.na(unitvar_F) ~ "Always-taker", 
                             is.na(refmarker_F) & is.na(unitvar_F) ~ "Never-taker",
                             !is.na(refmarker_F) & is.na(unitvar_F) ~ "Defier (Lref)",
                             is.na(refmarker_F) & !is.na(unitvar_F) ~ "Complier (Lunit)",),
         cfi_refmarker = if_else(is.na(refmarker_F), refmarker_T, refmarker_F),
         
         cfi_unitvar = if_else(is.na(unitvar_F), unitvar_T, unitvar_F)) 
## Compute stratum specific FP results
fourcond_cfi %>% group_by(stratum) %>%
  summarize(
    
    freq = n(),
    
    m_refmarker = mean(cfi_refmarker),
    
    m_unitvar = mean(cfi_unitvar),
    
    diff_m = m_unitvar - m_refmarker,
    
    fp_refmarker = sum(cfi_refmarker < .15) / n(),
    
    fp_unitvar = sum(cfi_unitvar < .15) / n(),
    
    diff_fp = fp_unitvar - fp_refmarker) -> cfi_all
#-------------------------------------------------------------------------------

# loading the dataset
load(here::here("Online-Supplement","Data-Files", "output_poscor_3_all_conditions.RData"))

study3_pos_fits <- output_conditions

# removing these objects from the environment
rm(output_conditions)

# Use the above code for the positive subspace
