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
comparison <- study3_conditions %>% filter(bounds == "none", mlit == 150, model == 1)

## Extract relevant LOGL values-------------------------------------------------
for(i in 1:nrow(comparison)) {
  element <- comparison$condition[i]
  
  tmp <- data.frame(study3_all_fits[[element]][, "logl"])
  names(tmp) <- paste0("logl_", element)
  
  if(i == 1) {
    fourcond_logl <- tmp
  } else {
    fourcond_logl <- bind_cols(fourcond_logl, tmp)
  }
}

# column names
colnames(fourcond_logl) <- c("refmarker_F","unitvar_F", "refmarker_T", "unitvar_T")

## Assign stratum, combine LOGL across converged and forced converged iterations
fourcond_logl <- fourcond_logl %>%
  mutate(stratum = case_when(!is.na(refmarker_F) & !is.na(unitvar_F) ~ "Always-taker", 
                             is.na(refmarker_F) & is.na(unitvar_F) ~ "Never-taker",
                             !is.na(refmarker_F) & is.na(unitvar_F) ~ "Defier (Lref)",
                             is.na(refmarker_F) & !is.na(unitvar_F) ~ "Complier (Lunit)",),
         # 
         logl_refmarker = if_else(is.na(refmarker_F), refmarker_T, refmarker_F),
         
         logl_unitvar = if_else(is.na(unitvar_F), unitvar_T, unitvar_F)) 

## Compute stratum specific FP results
fourcond_logl %>% group_by(stratum) %>%
  summarize(
    freq = n(),
    
    m_refmarker = mean(logl_refmarker),
    
    m_unitvar = mean(logl_unitvar),
    
    diff_m = m_unitvar - m_refmarker,
    
    fp_refmarker = sum(logl_refmarker > -11460) / n(),
    
    fp_unitvar = sum(logl_unitvar > -11460) / n(),
    
    diff_fp = fp_unitvar - fp_refmarker) -> logl_all_3

# A tibble: 4 × 8
# stratum            freq m_refmarker m_unitvar diff_m fp_refmarker fp_unitvar  diff_fp
# <chr>             <int>       <dbl>     <dbl>  <dbl>        <dbl>      <dbl>    <dbl>
# 1 Always-taker      72308     -11395.   -11411.  -15.8        0.520      0.502 -0.0175 
# 2 Complier (Lunit)  47157     -11269.   -11358.  -88.9        0.616      0.543 -0.0725 
# 3 Defier (Lref)    155629     -11318.   -11361.  -42.6        0.601      0.545 -0.0561 
# 4 Never-taker      724906     -11183.   -11164.   18.3        0.652      0.654  0.00250


save(logl_all_3, file = here::here("Online-Supplement","Data-Files","logl_all_3A.RDS"))
#-------------------------------------------------------------------------------


# setting the theme
my_theme <- theme(panel.background = element_blank(),
                  legend.position = 'bottom',
                  
                  axis.title = element_text(face = "bold"),
                  
                  panel.grid.major = element_blank(),
                  
                  panel.grid.minor = element_blank(), strip.text = element_text(face = 'bold', size = 6.0),
                  
                  strip.background = element_rect(colour = 'grey79'), axis.line = element_line(colour = 'grey79', linewidth = .35),
                  
                  axis.text = element_text(size = 9.5)) 

theme_set(my_theme)


# Cumulative plots--------------------------------------------------------------
fourcond_logl %>%
  select(stratum, logl_refmarker, logl_unitvar) %>%
  pivot_longer(cols = logl_refmarker:logl_unitvar,
               values_to = 'logl', names_to = 'identification') %>%
  separate(col = identification, into = c('index', 'scaling'), sep = '_') %>%
  select(-index) %>%
  mutate(
    stratum = factor(stratum),
    
    scaling = factor(scaling)
  ) %>%
  filter(stratum == 'Always-taker') %>%
  ggplot(aes(x = logl, colour = scaling)) +  stat_ecdf(geom = 'step') + 
  labs(x = 'Log-likelihood', y = 'Cumulative proportion', color = 'Factor Scaling') 

#-------------------------------------------------------------------------------
fourcond_logl %>%
  select(stratum, logl_refmarker, logl_unitvar) %>%
  pivot_longer(cols = logl_refmarker:logl_unitvar,
               values_to = 'logl', names_to = 'identification') %>%
  separate(col = identification, into = c('index', 'scaling'), sep = '_') %>%
  select(-index) %>%
  mutate(
    stratum = factor(stratum),
    
    scaling = factor(scaling)
  ) %>%
  filter(stratum == 'Complier (Lunit)') %>%
  ggplot(aes(x = logl, colour = scaling)) +  stat_ecdf(geom = 'step') + 
  labs(x = 'Log-likelihood', y = 'Cumulative proportion', color = 'Factor Scaling') 
#-------------------------------------------------------------------------------
fourcond_logl %>%
  select(stratum, logl_refmarker, logl_unitvar) %>%
  pivot_longer(cols = logl_refmarker:logl_unitvar,
               values_to = 'logl', names_to = 'identification') %>%
  separate(col = identification, into = c('index', 'scaling'), sep = '_') %>%
  select(-index) %>%
  mutate(
    stratum = factor(stratum),
    
    scaling = factor(scaling)
  ) %>%
  filter(stratum == 'Defier (Lref)') %>%
  ggplot(aes(x = logl, colour = scaling)) +  stat_ecdf(geom = 'step') + 
  labs(x = 'Log-likelihood', y = 'Cumulative proportion', color = 'Factor Scaling') 
#-------------------------------------------------------------------------------
fourcond_logl %>%
  select(stratum, logl_refmarker, logl_unitvar) %>%
  pivot_longer(cols = logl_refmarker:logl_unitvar,
               values_to = 'logl', names_to = 'identification') %>%
  separate(col = identification, into = c('index', 'scaling'), sep = '_') %>%
  select(-index) %>%
  mutate(
    stratum = factor(stratum),
    
    scaling = factor(scaling)
  ) %>%
  filter(stratum == 'Never-taker') %>%
  ggplot(aes(x = logl, colour = scaling)) +  stat_ecdf(geom = 'step') + 
  labs(x = 'Log-likelihood', y = 'Cumulative proportion', color = 'Factor Scaling') 
