#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 09/10/2025
#-------------------------------------------------------------------------------
# This program visualizes the distributions of the log-likelihoods from
# both data spaces.
#------------------------------------------------------------------------------

# loading packages
pkgs <- c('tidyverse', 'ggdist', 'ggthemes')

# check if packages are available
rlang::check_installed(pkgs)

lapply(pkgs, library, character.only = T)

## Load full data space(s) ----------------------------------------------------
study1_pos <- readRDS(here::here("Online-Supplement","Data-Files" ,"study1_poscor.RDS"))

study1_all <- readRDS(here::here("Online-Supplement","Data-Files" ,"study1_allcor.RDS"))

## Extract fit information (list with 1 element for each model)
study1_pos_fits <- study1_pos$fit_list

study1_all_fits <- study1_all$fit_list

# removing the first two objects from my environment since they are large files.
rm(study1_all, study1_pos)

# creating a function that extracts the log-likelihood--------------------------

logl_extract <- function(fit_data, fit_index = c('logl', 'cfi', 'srmr', 'rmsea')){
  
  
  # number of data objects stored in the list
  n <- length(fit_data)
  
  # container
  logl <- list()
  
  for(i in 1:n){
    
    logl[[i]] <- fit_data[[i]][, fit_index]
    
  }
  
  # convert list to data frame
  logl <- do.call(cbind, logl)
  
  return(logl)
  
}

# positive correlations data space
logl_pos <- logl_extract(fit_data = study1_pos_fits, fit_index = 'logl')

# labeling the columns
colnames(logl_pos) <- c('M1_pos', 'M2_pos')

# save as a data frame
logl_pos <- logl_pos %>% as.data.frame()

# all correlations data space
logl_all <- logl_extract(fit_data = study1_all_fits, fit_index = 'logl')

# labeling the columns
colnames(logl_all) <- c('M1_all', 'M2_all')

# save as a data frame
logl_all <- logl_all %>% as.data.frame()


# creating a new columns called space in each data frame
logl_all$space <- 1:1e6

logl_pos$space <- 1:1e6


# merging the two datasets
logl_combined <- full_join(logl_pos, logl_all)


# deleting the space column
logl_combined <- logl_combined %>% select(-space)

# reshaping the data------------------------------------------------------------
logl_long <- logl_combined %>%
  pivot_longer(cols = everything(),
  names_to = c('model', 'correlations'),
  values_to = 'logl', names_sep = '_')

# data wrangling
logl_long <- logl_long %>%
  mutate(
    
    # models
    model = factor(model, levels = c('M1', 'M2'), labels = c('Model 1', 'Model 2')),
    
    # correlations
    correlations = factor(correlations, levels = c('pos', 'all'), labels = c('Positive Correlations', 'All Correlations'))
    
  )


# setting the theme for our figures----------------------------------------------
my_theme <- theme(
  
  axis.text = element_text(size = 7.5),
  
  axis.title = element_text(size = 9.9, face = 'bold', colour = 'grey7'),
  
  legend.position = 'bottom', panel.background = element_blank(),
  
  panel.grid = element_blank(), axis.line = element_line(colour = 'grey0', linewidth = .65),
  
  strip.background = element_rect(colour = 'grey79'), strip.text = element_text(face = 'bold', color = 'grey7')
)

# setting the theme-------------------------------------------------------------
theme_set(my_theme)

# density plots
logl_dens <- logl_long %>%
  ggplot(aes(x = logl, fill = model)) + geom_density(alpha = .25) +
  facet_grid(cols = vars(correlations), scales = 'free') +
  scale_fill_manual(values = c('grey25', 'grey76')) + labs(x = 'Loglikelihood', fill = 'Model') +
  scale_x_continuous(breaks = seq(-5000, 3000, 550)) 
 
# saving the density plots
ggsave(plot = logl_dens, path = here::here('Online-Supplement','Plots' , 'Study One'), 
filename = 'loglikelihood-density-plots.png', width = 10, height = 8, units = 'in')

