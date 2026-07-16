#-------------------------------------------------------------------------------
# Authors: Sonja D. Winter & Tive B.P. Khumalo
# Date: 09/10/2025
#-------------------------------------------------------------------------------
# Combine result files for Study 2 positive correlations data space

# using the here package


# list all output files
files <- list.files("Data-Files/output_all_2",full.names = T)

# create en empty list with 48 elements
output <- vector(mode = "list", length = 48)

# for each output file, load file (temp), save model number, extract fit values
# once a model has been added to its output element once, we use rbind to add the
# later seed results
for(i in seq_along(files)) {
  temp <- readRDS(files[i])
  
  mod_num <- temp[[1]]
  
  if(is.null(output[[mod_num]])) {
    output[[mod_num]] <- as.data.frame(temp[[3]])
  } else {
    output[[mod_num]] <- rbind(output[[mod_num]], as.data.frame(temp[[3]]))
  }
  print(i)
}

# save as RData file (quite large!)
save(output, file = "Data-Files/output_all_2.RData")

files