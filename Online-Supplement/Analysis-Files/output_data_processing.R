# Authors:  Sonja Winter& Tive Khumalo
# Date: 11.05.2025
# ------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Data pre-processing function
#-------------------------------------------------------------------------------

output <- function(data_files, num_list){
  
  #'@description the function takes list of data files and combines them
  #'@param num_list is the number of list files
  #'@returns the function returns combined files
  
  
  
  # container
  results <- vector(mode = 'list', length = num_list)
  
  # combining the files
  for(i in  seq_along(data_files)){
    
    # reading the files
    temp <- readRDS(data_files[i])
    
    mod_num <- temp[[1]]
    
    if(is.null(results[[mod_num]])){
      
      results[[mod_num]] <- as.data.frame(temp[[3]])
      
    }else{
      
      results[[mod_num]] <- rbind(results[[mod_num]], as.data.frame(temp[[3]]))
      
    }
    
    print(i)
  }
  
  return(results)
  
}




