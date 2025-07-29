#' Helper function to get started 
#'
#' @returns simply prints code to copy and use 
#' @export
#'
#' @examples 
#' \dontrun{
#' .hello_I_want_to_analyse_my_PELSA_data()
#' }
.hello_I_want_to_analyse_my_PELSA_data <- function() {
  
  message("\nSure, we'll get this done in no time.")
  
  Sys.sleep(3)
  
  message("\nLet's start with selecting a data file.")
  
  Sys.sleep(3.5)
  
  message("\nWe'll start in")
  
  Sys.sleep(1.5)
  
  message("3")
  
  Sys.sleep(1)
  
  message("2")
  
  Sys.sleep(0.75)
  
  message("1")
  
  Sys.sleep(0.5)
  
  file <- choose.files(multi = F)
  
  if (length(file) == 1) {
    
    message("\nGreat! Now we need some info about your data.\n")
    
    # Data origin 
    data.origin <- menu(c("DIA-NN", "other"), 
                        title = "What's the origin of this file?")
    
    if (data.origin == 0)
      return(invisible(F))
    else 
      data.origin <- c("diann", "other")[data.origin]
    
    
    # Data type 
    data.type <- menu(c("long", "wide"), 
                        title = "Is the data in long or wide format? (wide = samples are columns)")
    
    if (data.type == 0)
      return(invisible(F))
    else 
      data.type <- c("long", "wide")[data.type]
    
    
    # Experiment type
    experiment.type <- menu(c("grouped", "dose_response"), 
                      title = "And finally, are you comparing two sample groups (1) or did you perform a dose-response experiment (2)?")
    
    if (experiment.type == 0)
      return(invisible(F))
    else 
      experiment.type <- c("grouped", "dose_response")[experiment.type]
    
    
    message("\nThanks! Now you can copy the code below and continue with importing your data.\n")
    
    
    Sys.sleep(2)
    to_cat <- paste0('# Import data\n', 
                     'data_raw <- import_data_', data.type, '("', 
                     gsub("\\\\", "/", file), 
                     '",\n',
                     '\tdata.origin = "', data.origin, '",\n', 
                     '\texperiment.type = "', experiment.type, '",\n', 
                     'silent = F)')
    
    for (i in unlist(strsplit(to_cat, split = ""))) {
      cat(i)
      Sys.sleep(0.001)
    }
    # cat(paste0('data_raw <- import_diann_report("', 
    #            gsub("\\\\", "\\\\\\\\", file), 
    #            '")'))
    
    return(invisible(T))
  } else {
    return(invisible(F))
  }
  
}


#' Identify the longest common prefix of a character vector 
#'
#' @param strs character vector 
#'
#' @returns string of longest common prefix ("" if nothing in common)
#'
#' @examples
#' .longest_common_prefix(c("prefix_A", "prefix_B"))
.longest_common_prefix <- function(strs) {
  if (length(strs) == 0) return("")
  
  min_len <- min(nchar(strs))
  if (min_len == 0) return("")
  
  prefix <- character(0)
  
  for (i in 1:min_len) {
    current_char <- substr(strs[1], i, i)
    if (all(substr(strs, i, i) == current_char)) {
      prefix <- c(prefix, current_char)
    } else {
      break
    }
  }
  
  return(paste(prefix, collapse = ""))
}

