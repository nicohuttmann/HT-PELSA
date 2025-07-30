#' Import long-format data and (optionally) prepare the processing command 
#'
#' @param file file location (select when empty)
#' @param data.origin by which software the data was searched 
#' @param experiment.type whether the experimental design was group-wise ("grouped") 
#' or dose-response ("dose_response") 
#' @param silent Turn off feedback from the function?
#' @param ... further arguments if data.origin is "other" like <sample_names>
#' 
#' @importFrom magrittr %>%
#' @returns A tibble of the imported data 
#' @export
#'
#' @examples
#' import_data_long("Data/report.tsv", 
#'   data.origin = "diann", 
#'   experiment.type = "grouped")
import_data_long <- function(file, 
                             data.origin = c("diann", "other"), 
                             experiment.type = c("grouped", "dose_response"), 
                             silent = F, 
                             ...) {
  
  # Check file argument 
  if (!hasArg(file)) 
    file <- choose.files(multi = F)
  
  # Check if file exists 
  if (!file.exists(file)) 
    stop(paste0('File "', file, '" could not be found.'))
  
  
  # Read file 
  if (tools::file_ext(file) == "parquet") 
    data_import <- arrow::read_parquet(file)
  else if (tools::file_ext(file) %in% c("tsv", "csv", "txt"))
    data_import <- vroom::vroom(file, col_types = vroom::cols()) 
  else {
    data_import <- vroom::vroom(file, col_types = vroom::cols())
    message("Please check that the data was imported correctly. If not, please import the data manually and proceed with the process_data_... function.")
  }
  
  
  
  # Output processing command (if silent = F)
  if (silent) {
    
    return(data_import)
    
    # Check data and suggest next command 
  } else {
    
    
    # Check data origin
    if (length(data.origin) > 1) {
      message('More than one argument was provided for <data.origin>.')
      data.origin <- 
        select.list(choices = c("diann", "other"), 
                    title = "What is the source of your raw data?")
    } else {
      data.origin <- match.arg(data.origin, c("diann", "other"))
    }
    
    # Check experiment type 
    if (length(experiment.type) > 1) {
      message('More than one argument was provided for <experiment.type>.')
      
      experiment.type <- 
        select.list(choices = c("grouped", "dose_response"), 
                    title = "Are you comparing two groups or did you do a dose-response assay?")
    } else {
      experiment.type <- match.arg(experiment.type, c("grouped", "dose_response"))
    }
    
    additional.args <- list(...)
    
    
    # Identify file names 
    if (data.origin == "diann") {
      
      sample_names <- data_import$Run %>% 
        unique() %>% 
        sort()
      
    } else if ("sample_names" %in% names(additional.args)) {
      
      sample_names <- data_import[[additional.args$sample_names]] %>% 
        unique() %>% 
        sort()
      
    } else {
      
      sample_names <- c("Samples_Control_A", 
                        "Samples_Control_B", 
                        "Samples_Control_C", 
                        "Samples_Treatment_A", 
                        "Samples_Treatment_B", 
                        "Samples_Treatment_C")
      message("No sample names were found. Please specify during further processing.")
      
    }
    
    
    
    # Prepare vector for sample names 
    sample_names <- sample_names %>% 
      setNames(stringr::str_remove(., .longest_common_prefix(.)), .)
    
    samples_n <- paste0(names(sample_names), '" = "', sample_names, '"')
    samples_v <- paste0('c("', paste(samples_n, collapse = ',\n\t"'), ')')
    
    
    # Prepare vector for sample groups  
    sample_groups <- names(sample_names) %>% 
      setNames(rep("Control", length(sample_names)), .)
    
    groups_n <- paste0(names(sample_groups), '" = "', sample_groups, '"')
    groups_v <- paste0('c("', paste(groups_n, collapse = ',\n\t"'), ')')
    
    
    # Prepare vector for sample concentrations  
    sample_concentrations <- names(sample_names) %>% 
      setNames(rep("1e-0", length(sample_names)), .)
    
    concentrations_n <- paste0(names(sample_concentrations), '" = ', sample_concentrations)
    concentrations_v <- paste0('c("', paste(concentrations_n, collapse = ',\n\t"'), ')')
    
    
    # Prepare vector for sample repicates annotation 
    sample_replicates <- names(sample_names) %>% 
      setNames(rep("rep1", length(sample_names)), .)
    
    replicates_n <- paste0(names(sample_replicates), '" = "', sample_replicates, '"')
    replicates_v <- paste0('c("', paste(replicates_n, collapse = ',\n\t"'), ')')
    
    
    
    # Different commands for gr and dr experiments 
    if (experiment.type == "grouped") {
      
      # Suggest process function
      
      to_cat1 <- paste0('# Preprocess grouped data\n', 
                        'data_processed <- process_data_grouped(data_raw,\n', 
                        "\t# You can change the name for each sample below\n", 
                        '\tsample_names = ', samples_v, ',\n', 
                        "\t# You need to change the groups for each sample below\n", 
                        '\tsample_groups = ', groups_v, ',\n', 
                        '\tmin_fraction_per_group = 1,\n')
      
      
    } else if (experiment.type == "dose_response") {
      
      to_cat1 <- paste0('# Preprocess dose-response data\n', 
                        'data_processed <- process_data_dose_response(data_raw,\n', 
                        "\t# You can change the name for each sample below\n", 
                        '\tsample_names = ', samples_v, ',\n', 
                        "\t# You need to change the treatment concentrations for each sample below\n", 
                        '\tsample_concentrations = ', concentrations_v, ',\n', 
                        "\t# You need to change the treatment replicates for each sample below\n", 
                        '\tsample_replicates = ', replicates_v, ',\n')
      
    }
    
    # Different arguments for diann and other result origins 
    if (data.origin == "diann") {
      
      to_cat2 <- paste0('\tpeptide.column = "Stripped.Sequence",\n',
                        '\tsample.column = "Run",\n', 
                        '\tquant.column = "Precursor.Normalised",\n', 
                        '\tfilter_rows = Lib.Q.Value < 0.01 & Lib.PG.Q.Value < 0.01,\n')
      
    } else if (data.origin == "other") {
      
      to_cat2 <- paste0('\tpeptide.column = "Peptides",\n',
                        '\tsample.column = "Samples",\n', 
                        '\tquant.column = "Quant",\n')
      
    }
    
    
    if (experiment.type == "grouped") {
      
      # Suggest process function
      
      to_cat3 <- paste0('\tmin_fraction_per_group = 1,\n', 
                        '\tnorm.method = "none",\n', 
                        '\tsilent = F)\n')
      
      
    } else if (experiment.type == "dose_response") {
      
      to_cat3 <- paste0('\tmin_fraction_per_rep = 1,\n', 
                        '\tmin_reps_per_peptide = 3,\n', 
                        "\tscale_by_0 = T,\n", 
                        '\tsilent = F)\n')
      
    }
    
    
    # Output command
    
    message("Your data was imported successfully.")
    
    Sys.sleep(2.5)
    
    message("\nPlease copy the code below to prepare the precursor data for the analysis.\n")
    
    Sys.sleep(1.5)
    
    
    to_cat <- paste0(to_cat1, to_cat2, to_cat3)
    
    cat_chars <- paste0(unlist(strsplit(to_cat, split = "\n")), "\n")
    
    cat_final <- c(unlist(strsplit(cat_chars[c(1:4)], split = "")), 
                   cat_chars[-c(1:4)])
    
    cat_speed <- rep(0.05, nchar(to_cat))
    cat_speed[1:which(cat_final == "\n")[4]] <- 0.002
    
    
    for (i in seq_along(cat_final)) {
      cat(cat_final[i])
      Sys.sleep(cat_speed[i])
    }
    
    
    message("Please copy the code above and adjust names and sample information to match your experimental design.")
    
  }
  
  return(data_import)
  
}


#' Import wide-format data and (optionally) prepare the processing command 
#'
#' @param file file location (select when empty)
#' @param data.origin by which software the data was searched 
#' @param experiment.type whether the experimental design was group-wise ("grouped") 
#' or dose-response ("dose_response") 
#' @param silent Turn off feedback from the function?
#' @param ... further arguments if data.origin is "other" like <sample_names>
#' 
#' @importFrom magrittr %>%
#' @returns A tibble of the imported data 
#' @export
#'
#' @examples
#' import_data_wide("Data/report.pr_matrix.tsv", 
#'   data.origin = "diann", 
#'   experiment.type = "grouped")
import_data_wide <- function(file, 
                             data.origin = c("diann", "other"), 
                             experiment.type = c("grouped", "dose_response"), 
                             silent = F, 
                             ...) {
  
  # Check file argument 
  if (!hasArg(file)) 
    file <- choose.files()
  
  # Check if file exists 
  if (!file.exists(file)) 
    stop(paste0('File "', file, '" could not be found.'))
  
  
  # Read file 
  if (tools::file_ext(file) == "parquet") 
    data_import <- arrow::read_parquet(file)
  else if (tools::file_ext(file) %in% c("tsv", "csv", "txt"))
    data_import <- vroom::vroom(file, col_types = vroom::cols()) 
  else {
    data_import <- vroom::vroom(file, col_types = vroom::cols())
    message("Please check that the data was imported correctly. If not, please import the data manually and proceed with the process_data_... function.")
  }
  
  
  
  # Output processing command (if silent = F)
  if (silent) {
    
    if (any(stringr::str_detect(names(data_import), "\\\\"))) {
      names(data_import) <- stringr::str_replace_all(names(data_import), "\\\\", "/")
    }
    
    return(data_import)
    
    # Check data and suggest next command 
  } else {
    
    
    # Check sample names 
    if (any(stringr::str_detect(names(data_import), "\\\\"))) {
      names(data_import) <- stringr::str_replace_all(names(data_import), "\\\\", "/")
      message('"\\" detected in the sample names. Please avoid this in the future. Replacing with "/".')
    }
    
    
    # Check data origin
    if (length(data.origin) > 1) {
      message('More than one argument was provided for <data.origin>.')
      data.origin <- 
        select.list(choices = c("diann", "other"), 
                    title = "What is the source of your raw data?")
    } else {
      data.origin <- match.arg(data.origin, c("diann", "other"))
    }
    
    # Check experiment type 
    if (length(experiment.type) > 1) {
      message('More than one argument was provided for <experiment.type>.')
      
      experiment.type <- 
        select.list(choices = c("grouped", "dose_response"), 
                    title = "Are you comparing two groups or did you do a dose-response assay?")
    } else {
      experiment.type <- match.arg(experiment.type, c("grouped", "dose_response"))
    }
    
    
    additional.args <- list(...)
    
    
    
    # Identify file names 
    if (data.origin == "diann") {
      sample_names <- setdiff(names(data_import), c("Protein.Group",
                                                    "Protein.Ids",
                                                    "Protein.Names",
                                                    "Genes",
                                                    "First.Protein.Description",
                                                    "Proteotypic",
                                                    "Stripped.Sequence",
                                                    "Modified.Sequence",
                                                    "Precursor.Charge",
                                                    "Precursor.Id")) %>% 
        sort()
    } else {
      sample_names <- c("Samples_Control_A", 
                        "Samples_Control_B", 
                        "Samples_Control_C", 
                        "Samples_Treatment_A", 
                        "Samples_Treatment_B", 
                        "Samples_Treatment_C")
      message("No sample names were found. Please specify during further processing.")
    }
    
    
    
    # Prepare vector for sample names 
    sample_names <- sample_names %>% 
      setNames(stringr::str_remove(., .longest_common_prefix(.)), .)
    
    samples_n <- paste0(names(sample_names), '" = "', sample_names, '"')
    samples_v <- paste0('c("', paste(samples_n, collapse = ',\n\t"'), ')')
    
    
    # Prepare vector for sample groups  
    sample_groups <- names(sample_names) %>% 
      setNames(rep("Control", length(sample_names)), .)
    
    groups_n <- paste0(names(sample_groups), '" = "', sample_groups, '"')
    groups_v <- paste0('c("', paste(groups_n, collapse = ',\n\t"'), ')')
    
    
    # Prepare vector for sample concentrations  
    sample_concentrations <- names(sample_names) %>% 
      setNames(rep("1e-0", length(sample_names)), .)
    
    concentrations_n <- paste0(names(sample_concentrations), '" = ', sample_concentrations)
    concentrations_v <- paste0('c("', paste(concentrations_n, collapse = ',\n\t"'), ')')
    
    
    # Prepare vector for sample repicates annotation 
    sample_replicates <- names(sample_names) %>% 
      setNames(rep("rep1", length(sample_names)), .)
    
    replicates_n <- paste0(names(sample_replicates), '" = "', sample_replicates, '"')
    replicates_v <- paste0('c("', paste(replicates_n, collapse = ',\n\t"'), ')')
    
    
    
    # Different commands for gr and dr experiments 
    if (experiment.type == "grouped") {
      
      # Suggest process function
      
      to_cat1 <- paste0('# Preprocess grouped data\n', 
                        'data_processed <- \n\tprocess_data_grouped(\n\tdata_raw,\n', 
                        "\t# You can change the name for each sample below\n", 
                        '\tsample_names = ', samples_v, ',\n', 
                        "\t# You need to change the groups for each sample below\n", 
                        '\tsample_groups = ', groups_v, ',\n')
      
      
    } else if (experiment.type == "dose_response") {
      
      to_cat1 <- paste0('# Preprocess dose-response data\n', 
                        'data_processed <- \n\tprocess_data_dose_response(\n\tdata_raw,\n', 
                        "\t# You can change the name for each sample below\n", 
                        '\tsample_names = ', samples_v, ',\n', 
                        "\t# You need to change the treatment concentrations for each sample below\n", 
                        '\tsample_concentrations = ', concentrations_v, ',\n', 
                        "\t# You need to change the treatment replicates for each sample below\n", 
                        '\tsample_replicates = ', replicates_v, ',\n')
      
    }
    
    # Different arguments for diann and other result origins 
    if (data.origin == "diann") {
      
      to_cat2 <- paste0('\tpeptide.column = "Stripped.Sequence",\n')
      
    } else if (data.origin == "other") {
      
      to_cat2 <- paste0('\tpeptide.column = "Peptides",\n')
      
    }
    
    
    if (experiment.type == "grouped") {
      
      # Suggest process function
      
      to_cat3 <- paste0('\tmin_fraction_per_group = 1,\n', 
                        '\tnorm.method = "none",\n', 
                        '\tsilent = F)\n')
      
      
    } else if (experiment.type == "dose_response") {
      
      to_cat3 <- paste0('\tmin_fraction_per_rep = 1,\n', 
                        '\tmin_reps_per_peptide = 3,\n', 
                        "\tscale_by_0 = T,\n", 
                        '\tsilent = F)\n')
      
    }
    
    
    # Output command
    message("Your data was imported successfully.")
    
    Sys.sleep(2.5)
    
    message("\nPlease copy the code below to prepare the precursor data for the analysis.\n")
    
    Sys.sleep(1.5)
    
    
    to_cat <- paste0(to_cat1, to_cat2, to_cat3)
    
    cat_chars <- paste0(unlist(strsplit(to_cat, split = "\n")), "\n")
    
    cat_final <- c(unlist(strsplit(cat_chars[c(1:4)], split = "")), 
                   cat_chars[-c(1:4)])
    
    cat_speed <- rep(0.05, nchar(to_cat))
    cat_speed[1:which(cat_final == "\n")[4]] <- 0.002
    
    
    for (i in seq_along(cat_final)) {
      cat(cat_final[i])
      Sys.sleep(cat_speed[i])
    }
    
    
    message("Please copy the code above and adjust names and sample information to match your experimental design.")
    
    return(data_import)
  }
  
}

