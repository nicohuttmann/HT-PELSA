#' Process the imported DIA-NN precursor data 
#'
#'
#' @param data_raw tibble in long or wide format (e.g. from 
#' import_data_long()/_wide())
#' @param sample_names named vector containing sample names as value and 
#' original sample names as names (specifies which samples to use and renames 
#' them)
#' @param sample_concentrations named vector containing sample concentrations as
#'  value and original sample names as names 
#' @param sample_replicates named vector containing replicate annotation as
#'  value and original sample names as names 
#' @param peptide.column column containing peptide IDs to which the data should 
#' be summarized/collapse (by sum())
#' @param sample.column column containing sample names (when working 
#' with long format data)
#' @param quant.column column containing quantitative information (when working 
#' with long format data)
#' @param filter_rows tidy expression to filter the data e.g. filter_rows = 
#' Lib.Q.Value < 0.01
#' @param min_fraction_per_rep minimum fraction of values per replicate to use 
#' (e.g. 0.8 = min. 4 of 5 concentrations with valid value)
#' @param min_reps_per_peptide minimum number of complete replicates (see 
#' min_fraction_per_rep) to keep a peptide (default = 3; when using 4 
#' replicates)
#' @param scale_by_0 Scale data by values at concentration 0?
#' @param silent Turn off feedback from the function?
#'
#' @importFrom magrittr %>%
#' @returns list of processed peptide data ($data_processed) and used raw data 
#' ($data_raw)
#' @export
#'
#' @examples
#' # see import_data_long()/_wide() output 
process_data_dose_response <- function(data_raw, 
                                       sample_names, 
                                       sample_concentrations, 
                                       sample_replicates, 
                                       peptide.column = "Stripped.Sequence", 
                                       sample.column = "Run", 
                                       quant.column = "Precursor.Normalised", 
                                       filter_rows, 
                                       min_fraction_per_rep = 1, 
                                       min_reps_per_peptide = 3, 
                                       scale_by_0 = T, 
                                       silent = F) {
  
  # Check input
  if (!hasArg(data_raw)) 
    stop("Please provide the imported data with <data_raw>.")
  
  # Check names 
  if (!hasArg(sample_names)) 
    stop("Please define the samples to use and the final sample names with <sample_names>.")
  
  # Check concentrations 
  if (!hasArg(sample_concentrations)) 
    stop("Please add the treatment concentrations per sample with <sample_concentrations>.")
  
  # Check replicates 
  if (!hasArg(sample_replicates)) 
    stop("Please annotate sample replicates with <sample_replicates>.")
  
  
  # Check sample.column and pivot longer if wide
  if (!sample.column %in% names(data_raw)) {
    
    if (all(names(sample_names) %in% names(data_raw))) {
      data_raw <- tidyr::pivot_longer(data_raw, 
                                      cols = names(sample_names), 
                                      values_to = quant.column, 
                                      names_to = sample.column)
    } else {
      stop(paste0('The <sample.column> "', 
                  sample.column, 
                  '" was not found and the given sample names were not found in the column names as well. Please check if you are working with a long or wide data format and the function arguments.'))
    }
  }
  
  # Check quant column 
  if (!peptide.column %in% names(data_raw)) {
    stop(paste0('The <quant.column> "', 
                quant.column, 
                '" was not found in the data. Please check your inputs.'))
  }
  
  # Check peptide.column column 
  if (!peptide.column %in% names(data_raw)) {
    stop(paste0('The <peptide.column> "', 
                peptide.column, 
                '" was not found in the data. Please check your inputs.'))
  }
  
  
  
  # Filter data 
  if (hasArg(filter_rows)) {
    data_filtered <- dplyr::filter(data_raw, !!rlang::enquo(filter_rows)) 
  } else {
    data_filtered <- data_raw
  }
  
  
  
  
  # Rename sample names 
  data_filtered_r <- data_filtered %>% 
    dplyr::mutate(Samples = sample_names[!!rlang::sym(sample.column)], 
                  .before = 1) %>% 
    dplyr::arrange(Samples) %>% 
    dplyr::mutate(Peptides = !!rlang::sym(peptide.column), 
                  .after = 1) %>% 
    dplyr::mutate(Quant = !!rlang::sym(quant.column), 
                  .after = 2)
  
  
  # Also rename concentrations and replicates 
  names(sample_concentrations) <- sample_names[names(sample_concentrations)]
  names(sample_replicates) <- sample_names[names(sample_replicates)]
  
  
  
  # Summarise precursors to modified peptides 
  message(paste0('Your peptides will be summarized at the level of "', 
                 peptide.column, 
                 '".'))
  
  data_quant <- data_filtered_r %>%  
    dplyr::summarise(Quant = sum(Quant, na.rm = T), 
                     .by = c("Samples", "Peptides"))
  
  
  # Filter for number of replicates with valid values 
  data_output <- data_quant %>% 
    # Add column for fits 
    dplyr::mutate(do.fit = T) %>% 
    # Add group information 
    dplyr::mutate(Concentration = sample_concentrations[Samples],
                  Replicate = sample_replicates[Samples],
                  .after = "Samples") %>%
    # Remove missing values 
    dplyr::mutate(Quant = tidyr::replace_na(Quant, 0)) %>%  # NA should not exist 
    # dplyr::filter(!is.na(Quant)) %>%  
    # dplyr::filter(Quant > 0) %>% 
    # dplyr::filter(Quant != Inf) %>% 
    # Count values per peptide and replicate group  
    dplyr::mutate(n = sum(Quant > 0), 
                  .by = c("Peptides", "Replicate")) %>% 
    # Calculate fraction of each peptide in each replicate group 
    dplyr::mutate(p = n / max(n), .by = c("Replicate")) %>% 
    # Remove peptides in reps below min_fraction_per_rep threshold 
    dplyr::mutate(do.fit = p >= min_fraction_per_rep) %>% 
    
    # Remove peptides entirely below min number of reps 
    dplyr::mutate(reps.per.peptide = length(unique(subset(Replicate, do.fit))), 
                  .by = "Peptides") %>% 
    dplyr::mutate(do.fit = ifelse(reps.per.peptide >= min_reps_per_peptide, 
                                  do.fit, 
                                  FALSE)) %>% 
    dplyr::select(-c(n, p, reps.per.peptide)) %>% 
    dplyr::arrange(Replicate, Concentration) %>% 
    dplyr::mutate(found_in = paste(Concentration[Quant > 0], collapse = "/"), 
                  .by = c("Peptides", "Replicate"))
  
  
  # Scale to 0 concentration 
  if (scale_by_0) {
    
    data_output <- data_output %>% 
      dplyr::mutate(quant_0 = ifelse(Concentration == 0, Quant, -1)) %>% 
      dplyr::mutate(Quant = Quant / max(quant_0), 
                    .by = c("Peptides", "Replicate")) %>% 
      dplyr::select(-quant_0)
  }
  
  
  # Add peptide positions 
  
  
  # Prepare output list 
  list_output <- list(data_processed = data_output, 
                      #sample_info = sample_info, 
                      data_raw = data_filtered_r)
  
  
  # Propose next command 
  if (silent) {
    return(list_output)
  } else {
    
    # Propose columns to keep 
    potential_columns <- 
      intersect(purrr::map_lgl(list_output$data_raw, purrr::is_character) %>% 
                  which() %>% 
                  names(), 
                names(list_output$data_raw)
                [stringr::str_detect(names(list_output$data_raw), 
                                     "(P|p)rotein|(G|g)ene")])
    
    if (length(potential_columns) == 0) {
      potential_columns <- ""
    } else {
      potential_columns <- 
        paste0("\tkeep.cols = ", 
               paste0("c(\"", paste(potential_columns, 
                                    collapse = "\",\n\t\t\""), "\")"), 
               ",\n")
    }
    
    
    # Messages
    message()
    message(paste0(sum(dplyr::distinct(data_output, 
                                       Peptides, 
                                       Replicate, 
                                       .keep_all = T)$do.fit), 
                   " of ", 
                   length(dplyr::distinct(data_output, 
                                          Peptides, 
                                          Replicate, 
                                          .keep_all = T)$do.fit), 
                   " peptide replicates will be used for fitting dose-response curves."))
    message()
    message("Your data was successfully processed and can be used for dose-response analysis. Please copy the code below.")
    
    message()
    
    to_cat <- paste0('# Optional prefiltering\n', 
                     'data_processed <- prefilter_data_dose_response(\n', 
                     '\tdata_processed = data_processed,\n', 
                     '\tfilter.thresholds = list("1" = c(0.7, 1.3), \n\t"2" = c(1, 1)))\n\n', 
                     # Analysis part 
                     '# Analyze dose-response data\n', 
                     'data_EC50 <- analyze_data_dose_response(data_processed,\n', 
                     '\tc0.as = log10(1e-12),\n', 
                     '\tR2.threshold = 0.9,\n', 
                     '\tmin_rep_hits = 3,\n', 
                     potential_columns, 
                     '\tmulticore = T,\n', 
                     '\tsilent = F)\n')
    
    for (i in unlist(strsplit(to_cat, split = ""))) {
      cat(i)
      Sys.sleep(0.001)
    }
    
    # Return output list 
    return(list_output)
    
  }
}


#' Process the imported DIA-NN precursor data with name matching for inputs 
#'
#'
#' @param data_raw tibble in long or wide format (e.g. from 
#' import_data_long()/_wide())
#' @param sample_names a list with the arguments pattern and (optional) replace
#' @param sample_concentrations a list with the arguments pattern and (optional) replace
#' @param sample_replicates a list with the arguments pattern and (optional) replace
#' @param peptide.column column containing peptide IDs to which the data should 
#' be summarized/collapse (by sum())
#' @param sample.column column containing sample names (when working 
#' with long format data)
#' @param quant.column column containing quantitative information (when working 
#' with long format data)
#' @param filter_rows tidy expression to filter the data e.g. filter_rows = 
#' Lib.Q.Value < 0.01
#' @param min_fraction_per_rep minimum fraction of values per replicate to use 
#' (e.g. 0.8 = min. 4 of 5 concentrations with valid value)
#' @param min_reps_per_peptide minimum number of complete replicates (see 
#' min_fraction_per_rep) to keep a peptide (default = 3; when using 4 
#' replicates)
#' @param scale_by_0 Scale data by values at concentration 0?
#' @param silent Turn off feedback from the function?
#'
#' @importFrom magrittr %>%
#' @returns list of processed peptide data ($data_processed) and used raw data 
#' ($data_raw)
#' @export
#'
#' @examples
#' # see import_data_long()/_wide() output 
process_data_dose_response_m <- function(data_raw, 
                                         sample_names, 
                                         sample_concentrations, 
                                         sample_replicates, 
                                         peptide.column = "Stripped.Sequence", 
                                         sample.column = "Run", 
                                         quant.column = "Precursor.Normalised", 
                                         filter_rows, 
                                         min_fraction_per_rep = 1, 
                                         min_reps_per_peptide = 3, 
                                         scale_by_0 = T, 
                                         silent = F) {
  
  # Check input
  if (!hasArg(data_raw)) 
    stop("Please provide the imported data with <data_raw>.")
  
  # Check names 
  if (!hasArg(sample_names)) 
    stop("Please define the samples to use and the final sample names with <sample_names>.")
  
  # Check concentrations 
  if (!hasArg(sample_concentrations)) 
    stop("Please add the treatment concentrations per sample with <sample_concentrations>.")
  
  # Check replicates 
  if (!hasArg(sample_replicates)) 
    stop("Please annotate sample replicates with <sample_replicates>.")
  
  
  # Check sample.column and pivot longer if wide
  if (!sample.column %in% names(data_raw)) {
    
    if (all(names(sample_names) %in% names(data_raw))) {
      data_raw <- tidyr::pivot_longer(data_raw, 
                                      cols = names(sample_names), 
                                      values_to = quant.column, 
                                      names_to = sample.column)
    } else {
      stop(paste0('The <sample.column> "', 
                  sample.column, 
                  '" was not found and the given sample names were not found in the column names as well. Please check if you are working with a long or wide data format and the function arguments.'))
    }
  }
  
  # Check quant column 
  if (!peptide.column %in% names(data_raw)) {
    stop(paste0('The <quant.column> "', 
                quant.column, 
                '" was not found in the data. Please check your inputs.'))
  }
  
  # Check peptide.column column 
  if (!peptide.column %in% names(data_raw)) {
    stop(paste0('The <peptide.column> "', 
                peptide.column, 
                '" was not found in the data. Please check your inputs.'))
  }
  
  
  
  # Filter data 
  if (hasArg(filter_rows)) {
    data_filtered <- dplyr::filter(data_raw, !!rlang::enquo(filter_rows)) 
  } else {
    data_filtered <- data_raw
  }
  
  
  # Check input 
  if (!is.list(sample_names) || names(sample_names)[1] != "pattern") 
    stop("Sample names must be declared by a list with the arguments <pattern> and (optional) <replace>.")
  
  
  # Rename columns 
  data_filtered_r <- data_filtered %>% 
    dplyr::mutate(Samples = !!rlang::sym(sample.column), 
                  .before = 1) %>% 
    dplyr::arrange(Samples) %>% 
    dplyr::mutate(Peptides = !!rlang::sym(peptide.column), 
                  .after = 1) %>% 
    dplyr::mutate(Quant = !!rlang::sym(quant.column), 
                  .after = 2)
  
  
  # Summarise precursors to modified peptides 
  message(paste0('Your peptides will be summarized at the level of "', 
                 peptide.column, 
                 '".'))
  # sum of Quant 
  data_quant <- data_filtered_r %>%  
    dplyr::summarise(Quant = sum(Quant, na.rm = T), 
                     .by = c("Samples", "Peptides"))
  
  
  
  # Extract concentrations
  if (!is.list(sample_concentrations) || 
      names(sample_concentrations)[1] != "pattern" || 
      !"replace" %in% names(sample_concentrations)) 
    stop("Sample concentrations must be declared by a list with the arguments <pattern> and (optional) <replace>.")
  
  
  data_quant <- data_quant %>% 
    dplyr::mutate(Concentration = stringr::str_extract(Samples, 
                                                       pattern = sample_concentrations$pattern),
                  .after = "Samples") 
  
  # Replace matched concentrations by sample_concentrations$replace
  for (i in seq_along(sample_concentrations$replace))
    data_quant <- data_quant %>% 
    dplyr::mutate(Concentration = str_replace(Concentration, 
                                              names(sample_concentrations$replace)[i], 
                                              sample_concentrations$replace[i]))
  
  
  # Extract replicates 
  if (!is.list(sample_replicates) || 
      names(sample_replicates)[1] != "pattern" || 
      "replace" %in% names(sample_replicates)) 
    stop("Sample replicates must be declared by a list with the arguments <pattern> and (optional) <replace>.")
  
  data_quant <- data_quant %>% 
    dplyr::mutate(Replicate = stringr::str_extract(Samples, 
                                                   pattern = sample_replicates$pattern),
                  .after = "Concentration")  
  
  # Replace matched replicates by sample_replicates$replace
  if ("replace" %in% names(sample_replicates)) 
    data_quant <- data_quant %>% 
    dplyr::mutate(Replicate = sample_replicates$replace[Replicate])
  
  
  # Rename sample names 
  data_quant <- data_quant %>% 
    dplyr::mutate(Samples = 
                    setNames(stringr::str_extract(Samples,
                                         pattern = sample_names$pattern), 
                             Samples))
  
  
  # Replace matched names by sample_names$replace
  if ("replace" %in% names(sample_names)) 
    data_filtered_r <- data_filtered_r %>% 
    dplyr::mutate(Samples = sample_names$replace[Samples]) %>% 
    dplyr::arrange(Samples)
  
  
  
  # Filter for number of replicates with valid values 
  data_output <- data_quant %>% 
    # Make concentration numeric
    dplyr::mutate(Concentration = as.numeric(Concentration)) %>% 
    # Add column for fits 
    dplyr::mutate(do.fit = T) %>% 
    # Remove missing values 
    dplyr::mutate(Quant = tidyr::replace_na(Quant, 0)) %>%  # NA should not exist 
    # Count values per peptide and replicate group  
    dplyr::mutate(n = sum(Quant > 0), 
                  .by = c("Peptides", "Replicate")) %>% 
    # Calculate fraction of each peptide in each replicate group 
    dplyr::mutate(p = n / max(n), .by = c("Replicate")) %>% 
    # Remove peptides in reps below min_fraction_per_rep threshold 
    dplyr::mutate(do.fit = p >= min_fraction_per_rep) %>% 
    
    # Remove peptides entirely below min number of reps 
    dplyr::mutate(reps.per.peptide = length(unique(subset(Replicate, do.fit))), 
                  .by = "Peptides") %>% 
    dplyr::mutate(do.fit = ifelse(reps.per.peptide >= min_reps_per_peptide, 
                                  do.fit, 
                                  FALSE)) %>% 
    dplyr::select(-c(n, p, reps.per.peptide)) %>% 
    dplyr::arrange(Replicate, Concentration) %>% 
    dplyr::mutate(found_in = paste(Concentration[Quant > 0], collapse = "/"), 
                  .by = c("Peptides", "Replicate"))
  
  
  # Scale to 0 concentration 
  if (scale_by_0) {
    
    data_output <- data_output %>% 
      dplyr::mutate(quant_0 = ifelse(Concentration == 0, Quant, -1)) %>% 
      dplyr::mutate(Quant = Quant / max(quant_0), 
                    .by = c("Peptides", "Replicate")) %>% 
      dplyr::select(-quant_0)
  }
  
  
  # Add peptide positions 
  
  
  # Prepare output list 
  list_output <- list(data_processed = data_output, 
                      #sample_info = sample_info, 
                      data_raw = data_filtered_r)
  
  
  # Propose next command 
  if (silent) {
    return(list_output)
  } else {
    
    # Propose columns to keep 
    potential_columns <- 
      intersect(purrr::map_lgl(list_output$data_raw, purrr::is_character) %>% 
                  which() %>% 
                  names(), 
                names(list_output$data_raw)
                [stringr::str_detect(names(list_output$data_raw), 
                                     "(P|p)rotein|(G|g)ene")])
    
    if (length(potential_columns) == 0) {
      potential_columns <- ""
    } else {
      potential_columns <- 
        paste0("\tkeep.cols = ", 
               paste0("c(\"", paste(potential_columns, 
                                    collapse = "\",\n\t\t\""), "\")"), 
               ",\n")
    }
    
    
    # Messages
    message()
    message(paste0(sum(dplyr::distinct(data_output, 
                                       Peptides, 
                                       Replicate, 
                                       .keep_all = T)$do.fit), 
                   " of ", 
                   length(dplyr::distinct(data_output, 
                                          Peptides, 
                                          Replicate, 
                                          .keep_all = T)$do.fit), 
                   " peptide replicates will be used for fitting dose-response curves."))
    message()
    message("Your data was successfully processed and can be used for dose-response analysis. Please copy the code below.")
    
    message()
    
    to_cat <- paste0('# Optional prefiltering\n', 
                     'data_processed <- prefilter_data_dose_response(\n', 
                     '\tdata_processed = data_processed,\n', 
                     '\tfilter.thresholds = list("1" = c(0.7, 1.3), \n\t"2" = c(1, 1)))\n\n', 
                     # Analysis part 
                     '# Analyze dose-response data\n', 
                     'data_EC50 <- analyze_data_dose_response(data_processed,\n', 
                     '\tc0.as = log10(1e-12),\n', 
                     '\tR2.threshold = 0.9,\n', 
                     '\tmin_rep_hits = 3,\n', 
                     potential_columns, 
                     '\tmulticore = T,\n', 
                     '\tsilent = F)\n')
    
    for (i in unlist(strsplit(to_cat, split = ""))) {
      cat(i)
      Sys.sleep(0.001)
    }
    
    # Return output list 
    return(list_output)
    
  }
}


#' Prefilter dose-response data based on min. change at specified concentrations
#'
#' @param data_processed output list or tibble from process_data_dose_response()
#' @param filter.thresholds list of thresholds to filter peptides in the format 
#' list("rank of concentration" = c(ratio threshold stabilized, ratio threshold 
#' destabilized) with "1" = (1.) highest concentration etc. 
#'
#' @importFrom magrittr %>%
#' @returns list or tibble of prefiltered peptide data for 
#' analyze_data_dose_response()
#' @export
#'
#' @examples 
#' \dontrun{
#' data_processed <- prefilter_data_dose_response(
#'   data_processed = data_processed,
#'   filter.thresholds = list("1" = c(0.7, 1.3), 
#'     "2" = c(1, 1)))
#' }
prefilter_data_dose_response <- function(
    data_processed, 
    filter.thresholds = list("1" = c(0.7, 1.3), 
                             "2" = c(1, 1))) {
  
  # Check if list or data frame input 
  if (tibble::is_tibble(data_processed))
    data_quant <- data_processed
  else if (tibble::is_tibble(data_processed[[1]]))
    data_quant <- data_processed[[1]]
  
  
  
  # Do nothing if no thresholds given 
  if (length(filter.thresholds) > 0) {
    
    if ("do.fit" %in% names(data_quant)) {
      data_quant <- data_quant %>% 
        dplyr::rename(do.fit2 = do.fit)
    }
    
    data_filtered <- data_quant %>% 
      # Set do.fit counter to 0
      dplyr::mutate(do.fit = 0, 
                    # Add rank of concentration 
                    r = match(Concentration, 
                              sort(unique(Concentration), 
                                   decreasing = T))) %>% 
      split(.$r)
    
    # Iterate through given thresholds per concentration ranked by r 
    for (i in names(filter.thresholds)) {
      
      data_filtered[[i]] <- data_filtered[[i]] %>% 
        dplyr::mutate(do.fit = dplyr::case_when(
          Quant < filter.thresholds[[i]][1] ~ 1, 
          Quant > filter.thresholds[[i]][2] ~ -1, 
          .default = 0
        ))
      
    }
    
    # Combine concentrations and threshold peptides per replicate 
    data_filtered_c <- data_filtered %>% 
      dplyr::bind_rows() %>% 
      dplyr::mutate(do.fit = abs(sum(do.fit)) == length(filter.thresholds), 
                    .by = c("Peptides", "Replicate")) %>% 
      dplyr::select(-r) %>% 
      dplyr::arrange(Replicate, Concentration)
    
    if ("do.fit2" %in% names(data_filtered_c)) {
      data_filtered_c <- data_filtered_c %>% 
        dplyr::mutate(do.fit = do.fit & do.fit2) %>% 
        dplyr::select(-do.fit2)
    }
  }
  
  
  # Message about number of peptides 
  message(paste0(sum(dplyr::distinct(data_filtered_c, 
                                     Peptides, 
                                     Replicate, 
                                     .keep_all = T)$do.fit), 
                 " of ", 
                 length(dplyr::distinct(data_filtered_c, 
                                        Peptides, 
                                        Replicate, 
                                        .keep_all = T)$do.fit), 
                 " peptide replicates will be used for fitting dose-response curves."))
  
  # Return filtered data frame 
  if (tibble::is_tibble(data_processed))
    list_output <- data_filtered_c
  else if (tibble::is_tibble(data_processed[[1]])) {
    list_output <- data_processed
    list_output[[1]] <- data_filtered_c
  }
  
  return(list_output)
  
}


#' Fitting function wrapper around drc::drm(..., fct = drc::LL.4()) used by 
#' map(data_fit, this_function) in analyze_data_dose_response()
#'
#' @param x data.frame object with columns log10.Concentration and Quant 
#'
#' @returns tibble containing fit parameters 
#' @export
#'
#' 
.fit_LL.4 <- function(x) {
  
  # Fit model
  model <- tryCatch(suppressMessages(drc::drm(Quant ~ log10.Concentration,
                                              data = x,
                                              fct = drc::LL.4(),
                                              logDose = 10)),
                    error = function(cond) NULL)
  model_coeff <- coef(model)
  
  if (is.null(model))
    return(tibble::tibble(EC50 = NA,
                          R2 = NA,
                          lower = NA,
                          upper = NA,
                          direction = NA))
  else
    return(tibble::tibble(EC50 = model_coeff[[4]],
                          R2 = cor(c(x$Quant),
                                   fitted(model),
                                   method = "pearson")^2,
                          lower = model_coeff[[2]],
                          upper = model_coeff[[3]],
                          direction = ifelse(dplyr::last(x$Quant) > 1,
                                             "destabilized",
                                             "stabilized")))
}


#' Fit dose-response curves to individual peptide replicates
#'
#' @param data_processed output list or tibble from process_data_dose_response()
#'  or prefilter_data_dose_response()
#' @param c0.as replacement for log10(0)
#' @param R2.threshold threshold for individual curves to be accepted
#' @param min_rep_hits minimum number of valid fitted curves per peptide 
#' @param keep.cols Which columns should be copied from data_raw to 
#' data_peptides
#' @param multicore should curve fitting be parallelized (T for automatic 
#' detection of number of cores to use; number to specific manually; F to 
#' disable)
#' @param silent Turn off feedback from the function?
#'
#' @importFrom magrittr %>%
#' @returns list of EC50 results at peptide ($data_peptides) and replicates 
#' ($data_replicates) level 
#' @export
#'
#' @examples 
#' \dontrun{
#' data_EC50 <- analyze_data_dose_response(data_processed,
#'   c0.as = log10(1e-12),
#'   R2.threshold = 0.9,
#'   min_rep_hits = 3,
#'   keep.cols = c("Protein.Group",
#'     "Protein.Ids",
#'     "Protein.Names",
#'     "Genes",
#'     "First.Protein.Description"),
#'   multicore = T,
#'   silent = F)
#' }
analyze_data_dose_response <- function(data_processed, 
                                       c0.as = log10(1e-12), 
                                       R2.threshold = 0.9, 
                                       min_rep_hits = 3, 
                                       keep.cols, 
                                       multicore = T, 
                                       silent = F) {
  
  # Check packages 
  if (!suppressPackageStartupMessages(require(drc, quietly = T, warn.conflicts = F))) {
    if(menu(c("Yes, please install the drc package now.", 
              "No, I will install it myself."), 
            title = ("The drc package is not installed yet, but required for the analysis. Would you like to install it?")) == 1) {
      install.packages("drc")
    } else {
      return(invisible(FALSE))
    }
  }
  
  if (multicore != 0 && 
      !suppressPackageStartupMessages(require(furrr, quietly = T, warn.conflicts = F))) {
    if(menu(c("Yes, please install the furrr package now.", 
              "No, I will install it myself."), 
            title = ("The furrr package is not installed yet, but required for the analysis. Would you like to install it?")) == 1) {
      install.packages("furrr")
    } else {
      return(invisible(FALSE))
    }
  }
  
  
  # Check if list or data frame input 
  if (tibble::is_tibble(data_processed))
    data_quant <- data_processed
  else if (tibble::is_tibble(data_processed[[1]]))
    data_quant <- data_processed[[1]]
  else 
    stop("The <data_processed> must be either a tibble or a list with a tibble as first element. Please use process_data_dose_response() or prefilter_data_dose_response() to get the appropriate data type.")
  
  
  # Prepare data for fitting 
  data_results <- data_quant %>% 
    dplyr::arrange(Concentration) %>% 
    dplyr::mutate(log10.Concentration = 
                    ifelse(Concentration == 0, 
                           c0.as, 
                           log10(Concentration))) %>% 
    #dplyr::select(-Concentration) %>% 
    tidyr::nest(data = c(Quant, log10.Concentration), 
                .by = c("Peptides", "Replicate", "do.fit", "found_in")) %>% 
    dplyr::mutate(index = as.character(dplyr::row_number()))
  
  # Extract nested data as list 
  data_fit <- data_results %>% 
    dplyr::filter(do.fit) %>% 
    dplyr::pull(data, index) %>% 
    purrr::map(as.data.frame)
  
  # Prepare data for results 
  data_results <- data_results %>% 
    dplyr::select(-data)
  
  
  
  # Fit model 
  if (isTRUE(multicore) || multicore > 1) {
    
    
    # Set number of cores 
    if (isTRUE(multicore)) {
      multicore <- future::availableCores() - 2
      if (!silent)
        message(paste0("Setting number of cores to ", multicore, "."))
    }
    
    # Check number of cores again 
    if (multicore < 2 || multicore > future::availableCores()) 
      stop(paste0("Number of cores set: ", 
                  multicore, 
                  ". Choose a number between 2 and your available cores or ", 
                  "disable parallel processing with multicore = F."))
    
    # Set up cores for parallel processing
    future::plan(future::multisession, workers = multicore)
    
    # Start fitting in parallel 
    message(paste0("Begining to fit ", length(data_fit), " peptide curves.",
                   " (The progress bar may take some time to appear.)"))
    
    time <- Sys.time()
    
    data_model <- furrr::future_map(
      data_fit, .fit_LL.4, 
      .options = furrr::furrr_options(
        globals = list(.fit_LL.4 = .fit_LL.4, 
                       data_fit = data_fit)), 
      .progress = T)
    
    time <- difftime(Sys.time(), time, units = "min")
    
    #future:::ClusterRegistry("stop")
    future::plan(future::sequential)
    
    # Not parallel 
  } else {
    
    # Start fitting one one core 
    message(paste0("Begining to fit ", length(data_fit), " peptide curves."))
    
    time <- Sys.time()
    
    data_model <- purrr::map(data_fit, .fit_LL.4, .progress = T)
    
    time <- difftime(Sys.time(), time, units = "min")
    
  }
  
  data_model_c <- dplyr::bind_rows(data_model, .id = "index")
  
  
  # Message when fitting is done 
  message()
  if (!silent)
    message(paste0("Yay, the fitting is done. ", 
                   sum(!is.na(data_model_c$EC50)), "/", 
                   nrow(data_model_c), 
                   " peptide curves were fitted successfully in ", 
                   round(time, 1), " min. Starting post-processing."))
  
  
  # Post-processing 
  data_results_c <- data_results %>% 
    # Add fitting results to results data frame 
    dplyr::left_join(data_model_c, 
                     by = "index") %>% 
    dplyr::select(-index) %>% 
    # Add pEC50
    dplyr::mutate(pEC50 = -log10(EC50), 
                  .after = "EC50") %>% 
    # Evaluate fitting results 
    dplyr::mutate(n = sum(do.fit), 
                  n_fit = sum(!is.na(R2)), 
                  #n_hit = sum(replace_na(R2, 0) > R2.threshold), 
                  n_hit_stabilized = sum(tidyr::replace_na(R2, 0) > R2.threshold & 
                                           direction == "stabilized"), 
                  n_hit_destabilized = sum(tidyr::replace_na(R2, 0) > R2.threshold & 
                                             direction == "destabilized"), 
                  n_hits = max(n_hit_stabilized, n_hit_destabilized), 
                  # Assign regulation type using min_rep_hits
                  regulation_replicate = 
                    map2_chr(R2, direction, 
                             \(R2, direction) 
                             dplyr::case_when(
                               tidyr::replace_na(R2, 0) > R2.threshold & 
                                     direction == "stabilized" ~ "stabilized", 
                               tidyr::replace_na(R2, 0) > R2.threshold & 
                                 direction == "destabilized" ~ "destabilized", 
                               .default = "unresponsive")), 
                  regulation = dplyr::case_when(
                    n_hit_stabilized >= min_rep_hits ~ "stabilized", 
                    n_hit_destabilized >= min_rep_hits ~ "destabilized", 
                    .default = "unresponsive"), 
                  # Remove pEC50 values of replicates below threshold and of 
                  # inverse trend to regulation 
                  pEC50 = ifelse(regulation == direction & 
                                   tidyr::replace_na(R2, 0) > R2.threshold, 
                                 pEC50, 
                                 NA), 
                  EC50 = ifelse(regulation == direction & 
                                  tidyr::replace_na(R2, 0) > R2.threshold, 
                                EC50, 
                                NA), 
                  pEC50_mean = mean(pEC50, na.rm = T), 
                  pEC50_CV = sd(pEC50, na.rm = T) / mean(pEC50, na.rm = T), 
                  .by = "Peptides") %>% 
    # Add max hits 
    # dplyr::mutate(n_hits = max(n_hit_stabilized, n_hit_destabilized), 
    #               .by = c("Peptides")
    # )) %>% 
    # Set unresponsive to NA
    dplyr::mutate(pEC50_mean = ifelse(regulation == "unresponsive", 
                                      NA, 
                                      pEC50_mean), 
                  pEC50_CV = ifelse(regulation == "unresponsive", 
                                    NA, 
                                    pEC50_CV)) %>% 
    # Order by EC50 value 
    dplyr::arrange(dplyr::desc(pEC50_mean))
  
  
  # Extract peptide-centric results
  data_results_p <- data_results_c %>% 
    dplyr::select(Peptides, n_hits, pEC50_mean, pEC50_CV, regulation) %>% 
    dplyr::distinct(Peptides, .keep_all = T)
  
  
  # Transfer annotation columns 
  if (hasArg(keep.cols) && 
      "data_raw" %in% names(data_processed) && 
      any(keep.cols %in% names(data_processed$data_raw))) {
    
    data_results_p <- data_processed$data_raw %>% 
      dplyr::distinct(Peptides, .keep_all = T) %>% 
      dplyr::select(any_of(c("Peptides", keep.cols))) %>% 
      dplyr::right_join(data_results_p, 
                        by = "Peptides")
    
    data_results_c <- data_processed$data_raw %>% 
      dplyr::distinct(Peptides, .keep_all = T) %>% 
      dplyr::select(any_of(c("Peptides", keep.cols))) %>% 
      dplyr::right_join(data_results_c, 
                        by = "Peptides")
    
  }
  
  
  # Success message
  if (!silent) 
    message(paste0("pEC50 values were estimated for ", 
                   sum(data_results_p$regulation == "stabilized"), 
                   " stabilized and ", 
                   sum(data_results_p$regulation == "destabilized"), 
                   " destabilized peptides."))
  
  
  # Plot results 
  p <- data_results_p %>% 
    ggplot2::ggplot(ggplot2::aes(x = pEC50_mean, pEC50_CV, color = regulation)) + 
    ggplot2::geom_point(shape = 16, alpha = 0.5) + 
    ggplot2::theme_classic() + 
    ggplot2::scale_color_manual(values = c(stabilized = "blue", 
                                           destabilized = "red")) + 
    ggplot2::labs(x = expression(pEC[50]), 
                  y = expression(paste("CV ", pEC[50])))
  
  if (!silent) {
    print(p)
  }
  
  # Combine output 
  list_results <- list(data_peptides = data_results_p, 
                       data_replicates = data_results_c)
  
  return(list_results)
  
}


#' Analyze individual concentrations of dose-response experiment as groups 
#'
#' @param data_processed output list or tibble from process_data_dose_response()
#'  or prefilter_data_dose_response()
#' @param filter_rows tidy expression to filter the data e.g. filter_rows = 
#' Lib.Q.Value < 0.01
#' @param min_fraction_per_group minimum fraction of values per replicate to use 
#' (e.g. 0.75 = min. 3 of 4 replicates with valid value)
#' @param norm.method which method to use for normalization 
#' @param silent Turn off feedback from the function?
#' @param p.adjust.method method for p.adjust() (see p.adjust.methods)
#' @param p.threshold p-value threshold after adjustment 
#' @param fc.threshold log2. fold-change threshold 
#' @param eBayes.trend abundance-dependent variance estimation in the 
#' limma::eBayes function (default = TRUE)
#' @param keep.cols Which columns should be copied from data_raw to 
#' data_peptides
#' @param rename.output named vector of columns to rename in the output (e.g. 
#' #' c(name = "new_name"))
#' @param full.output Include full limma data in output?
#' @param translate.concentrations replace numeric concentrations with strings 
#' with unit 
#'
#' @returns list of concentration-grouped list results at peptide ($data_peptides) and replicates 
#' ($data_replicates) level 
#' @export
#'
#' @examples 
#' \dontrun{
#' data_limma <- analyze_data_dose_response_concentrations(data_processed,
#'   min_fraction_per_group = 1, 
#'   norm.method = "none", 
#'   p.adjust.method = "BH",
#'   p.threshold = 0.01,
#'   fc.threshold = log2(1.2), 
#'   keep.cols = c("Protein.Group", 
#'    "Genes", 
#'    "First.Protein.Description"), 
#'   rename.output = c(estimate = "log2.fc"))
#' }
analyze_data_dose_response_concentrations <- function(data_processed, 
                                                      filter_rows, 
                                                      min_fraction_per_group = 1, 
                                                      norm.method = "none", 
                                                      p.adjust.method = "BH", 
                                                      p.threshold = 0.01, 
                                                      fc.threshold = log2(1.2), 
                                                      eBayes.trend = TRUE, 
                                                      keep.cols = c(), 
                                                      rename.output = c(estimate = "log2.fc"), 
                                                      full.output = F, 
                                                      silent = T, 
                                                      translate.concentrations = T) {
  
  if (all(data_processed$data_processed$Concentration != 0)) 
    stop("There is no 0-concentration group in the data. Make sure the concentration is set to 0.")
  
  data_raw_new <- data_processed$data_processed %>% 
    dplyr::select(-c(do.fit, found_in, Quant, Replicate)) %>% 
    dplyr::left_join(data_processed$data_raw %>% 
                       dplyr::select(Peptides, 
                                     Samples, 
                                     Quant, 
                                     dplyr::any_of(keep.cols)) %>% 
                       mutate(Samples = (data_processed$data_processed %>% 
                                           dplyr::distinct(Samples) %>% 
                                           dplyr::pull(Samples))[Samples]), 
                     by = c("Peptides", "Samples"))
  
  data_observations <- data_processed$data_processed %>% 
    distinct(Samples, Concentration) %>% 
    dplyr::arrange(desc(Concentration), Samples)
  
  if (translate.concentrations) {
    data_observations <- data_observations %>% 
      dplyr::mutate(Concentration = .translate_concentrations(Concentration))
    control_group <- "0M"
  } else {
    data_observations <- data_observations %>% 
      dplyr::mutate(Concentration = as.character(Concentration))
    control_group <- 0
  }
  
  
  
  data_results <- map(
    data_observations %>% 
      dplyr::filter(Concentration != 0, Concentration != "0M") %>% 
      pull(Concentration) %>% 
      unique() %>% 
      subset(. > 0) %>% 
      setNames(., .), 
    
    \(x) {
      
      data_observations_filtered <- data_observations %>% 
        filter(Concentration %in% c(control_group, x)) %>% 
        arrange(Concentration)
      
      data_processed_grouped <- process_data_grouped(data_raw = data_raw_new, 
                                                     sample_names = data_observations_filtered %>% 
                                                       pull(Samples, Samples), 
                                                     sample_groups = data_observations_filtered %>% 
                                                       pull(Concentration, Samples), 
                                                     peptide.column = "Peptides", 
                                                     sample.column = "Samples", 
                                                     quant.column = "Quant", 
                                                     filter_rows = {{filter_rows}}, 
                                                     min_fraction_per_group = min_fraction_per_group, 
                                                     norm.method = norm.method, 
                                                     silent = silent)
      
      data_results <- analyze_data_grouped(data_processed = data_processed_grouped, 
                                           conditions = c(control_group, x), 
                                           p.adjust.method = p.adjust.method, 
                                           p.threshold = p.threshold, 
                                           fc.threshold = fc.threshold, 
                                           eBayes.trend = eBayes.trend, 
                                           keep.cols = keep.cols, 
                                           rename.output = rename.output, 
                                           full.output = full.output, 
                                           silent = silent) 
      
      return(data_results)
    })
  
  return(data_results)
  
}
