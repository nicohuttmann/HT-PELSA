#' Process the imported DIA-NN precursor data 
#'
#'
#' @param data_raw tibble in long or wide format (e.g. from 
#' import_data_long()/_wide())
#' @param sample_names named vector containing sample names as value and 
#' original sample names as names (specifies which samples to use and renames 
#' them)
#' @param sample_groups named vector containing group information as value and 
#' original sample names as names 
#' @param peptide.column column containing peptide IDs to which the data should 
#' be summarized/collapse (by sum())
#' @param sample.column column containing sample names (when working 
#' with long format data)
#' @param quant.column column containing quantitative information (when working 
#' with long format data)
#' @param filter_rows tidy expression to filter the data e.g. filter_rows = 
#' Lib.Q.Value < 0.01
#' @param min_fraction_per_group minimum fraction of values per replicate to use 
#' (e.g. 0.75 = min. 3 of 4 replicates with valid value)
#' @param norm.method which method to use for normalization 
#' @param silent Turn off feedback from the function?
#'
#' @importFrom magrittr %>%
#' @returns list of processed peptide data ($data_processed) and used raw data 
#' ($data_raw)
#' @export
#'
#' @examples
#' # see import_data_long()/_wide() output 
process_data_grouped <- function(data_raw, 
                                 sample_names, 
                                 sample_groups, 
                                 peptide.column = "Stripped.Sequence", 
                                 sample.column = "Run", 
                                 quant.column = "Precursor.Normalised", 
                                 filter_rows, 
                                 min_fraction_per_group = 1, 
                                 norm.method = "none", 
                                 silent = F) {
  
  # Check input
  if (!hasArg(data_raw)) 
    stop("Please provide the imported data with <data_raw>.")
  
  # Check names 
  if (!hasArg(sample_names)) 
    stop("Please define the samples to use and the final sample names with <sample_names>.")
  
  # Check groups 
  if (!hasArg(sample_groups) || length(unique(names(sample_groups))) < 2) 
    stop("Please define the sample groups by changing the names of the <sample_groups> argument vector. There must be minimum two groups.")
  
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
  
  
  # Rename data columns 
  data_filtered_r <- data_filtered %>% 
    dplyr::mutate(Samples = sample_names[!!rlang::sym(sample.column)], 
                  .before = 1) %>% 
    dplyr::filter(!is.na(Samples)) %>% 
    dplyr::arrange(Samples) %>% 
    dplyr::mutate(Peptides = !!rlang::sym(peptide.column), 
                  .after = 1) %>% 
    dplyr::mutate(Quant = !!rlang::sym(quant.column), 
                  .after = 2)
  
  
  # Collapse peptides with peptide.column
  message(paste0('Your peptides will be summarized at the level of "', 
                 peptide.column, 
                 '".'))
  
  data_quant <- data_filtered_r %>%  
    dplyr::summarise(Quant = sum(Quant, na.rm = T), 
                     .by = c("Samples", "Peptides"))
  
  
  
  # Add group information and filter for fraction within group 
  data_quant_c <- data_quant %>%
    # Remove missing values 
    dplyr::filter(!is.na(Quant)) %>% 
    # Add Group column 
    dplyr::mutate(Group = unname(setNames(sample_groups, 
                                          sample_names[names(sample_groups)])[Samples]), 
                  .after = "Samples") %>% 
    # Count values per peptide and group 
    dplyr::mutate(n = sum(Quant > 0), 
                  .by = c("Peptides", "Group")) %>% 
    # Calculate fraction of each peptide in group 
    dplyr::mutate(p = n / max(n), .by = "Group") %>% 
    # Remove peptides below min_fraction_per_group threshold 
    dplyr::filter(all(p >= min_fraction_per_group), 
                  .by = "Peptides") %>% ##### Check this 
    dplyr::mutate(n_peptide = length(Samples), 
                  .by = "Peptides") %>% 
    # Remove peptides again based on unequal group sizes 
    dplyr::filter(n_peptide == max(n_peptide)) %>% 
    # pivot to wide data format with samples as rows 
    tidyr::pivot_wider(id_cols = c("Samples", "Group"), 
                       names_from = "Peptides", 
                       values_from = "Quant")
  
  
  # Check sample groups 
  # sample_groups_m <- tibble::tibble(Samples = unname(sample_names), 
  #                                   Group = setNames(sample_groups, 
  #                                                    sample_names[names(sample_groups)])[Samples], 
  #                                   sample_name_original = names(sample_names)) 
  
  # Prepare output list 
  list_output <- list(data_processed = data_quant_c, 
                      #sample_groups = sample_groups_m, 
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
                                    collapse = "\",\n\t\""), "\")"), 
               ",\n")
    }
    
    
    # Message 
    message()
    message("Your data was successfully processed and can be used for the limma analysis. Please copy the code below.")
    message()
    
    to_cat <- paste0('# Analyze peptides with limma\n', 
                     'data_limma <- analyze_data_grouped(data_processed,\n', 
                     'conditions = c("', 
                     unique(sample_groups)[1], 
                     '", "', 
                     unique(sample_groups)[2], 
                     '"),\n', 
                     '\tp.adjust.method = "BH",\n', 
                     '\tp.threshold = 0.01,\n', 
                     '\tfc.threshold = log2(1.2),\n', 
                     potential_columns, 
                     '\trename.output = c(estimate = "log2.fc"))\n')
    
    for (i in unlist(strsplit(to_cat, split = ""))) {
      cat(i)
      Sys.sleep(0.001)
    }
    
    # Return output list 
    return(list_output)
    
  }
}


#' Peptide-based limma analysis
#'
#' @param data_processed output list or tibble from process_data_grouped()
#' @param conditions character vector specifying conditions to compare 
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
#' @param silent Turn off feedback from the function?
#'
#' @importFrom magrittr %>%
#' @returns tibble of limma results of full list of limma data incl. results 
#' @export
#'
#' @examples 
#' \dontrun{
#' data_limma <- analyze_data_grouped(data_processed,
#'   conditions = c("Control", "Treatment"),
#'   p.adjust.method = "BH",
#'   p.threshold = 0.01,
#'   fc.threshold = log2(1.2), 
#'   keep.cols = c("Protein.Group", 
#'    "Genes", 
#'    "First.Protein.Description"), 
#'   rename.output = c(estimate = "log2.fc"))
#' }
analyze_data_grouped <- function(data_processed, 
                                 conditions = c("Control", "Treatment"), 
                                 p.adjust.method = "BH", 
                                 p.threshold = 0.01, 
                                 fc.threshold = log2(1.2), 
                                 eBayes.trend = TRUE, 
                                 keep.cols = c(), 
                                 rename.output = c(estimate = "log2.fc"), 
                                 full.output = F, 
                                 silent = F) {
  
  # Check if dataset is given 
  if (hasArg(data_processed) && all(c("data_processed") %in% 
                                    names(data_processed))) {
    data_quant <- data_processed[["data_processed"]]
  } else {
    stop("Please use the output of the process_data_grouped() function as input for <data_processed> or reshape your data to match the input.")
  }
  
  
  # Check packages 
  
  ## limma 
  if (!suppressPackageStartupMessages(require(limma, quietly = T, warn.conflicts = F))) {
    if(menu(c("Yes, please install the limma package now.", 
              "No, I will install it myself."), 
            title = ("The limma package is not installed yet, but required for the analysis. Would you like to install it?")) == 1) {
      if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      
      BiocManager::install("limma")
    } else {
      return(invisible(FALSE))
    }
  }
  
  ## biobroom 
  if (!suppressPackageStartupMessages(require(biobroom, quietly = T, warn.conflicts = F))) {
    if(menu(c("Yes, please install the biobroom package now.", 
              "No, I will install it myself."), 
            title = ("The biobroom package is not installed yet, but required for the analysis. Would you like to install it?")) == 1) {
      if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      
      BiocManager::install("biobroom")
    } else {
      return(invisible(FALSE))
    }
  }
  
  
  
  # Check samples and groups 
  if (!"Group" %in% names(data_quant)) {
    stop('The data needs a sample "Group" column.')
  }
  
  
  
  # Add limma input 
  results_list <- list()
  
  
  ## Add base data frame 
  results_list[["data"]] <- data_quant %>% 
    tidyr::pivot_longer(-c(1, 2)) %>% 
    dplyr::mutate(value = log2(value)) %>% 
    tidyr::pivot_wider()
  
  ## Add eset 
  results_list[["eset"]] <- results_list[["data"]] %>%
    dplyr::select(-Group) %>% 
    dplyr::rename(rowname = 1) %>% 
    tibble::column_to_rownames() %>%
    as.matrix() %>% 
    t() %>% 
    Biobase::ExpressionSet()
  
  
  
  ## Describe experimental groups
  sample_groups_mod <- 
    setNames(c("Control", "Treatment"), conditions)[data_quant$Group]
  
  # Add to data 
  results_list$data <- results_list$data %>% 
    dplyr::mutate(Group = sample_groups_mod, .after = 1)
  
  results_list[["design"]] <- 
    model.matrix(
      ~0+factor(sample_groups_mod, 
                unique(sample_groups_mod)))
  
  ## Correct column names
  colnames(results_list[["design"]]) <- unique(sample_groups_mod)
  
  
  # Do first linear model fit
  results_list[["fit"]] <- limma::lmFit(results_list[["eset"]], 
                                        results_list[["design"]])
  
  # Describe comparisons (contrasts)
  results_list[["contrast.matrix"]] <- 
    limma::makeContrasts(Treatment-Control, 
                         levels=results_list[["design"]])
  
  # Compute contrasts between groups
  results_list[["fit2"]] <- 
    limma::contrasts.fit(fit = results_list[["fit"]], 
                         contrasts = results_list[["contrast.matrix"]])
  
  # Apply empirical ebayes to estimate t-statistic and calculate p-values
  results_list[["fit2_eBayes"]] <- limma::eBayes(results_list[["fit2"]], 
                                                 trend = eBayes.trend)
  
  
  # Define thresholds 
  results_list[["par"]] <- list(p.threshold = p.threshold, 
                                p.adjust.method = p.adjust.method, 
                                fc.threshold = fc.threshold)
  
  # Extract results 
  results_list[["results"]] <- suppressWarnings(
    biobroom::tidy.MArrayLM(results_list[["fit2_eBayes"]]) %>% 
      # rename contrasts 
      dplyr::mutate(term = paste(rev(conditions), collapse = " - ")) %>% 
      dplyr::rename(Peptides = gene) %>% 
      dplyr::arrange(p.value) %>% 
      # Adjust p-values
      dplyr::mutate(log10.p.value = -log10(p.value), 
                    p.adjust = 
                      p.adjust(p.value, 
                               method = results_list[["par"]]$p.adjust.method), 
                    .after = "p.value") %>% 
      # Apply thresholds
      dplyr::mutate(regulation = dplyr::case_when(
        p.adjust < results_list[["par"]]$p.threshold & 
          estimate >= results_list[["par"]]$fc.threshold ~ "destabilized", 
        p.adjust < results_list[["par"]]$p.threshold & 
          estimate <= - results_list[["par"]]$fc.threshold ~ "stabilized", 
        .default = "none")))
  
  
  # Plot results 
  results_list[["p"]] <- results_list[["results"]] %>%
    dplyr::filter(!is.na(regulation)) %>%
    ggplot2::ggplot(ggplot2::aes(x = estimate,
                                 y = -log10(p.value),
                                 color = regulation)) +
    ggplot2::geom_point(shape = 16, alpha = 0.5) +
    ggplot2::scale_color_manual(values = c(none = "grey",
                                           destabilized = "red",
                                           stabilized = "blue")) +
    ggplot2::theme_classic() +
    ggplot2::coord_cartesian(expand = F)
  
  
  if (!silent) {
    print(results_list[["p"]])
  }
  
  
  # Rename results columns 
  if (hasArg(rename.output) && 
      all(names(rename.output) %in% names(results_list[["results"]]))) {
    
    results_list[["results"]] <- results_list[["results"]] %>% 
      dplyr::rename_with(\(x) ifelse(x %in% names(rename.output), 
                                     rename.output[x], 
                                     x))
  }
  
  
  
  # Transfer annotation columns 
  if (length(keep.cols) > 0 && 
      "data_raw" %in% names(data_processed) && 
      any(keep.cols %in% names(data_processed$data_raw))) {
    
    results_list[["results"]] <- data_processed$data_raw %>% 
      dplyr::distinct(Peptides, .keep_all = T) %>% 
      dplyr::select(any_of(c("Peptides", keep.cols))) %>% 
      dplyr::right_join(results_list[["results"]], 
                        by = "Peptides")
  }
  
  
  # Success message
  if (!silent) 
    message(paste0("Good job, the limma analysis is done. ", 
                   sum(results_list$results$regulation == "stabilized"), 
                   " peptides were identified to be stabilized and ", 
                   sum(results_list$results$regulation == "destabilized"), 
                   " as destabilized by ", 
                   conditions[2], 
                   "."))
  
  
  # Remove intermediate data 
  if (!full.output) {
    results_list <- results_list[c("results", "data", "p")]
    names(results_list) <- c("data_peptides", "data_replicates", "p")
  } else {
    results_list <- results_list[c("results", 
                                   "data", 
                                   "p", 
                                   "eset",
                                   "design",
                                   "fit",
                                   "contrast.matrix",
                                   "fit2",
                                   "fit2_eBayes",
                                   "par")]
    names(results_list) <- c("data_peptides", 
                             "data_replicates", 
                             "p", 
                             "eset",
                             "design",
                             "fit",
                             "contrast.matrix",
                             "fit2",
                             "fit2_eBayes",
                             "par")
  }
  
  return(results_list)
  
}

