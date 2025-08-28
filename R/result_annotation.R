#' classify protein and peptide hits in group-wise PELSA experiment 
#'
#' @param data_peptides 
#' @param p.value_column 
#' @param fc_column 
#' @param p.threshold.protein 
#' @param fc.threshold.protein 
#' @param p.threshold.peptide 
#' @param fc.threshold.peptide 
#' @param protein.id_column 
#'
#' @returns
#' @export
#'
#' @examples
classify_peptide_protein_hits <- function(data_peptides, 
                                          p.value_column = "p.value", 
                                          fc_column = "log2.fc", 
                                          p.threshold.protein = 10e-4, 
                                          fc.threshold.protein = 0, 
                                          p.threshold.peptide = 10e-2, 
                                          fc.threshold.peptide = 0.3, 
                                          protein.id_column = "Protein.Group") {
  
  data_annotated <- data_peptides %>% 
    dplyr::mutate(
      # Classify protein regulation 
      regulation_protein = dplyr::case_when(
        !!rlang::sym(p.value_column) < p.threshold.protein & 
          !!rlang::sym(fc_column) > fc.threshold.protein ~ "destabilized", 
        !!rlang::sym(p.value_column) < p.threshold.protein & 
          !!rlang::sym(fc_column) < -fc.threshold.protein ~ "stabilized", 
        .default = "unchanged")) %>% 
    # Summarize protein level regulation to all peptides of a protein from most 
    # significant peptide (lowest p-value)
    dplyr::mutate(regulation_protein = regulation_protein[which.min(!!rlang::sym(p.value_column))], 
                  .by = dplyr::all_of(protein.id_column)) %>% 
    ## Assign peptide level regulation taking protein level into account 
    dplyr::mutate(regulation_peptide = dplyr::case_when(
      # Unchanged proteins 
      regulation_protein == "unchanged" ~ "unchanged", 
      # Peptides significant on protein level
      ## Destabilized protein 
      !!rlang::sym(p.value_column) < p.threshold.protein & 
        !!rlang::sym(fc_column) > fc.threshold.protein ~ "destabilized", 
      ## Stabilized protein 
      !!rlang::sym(p.value_column) < p.threshold.protein & 
        !!rlang::sym(fc_column) < -fc.threshold.protein ~ "stabilized", 
      # Destabilized peptides 
      !!rlang::sym(p.value_column) < p.threshold.peptide & 
        !!rlang::sym(fc_column) > fc.threshold.peptide & 
        regulation_protein %in% c("stabilized", "destabilized") ~ "destabilized", 
      # Stabilized peptides 
      !!rlang::sym(p.value_column) < p.threshold.peptide & 
        !!rlang::sym(fc_column) < -fc.threshold.peptide & 
        regulation_protein %in% c("stabilized", "destabilized") ~ "stabilized", 
      # Unchanged peptides 
      abs(!!rlang::sym(fc_column)) < fc.threshold.peptide ~ "unchanged", 
      # Unknown peptides 
      abs(!!rlang::sym(fc_column)) > fc.threshold.peptide & 
        !!rlang::sym(p.value_column) >= p.threshold.peptide ~ "unknown", 
      .default = "none")) %>% 
    dplyr::mutate(regulation = regulation_peptide)
  
  return(data_annotated)
  
}


#' Add peptide position based on stripped sequence column 
#'
#' @param x tibble containing sequence and protein Id columns 
#' @param sequences named vector containing protein sequences 
#' @param sequence.col column containing peptide sequence 
#' @param protein.col column containing protein Ids
#' @param keep.peptide_length Should the peptide length column be kept in the 
#' data? 
#' @param .after position of the sequence position columns 
#'
#' @returns
#' @export
#'
#' @examples
mutate_peptide_position <- function(x, 
                                    sequences = c(A = "ABC"), 
                                    sequence.col = "Stripped.Sequence", 
                                    protein.col = "Protein.Group", 
                                    keep.peptide_length = T, 
                                    .after = sequence.col) {
  
  # Start with removing possible modifications 
  data_pos <- x %>% 
    dplyr::mutate(plain_sequence = 
                    stringr::str_remove_all(!!rlang::sym(sequence.col), 
                                            "\\(.*?\\)")) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(from = NA_integer_, 
                  to = NA_integer_, 
                  peptide_length = nchar(plain_sequence), 
                  from = regexpr(plain_sequence, as.character(sequences[unlist(strsplit(!!rlang::sym(protein.col), ";"))[1]]))[1], 
                  to = from + peptide_length - 1, 
                  .after = dplyr::all_of(.after)) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-plain_sequence)
  
  if (!keep.peptide_length) 
    data_pos <- data_pos %>% 
      dplyr::select(-peptide_length)
  
  
  return(data_pos)
  
}


#' Download UniProt data for given protein accessions and data fields 
#' (see available fields with UniProt_fields())
#'
#' @param accession vector of UniProt accessions 
#' @param fields UniProt data fields to query
#' @param max.query maximum number of accessions to query at once; if the 
#' the number exceeds max.query, the query is split up in multiple parts 
#'
#' @returns
#' @export
#'
#' @examples
get_UniProt_data <- function(accession, 
                             fields = c("accession", 
                                        "gene_names", 
                                        "organism_name"), 
                             max.query = 1000) {
  
  # Add accession as field
  if (!"accession" %in% fields) fields <- c("accession", fields)
  
  # Check accession for ;
  if (any(stringr::str_detect(accession, ";"))) 
    stop("There are accessions contanining a semicolon (;);, please remove or correct protein groups with multiple Ids.")
  
  # Only query unique accessions 
  accession_query <- unique(accession)
  
  
  # Formulate query/ies and download data 
  if (length(accession_query) > max.query) {
    
    l <- length(accession_query)
    from <- seq(1, l, max.query)
    to <- c(seq(1, l, max.query)[-1] - 1, l)
    
    data_download <- purrr::map2(from, to, 
                                 \(from, to) get_UniProt_data(accession_query[from:to], 
                                                              fields = fields, 
                                                              max.query = max.query)) %>% 
      bind_rows()
    
  } else {
    
    query_url <- paste0("https://rest.uniprot.org/uniprotkb/stream?", 
                        "format=tsv", 
                        "&fields=", 
                        paste(fields, collapse = "%2C"), 
                        "&query=", 
                        paste(
                          paste0("accession%3A", accession_query), 
                          collapse = "+OR+"))
    
    data_download <- vroom::vroom(query_url, 
                                  delim = "\t", 
                                  col_types = readr::cols())
    
  }
  
  # Merge given accessions and downloaded data 
  data_output <- dplyr::left_join(tibble::tibble(Entry = accession), 
                                  data_download, 
                                  by = "Entry")
  
  # Return tibble with accessions as Entry and data columns 
  return(data_output)
  
}


#' Download UniProt data for given protein accessions, taxonomy identifiers and 
#' data fields (faster for 1000s of proteins; see available fields with 
#' UniProt_fields())
#'
#' @param accession 
#' @param fields 
#' @param taxon_id 
#'
#' @returns
#' @export
#'
#' @examples
get_UniProt_data_1o <- function(accession, 
                                fields = c("accession", 
                                           "gene_names", 
                                           "organism_name"), 
                                taxon_id = c(human = 9606, 
                                             mouse = 10900, 
                                             E.coliK12 = 83333)) {
  
  # Check organism identifier
  if (length(taxon_id) > 1) 
    warning("More than one taxon_id provided, only using the first one.")
  
  # Add accession as field
  if (!"accession" %in% fields) fields <- c("accession", fields)
  
  
  # Formulate query and download data 
  query_url <- paste0("https://rest.uniprot.org/uniprotkb/stream?", 
                      "format=tsv", 
                      "&fields=", 
                      paste(fields, collapse = "%2C"), 
                      "&query=%28model_organism%3A", 
                      taxon_id[1], "%29")
  
  data_download <- vroom::vroom(query_url, 
                                delim = "\t", 
                                col_types = readr::cols())
  
  
  # Merge given accessions and downloaded data 
  if (hasArg(accession)) {
    data_output <- dplyr::left_join(tibble::tibble(Entry = accession), 
                                    data_download, 
                                    by = "Entry")
  } else {
    data_output <- data_download
  }
  
  
  # Return tibble with accessions as Entry and data columns 
  return(data_output)
  
}


#' Add protein annotations to dose-response PELSA results 
#'
#' @param data results list or peptides data frame from 
#' analyze_data_dose_response() 
#' @param col_names named vector specifying which columns should be renamed 
#' @param fields UniProt fields to be downloaded (see UniProt_fields())
#' @param max.query maximum number of proteins per UniProt query 
#' @param taxon_id for single species data, taxon id to speed up UniProt 
#' download 
#' @param data_UniProt predownloaded UniProt annotation data 
#' @param add_anno2peptides Should protein annotations be added tot the 
#' peptides data frame?
#' @param keep.peptide_length Should the peptide length be kept in the data? 
#'
#' @returns
#' @export
#'
#' @examples
annotate_results_dr <- function(data, 
                                col_names = c(Peptide.Id = "Peptides", 
                                              Protein.Description = "First.Protein.Description", 
                                              pEC50 = "pEC50_mean", 
                                              CV_pEC50 = "pEC50_CV"), 
                                fields = c("length", 
                                           "ft_binding", 
                                           "xref_pdb_full", 
                                           "sequence"), 
                                max.query = 100, 
                                taxon_id = NULL, 
                                data_UniProt = NULL, 
                                add_anno2peptides = F, 
                                keep.peptide_length = T) {
  
  list_output <- 
    list(data_peptides = 
           tibble::tibble(Peptide.Id = NA_character_,
                          Protein.Group = NA_character_,
                          Genes = NA_character_,
                          Protein.Description = NA_character_,
                          regulation = NA_character_,
                          # this 
                          pEC50 = NA_real_,
                          pEC50_CV = NA_real_,
                          # or 
                          log2.fc = NA_real_, 
                          p.value = NA_real_, 
                          from = NA_integer_,
                          to = NA_integer_),
         data_raw = tibble::tibble(Peptide.Id = NA_character_,
                                   replicate = NA_character_,
                                   concentration = NA_real_,
                                   value = NA_real_),
         data_protein_annotation = tibble::tibble(Protein.Group = NA_character_))
  
  if (tibble::is_tibble(data)) {
    
    list_output$data_peptides <- data %>% 
      dplyr::rename(dplyr::any_of(col_names)) %>%
      dplyr::mutate(pEC50_signed = ifelse(regulation == "stabilized", 1, -1) * pEC50) 
    
  } else {
    
    list_output$data_peptides <- data$data_peptides %>% 
      dplyr::rename(dplyr::any_of(col_names)) %>% 
      dplyr::mutate(pEC50_signed = ifelse(regulation == "stabilized", 1, -1) * pEC50) 
    
    list_output$data_raw <- data$data_replicates
    
  }
  
  # Get UniProt annotations 
  if (is.null(data_UniProt)) {
    if (!is.null(taxon_id)) {
      data_UniProt <- list_output$data_peptides %>% 
        dplyr::pull(Protein.Group) %>% 
        unique() %>% 
        str_extract("^.+?(?=;|$)") %>% 
        get_UniProt_data_1o(fields = fields, 
                            taxon_id = taxon_id)
    } else {
      data_UniProt <- list_output$data_peptides %>% 
        dplyr::pull(Protein.Group) %>% 
        unique() %>% 
        str_extract("^.+?(?=;|$)") %>% 
        get_UniProt_data(fields = fields, 
                         max.query = max.query)
    }
  }
  
  
  list_output$data_protein_annotation <- list_output$data_peptides %>% 
    dplyr::distinct(Protein.Group) %>% 
    dplyr::left_join(data_UniProt, 
                     by = c("Protein.Group" = "Entry"))
  
  
  
  list_output$data_peptides <- list_output$data_peptides %>% 
    mutate_peptide_position(sequences = dplyr::pull(list_output$data_protein_annotation, 
                                                    Sequence, Protein.Group), 
                            sequence.col = "Peptide.Id", 
                            keep.peptide_length = keep.peptide_length)
  
  
  # Transfer protein annotation to peptide data frame 
  if (add_anno2peptides) {
    list_output$data_peptides <- list_output$data_peptides %>% 
      dplyr::left_join(list_output$data_protein_annotation, 
                       by = "Protein.Group")
  }
  
  return(list_output)
  
}


#' Add protein annotations to group-wise PELSA results 
#'
#' @param data results list or peptides data frame from analyze_data_grouped() 
#' @param col_names named vector specifying which columns should be renamed 
#' @param fields UniProt fields to be downloaded (see UniProt_fields())
#' @param max.query maximum number of proteins per UniProt query 
#' @param taxon_id for single species data, taxon id to speed up UniProt 
#' download 
#' @param data_UniProt predownloaded UniProt annotation data 
#' @param add_anno2peptides Should protein annotations be added tot the 
#' peptides data frame?
#' @param keep.peptide_length Should the peptide length be kept in the data? 
#'
#' @returns
#' @export
#'
#' @examples
annotate_results_gw <- function(data, 
                                col_names = c(Peptide.Id = "Peptides", 
                                              Protein.Description = "First.Protein.Description", 
                                              log2.fc = "Log2FC", 
                                              p.value = "P.Value"), 
                                fields = c("length", 
                                           "ft_binding", 
                                           "xref_pdb_full", 
                                           "sequence"), 
                                max.query = 100, 
                                taxon_id = NULL, 
                                data_UniProt = NULL, 
                                add_anno2peptides = F, 
                                keep.peptide_length = T) {
  
  list_output <- 
    list(data_peptides = 
           tibble::tibble(Peptide.Id = NA_character_,
                          Protein.Group = NA_character_,
                          Genes = NA_character_,
                          Protein.Description = NA_character_,
                          regulation = NA_character_,
                          # this 
                          pEC50 = NA_real_,
                          pEC50_CV = NA_real_,
                          # or 
                          log2.fc = NA_real_, 
                          p.value = NA_real_, 
                          from = NA_integer_,
                          to = NA_integer_),
         data_raw = tibble::tibble(Peptide.Id = NA_character_,
                                   replicate = NA_character_,
                                   concentration = NA_real_,
                                   value = NA_real_),
         data_protein_annotation = tibble::tibble(Protein.Group = NA_character_))
  
  if (tibble::is_tibble(data)) {
    
    list_output$data_peptides <- data %>% 
      dplyr::rename(dplyr::any_of(col_names))
  } else {
    
    list_output$data_peptides <- data$data_peptides %>% 
      dplyr::rename(dplyr::any_of(col_names))
    
    list_output$data_raw <- data$data_replicates
    
  }
  
  # Get UniProt annotations 
  if (is.null(data_UniProt)) {
    if (!is.null(taxon_id)) {
      data_UniProt <- list_output$data_peptides %>% 
        dplyr::pull(Protein.Group) %>% 
        unique() %>% 
        str_extract("^.+?(?=;|$)") %>% 
        get_UniProt_data_1o(fields = fields, 
                            taxon_id = taxon_id)
    } else {
      data_UniProt <- list_output$data_peptides %>% 
        dplyr::pull(Protein.Group) %>% 
        unique() %>% 
        str_extract("^.+?(?=;|$)") %>% 
        get_UniProt_data(fields = fields, 
                         max.query = max.query)
    }
  }
  
  list_output$data_protein_annotation <- list_output$data_peptides %>% 
    dplyr::distinct(Protein.Group) %>% 
    dplyr::left_join(data_UniProt, 
                     by = c("Protein.Group" = "Entry"))
  
  
  list_output$data_peptides <- list_output$data_peptides %>% 
    mutate_peptide_position(sequences = pull(list_output$data_protein_annotation, 
                                             Sequence, Protein.Group), 
                            sequence.col = "Peptide.Id", 
                            keep.peptide_length = keep.peptide_length)
  
  
  # Transfer protein annotation to peptide data frame 
  if (add_anno2peptides) {
    list_output$data_peptides <- list_output$data_peptides %>% 
      dplyr::left_join(list_output$data_protein_annotation, 
                       by = "Protein.Group")
  }
  
  return(list_output)
  
}

