


# Colors for peptides in plots 
peptide_colors <- c(underquantified = "#b0b0b0",  #"grey40", 
                    
                    unchanged = "#666666", #"grey30",
                    unresponsive = "#666666", #"grey30", 
                    none = "#666666", #"grey30", #4d4d4d
                    
                    regular = "#666666", 
                    unknown = "#734595", #"grey30",
                    
                    destabilized = "#FF0000", # "red", 
                    stabilized = "#0000FF", # "blue", 
                    highlight = "#18974C")  



#' Title
#'
#' @param p 
#' @param data_sites 
#' @param name 
#' @param start 
#' @param end 
#' @param merge_name_position 
#' @param protein_width 
#' @param base_size 
#' @param legend.position 
#' @param colors_manual 
#'
#' @returns
#' @export
#'
#' @examples
plot_protein_seq <- function(protein_range = c(1, 100), 
                             protein_width = 0.6, 
                             x_lab = "protein", 
                             y_lab = "", 
                             color = "#b0b0b0", 
                             alpha = 0.5, 
                             base_size = 16) {
  
  p_protein <- ggplot2::ggplot() + 
    ggplot2::annotate("rect", 
                      xmin = protein_range[1], 
                      xmax = protein_range[2], 
                      ymin = - protein_width / 2, 
                      ymax =   protein_width / 2, 
                      fill = color, 
                      alpha = alpha) + 
    ggplot2::theme_classic(base_size = base_size) + 
    ggplot2::scale_x_continuous(expand = c(0, 0)) + 
    ggplot2::coord_cartesian(xlim = protein_range) + 
    ggplot2::labs(x = x_lab, y = y_lab)
  
  return(p_protein)
  
} 



#' Title
#'
#' @param p 
#' @param data_sites 
#' @param name 
#' @param start 
#' @param end 
#' @param merge_name_position 
#' @param protein_width 
#' @param base_size 
#' @param legend.position 
#' @param colors_manual 
#'
#' @returns
#' @export
#'
#' @examples
add_plot_protein_sites <- function(p, 
                                   data_sites, 
                                   name = "site", 
                                   start = "start", 
                                   end = "end", 
                                   add.length = 0, 
                                   merge_name_position = F, 
                                   protein_width = 0.6, 
                                   base_size = 16, 
                                   legend.position = "bottom", 
                                   colors_manual) {
  
  data_sites_p <- data_sites %>% 
    mutate(ymin = -1, 
           ymax = 1) %>% 
    rename(name = all_of(name), 
           start = all_of(start), 
           end = all_of(end)) %>% 
    mutate(position = ifelse(start != end, 
                             paste0(start, 
                                    "-", 
                                    end), 
                             start)) %>%  
    arrange(start, end, name)
  
  # if (merge_name_position) 
  #   data_sites_p <- data_sites_p %>% 
  #     mutate(name = ifelse(start != end, 
  #                            paste0(name, 
  #                                   " (", 
  #                                   start, 
  #                                   "-", 
  #                                   end, 
  #                                   ")"), 
  #                            paste0(name, 
  #                                   " (", 
  #                                   start, 
  #                                   ")")), 
  #            .by = any_of(c("name", "start", "end")))
  
  
  if (merge_name_position) 
    data_sites_p <- data_sites_p %>% 
      mutate(name = paste0(unique(name), " (", 
                           paste(position, collapse = "/"), 
                           ")"), 
             .by = "name")
  
  
  for (i in 1:nrow(data_sites_p)) {
    if (i == 1) 
      y_height <- protein_width
    else 
      y_height <- c(y_height, 
                    protein_width * (1 - sum(data_sites_p$start[i] < 
                                               data_sites_p$end[1:(i-1)]) / 6))
  }
  
  data_sites_p <- data_sites_p %>% 
    mutate(ymin = - y_height / 2, 
           ymax =   y_height / 2, 
           name = factor(name, 
                         levels = unique(name))) #%>% 
  # add data for rounded corners 
  # mutate(x = mean(c(start, end)), 
  #        y = 0, 
  #        width = end - start, 
  #        height = 2 * ymax, 
  #        .by = c("name", "start", "end"))
  
  
  p_protein <- p + 
    # annotate("rect", 
    #          xmin = protein_range[1], 
    #          xmax = protein_range[2], 
    #          ymin = - protein_width / 2, 
    #          ymax =   protein_width / 2, 
    #          fill = "grey") + 
    geom_rect(aes(xmin = start - add.length, 
                  xmax = end + add.length, 
                  ymin = ymin, 
                  ymax = ymax, 
                  fill = name), 
              data_sites_p) + 
    #theme_void(base_size = base_size) + 
    #theme(legend.position = ) +
    guides(fill = guide_legend(position = legend.position, ncol = 1)) 
  #scale_fill_embl() + 
  
  if (hasArg(colors)) 
    p_protein <- p_protein + 
    scale_fill_manual(values = colors_manual)
  else 
    p_protein <- p_protein + 
    scale_fill_embl()
  
  
  return(p_protein)
  
} 



#' Title
#'
#' @param data_peptides 
#' @param p
#' @param yvalue 
#' @param color 
#' @param linewidth 
#' @param Start 
#' @param End 
#' @param title 
#' @param base_size 
#' @param min_y_range 
#' @param protein_width 
#' @param protein_range 
#'
#' @returns
#' @export
#'
#' @examples
plot_peptides_seq2 <- function(data_peptides, 
                               p, 
                               yvalue = "estimate", 
                               color = "regulation", 
                               color.scale = NULL, 
                               linewidth = "peptide_group", 
                               alpha = "peptide_group", 
                               Start = "Start", 
                               End = "End", 
                               title = "", 
                               base_size = 16, 
                               min_y_range = 1, 
                               protein_width = 0.6, 
                               protein_range, 
                               add.labels = F, 
                               label = "regulation", 
                               label.size = 5, 
                               nudge_x = 0, 
                               nudge_y = 0, 
                               hjust = 0,
                               vjust = 0, 
                               direction = "both", 
                               min.segment.length = 2) {
  
  
  # data_p <- data_peptides %>% 
  #   rename(yvalue = all_of(yvalue), 
  #          Start = all_of(Start), 
  #          End = all_of(End), 
  #          color = all_of(color), 
  #          alpha = all_of(alpha), 
  #          linewidth = all_of(linewidth), 
  #          label = all_of(label))
  
  data_p <- data_peptides %>% 
    mutate(yvalue = !!rlang::sym(yvalue), 
           Start = !!rlang::sym(Start), 
           End = !!rlang::sym(End), 
           color = !!rlang::sym(color), 
           alpha = !!rlang::sym(alpha), 
           linewidth = !!rlang::sym(linewidth), 
           label = !!rlang::sym(label))
  
  # add if
  if (!hasArg(protein_range))
    protein_range <- c(min(data_p$Start), 
                       max(data_p$End))
  
  # new plot if nothing givenb 
  if (!hasArg(p))
    p_peptides <- ggplot() + 
      annotate("rect", 
               xmin = protein_range[1], 
               xmax = protein_range[2], 
               ymin = - protein_width / 2, 
               ymax =   protein_width / 2, 
               fill = "grey")
  else 
    p_peptides <- p
  
  # Add peptides 
  p_peptides <- p_peptides + 
    geom_segment(aes(x = Start, 
                     xend = End, 
                     y = yvalue, 
                     yend = yvalue, 
                     color = color, 
                     alpha = alpha, 
                     linewidth = linewidth), 
                 data = data_p) + 
    theme_classic(base_size = base_size) + 
    scale_color_manual(values = color.scale) + 
    scale_alpha_manual(values = c(regular = 0.8, 
                                  highlight = 1)) + 
    scale_linewidth_manual(values = c(regular = 1.2, 
                                      highlight = 2.2)) + 
    scale_x_continuous(expand = c(0, 0)) + 
    coord_cartesian(xlim = protein_range, 
                    ylim = c(min(c(-min_y_range, data_p$yvalue), na.rm = T), 
                             max(c(min_y_range, data_p$yvalue), na.rm = T))) + 
    labs(color = color) + 
    guides(linewidth = "none", 
           alpha = "none")
  
  # Annotate peptides 
  if (add.labels)
    p_peptides <- p_peptides + 
    ggrepel::geom_text_repel(ggplot2::aes(x = End, 
                                          y = yvalue, 
                                          label = label), 
                             data = data_p, 
                             size = label.size, 
                             nudge_x = nudge_x,
                             nudge_y = nudge_y,
                             hjust = hjust,
                             vjust = vjust, 
                             direction = direction, 
                             min.segment.length = min.segment.length)
  
  return(p_peptides)
  
}



NGLVieweR_AFdb <- function (id) {
  
  file <- tempfile(fileext = ".pdb")
  
  ok <- tryCatch(
    download.file(paste0("https://alphafold.ebi.ac.uk/files/AF-", 
                         id, "-F1-model_v4.pdb"), destfile = file), 
    error = function(e) 1)
  
  if (!ok) NGLVieweR(file)
  else NA
}



NGLVieweR_pdb <- function(id) {
  
  file <- tempfile(fileext = ".pdb1")
  
  ok <- tryCatch(
    download.file(paste0("https://files.rcsb.org/download/", 
                         id, ".pdb1"), destfile = file), 
    error = function(e) 1)
  
  if (!ok) NGLVieweR(file)
  else NA
}



library(httr)
library(jsonlite)

# Function to retrieve PDB title by ID
get_pdb_title <- function(pdb_id) {
  
  if (length(pdb_id) > 1) {
    
    pdb_id_query <- unique(pdb_id)
    
    pdb_title <- map_chr(setNames(pdb_id_query, pdb_id_query), 
                         get_pdb_title)
    
  } else {
    
    # Construct the API URL for the PDB entry
    url <- paste0("https://data.rcsb.org/rest/v1/core/entry/", toupper(pdb_id))
    
    # Make the GET request
    response <- GET(url)
    
    # Check if the request was successful
    if (status_code(response) == 200) {
      # Parse JSON response
      data <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
      
      # Extract the title
      title <- data$struct$title
      
      # Return a list with PDB ID and title
      return(setNames(title, pdb_id))
    } else {
      # Return error message if request fails
      return(setNames(NA, pdb_id))
    }
    
  }
  
  return(pdb_title[pdb_id])
  
}



dynamic_mode_code <- 
  '
    // Debounce function (now active)
    function debounce(func, wait) {
      var timeout;
      return function() {
        var context = this, args = arguments;
        clearTimeout(timeout);
        timeout = setTimeout(function() {
          timeout = null;
          func.apply(context, args);
        }, wait);
      };
    }

    // Keydown handler
    var keyHandler = function(event) {
      var key = event.key;
      if (key === "ArrowLeft" || key === "ArrowRight"|| key === "ArrowUp"|| key === "ArrowDown") {
        Shiny.setInputValue("key_pressed_arrows", {key: key, time: Date.now()});
      }
    };

    // Debounced handler (active)
    var debouncedKeyHandler = debounce(keyHandler, 1000);  // 1000ms debounce delay

    // Add the debounced event listener (active)
    document.addEventListener("keydown", debouncedKeyHandler);

    // Non-debounced listener (commented out)
    // document.addEventListener("keydown", keyHandler);
  '

iframe_code <- "
    .recalculating {opacity: 1.0;}
    .navbar.navbar-static-top {
      background-color: #18974C !important;
    }
    .bslib-page-title.navbar-brand {
      color: white !important;
      background-color: #18974C !important;
    }
    /* Style for iframe card */
    .iframe-card {
      width: 100% !important;
      height: 100% !important;
      margin: 0 !important;
      padding: 0 !important;
    }
    /* Make iframe fill card body */
    .iframe-card .card-body {
      padding: 0 !important; /* Remove card padding */
      width: 100% !important;
      height: 100% !important;
    }
    /* Iframe fills container */
    .iframe-card iframe {
      width: 100% !important;
      height: 100% !important;
      border: none !important;
      display: block;
    }
    /* Tooltip box size and appearance */
    .tooltip-inner {
      max-width: 400px;   /* default is 200px */
      font-size: 16px;
      padding: 10px 15px;
    }

    /* Popover box size and appearance */
    .popover {
      max-width: 750px;   /* default is ~276px */
    }

    .popover-body {
      font-size: 16px;
      padding: 15px;
    }
  "



dataset_description <-
  "1. Staurosporine K562:<br>

This dataset is about using HT-PELSA to screen staurosporine targets in the K562 lysates with staurosporine concentration of 20 uM, related to Fig.1B and S1
<br><br>

2. Staurosporine K562 dose-response:<br>

This dataset is about using HT-PELSA to screen staurosporine targets and binding affinties in the K562 lysates with a series of staurosporine concentrations: 0.2 nM, 2 nM, 20 nM, 200 nM, 2 μM, 20 μM, and 100 μM, related to Fig.1C-1G,and S2.
<br><br>

3. ATP E.coli:<br>

This dataset is about using HT-PELSA to screen ATP targets in the Ecoli lysates with a series of ATP concentrations:

0 μM, 100 μM, 200 μM, 500 μM, 1 mM, 5 mM, and 10 mM, related to Fig.2 and Fig.S3 and S4.
<br><br>

4. ATP E.coli separately:<br>

This dataset is about using HT-PELSA to screen ATP targets in the Ecoli lysates with a series of ATP concentrations:

0. μM, 100 μM, 200 μM, 500 μM, 1 mM, 5 mM, and 10 mM, but by aliquoting lysates into different lysates before ATP treatment, related to Fig.S3A and S3B.
<br><br>

5. Dasatinib K562 Supernatant:<br>

This dataset involves using HT-PELSA to screen Dasatinib targets in cleared K562 lysates, prepared by centrifuging crude lysates at 20,000g for 10 minutes and collecting the supernatants, with a Dasatinib concentration of 5 μM, related to Fig.3B, S5A-S5C.
<br><br>

6. Dasatinib K562 Crude:<br>

This dataset involves using HT-PELSA to screen Dasatinib targets in crude K562 lysates, with a Dasatinib concentration of 5 μM, related to Fig.3C, S5A-S5F.
<br><br>

7. Ampicillin E.coli:<br>

This dataset involves using HT-PELSA to screen Ampicillin targets in crude E.coli lysates, with a Ampicillin concentration of 10 μM, related to Fig.3D,3E,S5G,and S5H.
<br><br>

8. Staurosporine Liver:<br>

This dataset involves using HT-PELSA to screen staurosporine targets in crude mouse liver lysates, with a staurosporine concentration of 20 μM, related to Fig.3F and 3G.
<br><br>

9. Staurosporine Liver Separately:<br>

This dataset involves using HT-PELSA to screen staurosporine targets in crude mouse liver lysates, with a staurosprine concentration of 20 μM, but by aliquoting lysates into different lysates before staurosporine treatment, related to Fig.3G and S6A.
<br><br>

10. Sunitinib Heart:<br>

This dataset involves using HT-PELSA to screen sunitinib targets in crude mouse heart lysates, with a series of sunitinib concentrations: 0.2 nM, 2 nM, 20 nM, 200 nM, 1 μM, 5 μM, and 20 μM, related to Fig.3H.
"

