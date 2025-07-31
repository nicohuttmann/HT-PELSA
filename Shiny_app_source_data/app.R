

library(shiny)
library(bslib)
library(bsicons)
library(tidyverse)
library(DT)
library(plotly)
library(crosstalk)
library(NGLVieweR)
library(drc)
library(ggtext)
library(pOmics3)

library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::first)

source("PELSA-Package-App/functions.R")


# Load data 
input_list <- readRDS("htpelsaapp_final/htpelsaapp_data.Rds")



# ---- UI ---- 
ui <- page_sidebar(
  theme = bs_theme(
    version = 5, 
    bootswatch = "flatly", 
    primary = "#18974C", 
    secondary = "#0A5032", 
    "bs-link-color" = "#18974C", 
    "bs-nav-link-color" = "#18974C", 
    "bs-link-hover-color" = "#18974C", 
    "bs-pagination-bg" = "#18974C", 
    "bs-pagination-hover-bg" = "#18974C", 
    "bs-pagination-disabled-bg" = "#18974C"
  ), 
  # Remove graying out of plots when loading and set colors 
  tags$style(type = "text/css", iframe_code),
  #   tags$style(HTML("
  #   .tooltip {
  #     border: 1px solid black;
  #   }
  #   
  #   .tooltip-inner {
  #     font-size: 14px;
  #     color: #222;
  #     background-color: #f9f9f9;
  #     font-family: 'Courier New', monospace;
  #   }
  # 
  #   .popover-body {
  #     font-size: 16px;
  #     color: #222;
  #     background-color: #f9f9f9;
  #     font-family: Arial, sans-serif;
  #   }
  # 
  #   .popover-header {
  #     font-weight: bold;
  #     background-color: #e9ecef;
  #     color: #333;
  #   }
  # ")), 
  # Add key control of the app in dynamic mode 
  tags$script(HTML(dynamic_mode_code)), 
  
  # App title 
  title = "HT-PELSA APP", 
  # Sidebar content 
  sidebar = sidebar(width = 300, 
                    bg = "#0000000d", 
                    fg = "#000000", 
                    tags$h4("Dataset"), 
                    selectizeInput("select_dataset", 
                                   "Select: ", 
                                   choices = NULL), 
                    hr(), 
                    tags$h4("Protein"), 
                    input_switch("filter_protein_selection", 
                                 tooltip(
                                   trigger = list("Show hits only?", 
                                                  bs_icon("info-circle")), 
                                   HTML("Should only stabilized and destabilized proteins be listed in the selection below or all protein identified?<br><br>
                                          If this option is enabled, non-hit proteins selected in the volcano plot won't be shown in the selection list below unless you disabled.")), 
                                 TRUE),  
                    input_switch("protein_selection_stabilized", 
                                 tooltip(
                                   trigger = list("Stabilized first?", 
                                                  bs_icon("info-circle")), 
                                   "The list of proteins is ordered by the log2 fold-change or the signed pEC50. Stabilized proteins are listed first by default"), 
                                 TRUE),
                    # Choose protein via gene, UniProt id or name 
                    selectizeInput("select_gene", 
                                   "Select by gene:", 
                                   choices = NULL,), 
                    selectizeInput("select_protein", 
                                   "Select by UniProt ID:", 
                                   choices = NULL), 
                    selectizeInput("select_protein_name", 
                                   "Select by protein name:", 
                                   choices = NULL), 
                    # Settings 
                    hr(), 
                    tags$h4("Settings"), 
                    # Further customize the app 
                    input_switch("collapse_proteins", 
                                 tooltip(
                                   trigger = list("Plots at protein-level?", 
                                                  bs_icon("info-circle")), 
                                   "Proteins are represented by the peptide with the lowest p-value or highest absolute pEC50."), 
                                 TRUE),  
                    # Dynamic mode to control the app via arrow keys 
                    input_switch("dynamic_mode", 
                                 tooltip(
                                   trigger = list("Dynamic Mode", 
                                                  bs_icon("info-circle")), 
                                   HTML("ArrowRight = Next protein<br>
                            ArrowLeft = Previous protein<br>
                            ArrowDown = Next peptide<br>
                            ArrowUP = Previous peptide")), 
                                 TRUE), 
                    # Enable explorer mode for more functions 
                    input_switch("explorer_mode", "Explorer Mode", FALSE), 
                    hr(), 
                    tags$h4("About"), 
                    popover(
                      actionLink("link", "Description"),
                      HTML(dataset_description),
                      placement = "right"),
                    popover(
                      actionLink("link", "Contact"),
                      "For questions, please contact Nico Hüttmann (nico.huettmann@embl.de) or Kejia Li (kejia.li@embl.de).", 
                      placement = "right"), 
                    popover(
                      actionLink("link", "Publication"),
                      "This app accompanies the HT-PELSA ", 
                      tags$a("publication",
                             target = "_blank",
                             href = "https://www.biorxiv.org/content/10.1101/2025.04.28.650974v1"), 
                      ".", 
                      placement = "right"), 
                    popover(
                      actionLink("link", "Source data"),
                      "You can find the source data for this app and more ", 
                      tags$a("here",
                             target = "_blank",
                             href = "https://github.com/nicohuttmann/HT-PELSA/"), 
                      ".", 
                      placement = "right")
  ), 
  # Main area content (will be rendered in the server due to explorer mode)
  uiOutput("dynamic_layout", fill = T), 
)



# ---- Server ----
server <- function(input, output, session) {
  
  # ---- Render main area layout ----
  output$dynamic_layout <- renderUI({
    # Depending on explorer mode enabled or disabled 
    if (!input$explorer_mode) {
      # Layout of default app (2 columns of equal width)
      layout_columns(
        # Layout parameters 
        col_widths = c(6, 6, 12), 
        row_heights = c(3, 3, 2), 
        
        # Left column: pEC50 overview plot, table and peptide curves 
        card(NGLVieweROutput("protein_structure"), 
             #actionButton("snapshot", "Snapshot"), 
             full_screen = T), 
        
        # Right column: Protein structure and sequence plot 
        card(plotOutput("overview_plot", 
                        click = "overview_plot_click", 
                        brush = brushOpts(id = "overview_plot_brush", 
                                          fill = "#18974C", 
                                          stroke = "#0A5032", 
                                          opacity = 0.25,
                                          delay = 500, 
                                          delayType = c("debounce", 
                                                        "throttle"), 
                                          clip = FALSE,
                                          direction = c("xy"), 
                                          resetOnNew = TRUE), 
                        dblclick = "overview_plot_dbclick"), 
             full_screen = T), 
        card(
          plotOutput("protein_sequence", 
                     brush = brushOpts(id = "protein_sequence_plot_brush", 
                                       fill = "#18974C", 
                                       stroke = "#0A5032", 
                                       opacity = 0.25,
                                       delay = 500, 
                                       delayType = c("debounce", 
                                                     "throttle"), 
                                       clip = FALSE,
                                       direction = c("x"), 
                                       resetOnNew = TRUE), 
                     dblclick = "protein_sequence_plot_dbclick"), 
          full_screen = T))
    } else {
      # Explorer mode 
      layout_columns(
        # Layout parameters 
        col_widths = c(7, 5), 
        row_heights = c(1), 
        
        layout_columns(
          # Layout parameters 
          col_widths = c(6, 6, 12), 
          row_heights = c(3, 3, 2), 
          
          # Left column: pEC50 overview plot, table and peptide curves 
          card(NGLVieweROutput("protein_structure"), 
               #actionButton("snapshot", "Snapshot"), 
               full_screen = T), 
          
          # Right column: Protein structure and sequence plot 
          card(plotOutput("overview_plot", 
                          click = "overview_plot_click", 
                          brush = brushOpts(id = "overview_plot_brush", 
                                            fill = "#18974C", 
                                            stroke = "#0A5032", 
                                            opacity = 0.25,
                                            delay = 500, 
                                            delayType = c("debounce", 
                                                          "throttle"), 
                                            clip = FALSE,
                                            direction = c("xy"), 
                                            resetOnNew = TRUE), 
                          dblclick = "overview_plot_dbclick"), 
               full_screen = T), 
          card(
            plotOutput("protein_sequence", 
                       brush = brushOpts(id = "protein_sequence_plot_brush", 
                                         fill = "#18974C", 
                                         stroke = "#0A5032", 
                                         opacity = 0.25,
                                         delay = 500, 
                                         delayType = c("debounce", 
                                                       "throttle"), 
                                         clip = FALSE,
                                         direction = c("x"), 
                                         resetOnNew = TRUE), 
                       dblclick = "protein_sequence_plot_dbclick"), 
            full_screen = T)), 
        
        # Explorer mode area 
        navset_card_tab(
          full_screen = TRUE,
          title = "Explorer Mode", 
          id = "expmode_navset", 
          # Table 
          nav_panel(
            "Table", 
            value = "table", 
            DTOutput("results_table")
          ), 
          # Raw data 
          nav_panel(
            "Peptides", 
            value = "raw", 
            selectizeInput("select_peptide", 
                           "Select peptide:", 
                           choices = NULL), 
            plotOutput("peptide_ds_curves"), 
            tableOutput("peptide_EC50s")
          ),
          # UniProt info tab
          # nav_panel(
          #   "Info", 
          #   value = "info", 
          #   "UniProt info - more to come"
          # ),
          # PDB structure tab
          nav_panel(
            "PDB",
            value = "pdb", 
            selectizeInput("expmode_pdb_id", 
                           "Select PDB entry:", 
                           choices = NULL), 
            textOutput("expmode_pdb_title"), 
            NGLVieweROutput(outputId = "expmode_pdb_structure")
          )#,
          # nav_panel(
          #   "GO",
          #   value = "go", 
          #   "GO content - more to come"
          # )
        )
      )
    }
  })
  
  
  
  # ---- Main data source ---- 
  
  ## Select dataset ----
  updateSelectizeInput(session, 
                       "select_dataset", 
                       choices = names(input_list), 
                       selected = names(input_list)[1], 
                       server = TRUE)
  
  dataset_selected <- reactiveVal(names(input_list)[1])
  
  data_list <- reactive(input_list[[dataset_selected()]])
  
  
  
  # Select dataset and update selections 
  observeEvent(input$select_dataset, {
    if (!is.null(input$select_dataset) && 
        input$select_dataset != "" &&
        input$select_dataset != dataset_selected()) {
      
      
      # Change protein if not present in next dataset/choices
      new_dataset_genes <- 
        input_list[[input$select_dataset]]$data_peptides %>% 
        filter(if (input$filter_protein_selection) regulation %in% c("stabilized", "destabilized") else T) %>% 
        arrange(
          (if ("log2.fc" %in% names(.)) 1 else -1)
          *
            (if (input$protein_selection_stabilized) 1 else -1)
          *
            (if("log2.fc" %in% names(.)) log2.fc else pEC50_signed)
        ) %>% 
        filter(!duplicated(Protein.Group)) %>%  
        pull(Protein.Group, Genes)
      
      new_dataset_protein.groups <- 
        input_list[[input$select_dataset]]$data_peptides %>% 
        filter(if (input$filter_protein_selection) regulation %in% c("stabilized", "destabilized") else T) %>% 
        arrange(
          (if ("log2.fc" %in% names(.)) 1 else -1)
          *
            (if (input$protein_selection_stabilized) 1 else -1)
          *
            (if("log2.fc" %in% names(.)) log2.fc else pEC50_signed)
        ) %>% 
        filter(!duplicated(Protein.Group)) %>% 
        pull(Protein.Group, Protein.Group)
      
      
      new_dataset_protein.description <- 
        input_list[[input$select_dataset]]$data_peptides %>% 
        filter(if (input$filter_protein_selection) regulation %in% c("stabilized", "destabilized") else T) %>% 
        arrange(
          (if ("log2.fc" %in% names(.)) 1 else -1)
          *
            (if (input$protein_selection_stabilized) 1 else -1)
          *
            (if("log2.fc" %in% names(.)) log2.fc else pEC50_signed)
        ) %>% 
        filter(!duplicated(Protein.Group)) %>% 
        pull(Protein.Group, Protein.Description)
      
      if (protein_group_selected() %in% new_dataset_genes) {
        new_protein_group_selected <- protein_group_selected()
      } else {
        new_protein_group_selected <- unname(new_dataset_genes)[1]
      }
      
      
      
      dataset_selected(input$select_dataset)
      
      protein_group_selected(new_protein_group_selected)
      
      updateSelectizeInput(
        session,
        "select_gene",
        choices = NULL,
        selected = NULL,
        server = TRUE
      )
      updateSelectizeInput(session,
                           "select_gene",
                           choices = new_dataset_genes, 
                           selected = new_protein_group_selected)
      
      updateSelectizeInput(
        session,
        "select_protein",
        choices = NULL,
        selected = NULL,
        server = TRUE
      )
      updateSelectizeInput(session,
                           "select_protein",
                           choices = new_dataset_protein.groups, 
                           selected = new_protein_group_selected)
      
      updateSelectizeInput(
        session,
        "select_protein_name",
        choices = NULL,
        selected = NULL,
        server = TRUE
      )
      updateSelectizeInput(session,
                           "select_protein_name", 
                           choices = new_dataset_protein.description, 
                           selected = new_protein_group_selected)
      
      
      # Update selected peptide
      updateSelectizeInput(
        session,
        "select_peptide",
        choices = NULL,
        selected = NULL,
        server = TRUE
      )
      
      peptide_choices <- input_list[[input$select_dataset]]$data_peptides %>%
        filter(Protein.Group == new_protein_group_selected) %>%
        arrange(from, to) %>%
        mutate(
          peptide_display = paste0(
            Peptide.Id, " (", from, "-", to, ")"
          )
        ) %>%
        pull(Peptide.Id, name = peptide_display)
      
      selected_peptide <- if (length(peptide_choices) > 0) {
        peptide_choices[[1]]  # Select the first peptide
      } else {
        NULL
      }
      peptide_selected(selected_peptide)
      updateSelectizeInput(
        session,
        "select_peptide",
        choices = peptide_choices,
        selected = selected_peptide,
        server = TRUE)
    }
    
  }, ignoreInit = TRUE)
  
  
  
  # Reload choices after filter/sort option
  observeEvent(c(input$filter_protein_selection, 
                 input$protein_selection_stabilized), {
                   new_dataset_genes <- 
                     input_list[[input$select_dataset]]$data_peptides %>% 
                     filter(if (input$filter_protein_selection) regulation %in% c("stabilized", "destabilized") else T) %>% 
                     arrange(
                       (if ("log2.fc" %in% names(.)) 1 else -1)
                       *
                         (if (input$protein_selection_stabilized) 1 else -1)
                       *
                         (if("log2.fc" %in% names(.)) log2.fc else pEC50_signed)
                     ) %>% 
                     filter(!duplicated(Protein.Group)) %>%  
                     pull(Protein.Group, Genes)
                   
                   new_dataset_protein.groups <- 
                     input_list[[input$select_dataset]]$data_peptides %>% 
                     filter(if (input$filter_protein_selection) regulation %in% c("stabilized", "destabilized") else T) %>% 
                     arrange(
                       (if ("log2.fc" %in% names(.)) 1 else -1)
                       *
                         (if (input$protein_selection_stabilized) 1 else -1)
                       *
                         (if("log2.fc" %in% names(.)) log2.fc else pEC50_signed)
                     ) %>% 
                     filter(!duplicated(Protein.Group)) %>% 
                     pull(Protein.Group, Protein.Group)
                   
                   
                   new_dataset_protein.description <- 
                     input_list[[input$select_dataset]]$data_peptides %>% 
                     filter(if (input$filter_protein_selection) regulation %in% c("stabilized", "destabilized") else T) %>% 
                     arrange(
                       (if ("log2.fc" %in% names(.)) 1 else -1)
                       *
                         (if (input$protein_selection_stabilized) 1 else -1)
                       *
                         (if("log2.fc" %in% names(.)) log2.fc else pEC50_signed)
                     ) %>% 
                     filter(!duplicated(Protein.Group)) %>% 
                     pull(Protein.Group, Protein.Description)
                   
                   if (protein_group_selected() %in% new_dataset_genes) {
                     new_protein_group_selected <- protein_group_selected()
                   } else {
                     new_protein_group_selected <- unname(new_dataset_genes)[1]
                     protein_group_selected(new_protein_group_selected)
                   }
                   
                   
                   updateSelectizeInput(
                     session,
                     "select_gene",
                     choices = NULL,
                     selected = NULL,
                     server = TRUE
                   )
                   updateSelectizeInput(session,
                                        "select_gene",
                                        choices = new_dataset_genes, 
                                        selected = new_protein_group_selected)
                   
                   updateSelectizeInput(
                     session,
                     "select_protein",
                     choices = NULL,
                     selected = NULL,
                     server = TRUE
                   )
                   updateSelectizeInput(session,
                                        "select_protein",
                                        choices = new_dataset_protein.groups, 
                                        selected = new_protein_group_selected)
                   
                   updateSelectizeInput(
                     session,
                     "select_protein_name",
                     choices = NULL,
                     selected = NULL,
                     server = TRUE
                   )
                   updateSelectizeInput(session,
                                        "select_protein_name", 
                                        choices = new_dataset_protein.description, 
                                        selected = new_protein_group_selected)
                   
                   
                   # Update selected peptide
                   updateSelectizeInput(
                     session,
                     "select_peptide",
                     choices = NULL,
                     selected = NULL,
                     server = TRUE
                   )
                   
                   peptide_choices <- input_list[[input$select_dataset]]$data_peptides %>%
                     filter(Protein.Group == new_protein_group_selected) %>%
                     arrange(from, to) %>%
                     mutate(
                       peptide_display = paste0(
                         Peptide.Id, " (", from, "-", to, ")"
                       )
                     ) %>%
                     pull(Peptide.Id, name = peptide_display)
                   
                   selected_peptide <- if (length(peptide_choices) > 0) {
                     peptide_choices[[1]]  # Select the first peptide
                   } else {
                     NULL
                   }
                   peptide_selected(selected_peptide)
                   updateSelectizeInput(
                     session,
                     "select_peptide",
                     choices = peptide_choices,
                     selected = selected_peptide,
                     server = TRUE)
                   
                 }, 
               ignoreInit = TRUE)
  
  
  # Define colors for peptides in all plots 
  peptide_colors <- c(underquantified = "#b0b0b0",  #"grey40", 
                      
                      unchanged = "#666666", #"grey30",
                      unresponsive = "#666666", #"grey30", 
                      none = "#666666", #"grey30", #4d4d4d
                      
                      regular = "#666666", 
                      unknown = "#734595", #"grey30",
                      
                      destabilized = "#FF0000", # "red", 
                      stabilized = "#0000FF", # "blue", 
                      highlight = "#18974C")  
  
  
  
  ## ---- Select Protein to show ---- 
  
  # Reactive protein name
  protein_default <- input_list[[1]]$data_peptides %>% 
    filter(regulation %in% c("stabilized", "destabilized")) %>% 
    filter(!duplicated(Protein.Group)) %>% 
    mutate(protein_id = paste0(Genes, " (", Protein.Group, ")")) %>% 
    arrange(if("log2.fc" %in% names(.)) log2.fc else -pEC50_signed) %>% 
    pull(Protein.Group) %>% 
    first()
  
  protein_group_selected <- reactiveVal(protein_default)
  
  peptide_selected <- reactiveVal(input_list[[1]]$data_peptides %>% 
                                    filter(Protein.Group == 
                                             protein_default) %>% 
                                    arrange(from, to) %>% 
                                    pull(Peptide.Id) %>% 
                                    first())
  
  
  # ---- Initialize selections ---- 
  
  # Select by gene 
  updateSelectizeInput(session, 
                       "select_gene", 
                       choices = input_list[[1]]$data_peptides %>% 
                         filter(regulation %in% c("stabilized", "destabilized")) %>% 
                         arrange(if("log2.fc" %in% names(.)) log2.fc else -pEC50_signed) %>% 
                         filter(!duplicated(Protein.Group)) %>% 
                         pull(Protein.Group, Genes), 
                       selected = protein_default, 
                       server = TRUE)
  
  # Select by UniProt ID 
  updateSelectizeInput(session, 
                       "select_protein", 
                       choices = input_list[[1]]$data_peptides %>% 
                         filter(regulation %in% c("stabilized", "destabilized")) %>% 
                         arrange(if("log2.fc" %in% names(.)) log2.fc else -pEC50_signed) %>% 
                         filter(!duplicated(Protein.Group)) %>% 
                         pull(Protein.Group, Protein.Group), 
                       selected = protein_default, 
                       server = TRUE)
  
  # Select by UniProt ID 
  updateSelectizeInput(session, 
                       "select_protein_name", 
                       choices = input_list[[1]]$data_peptides %>% 
                         filter(regulation %in% c("stabilized", "destabilized")) %>% 
                         arrange(if("log2.fc" %in% names(.)) log2.fc else -pEC50_signed) %>% 
                         filter(!duplicated(Protein.Group)) %>% 
                         pull(Protein.Group, Protein.Description), 
                       selected = protein_default, 
                       server = TRUE)
  
  
  # Select peptide
  updateSelectizeInput(session, 
                       "select_peptide", 
                       choices = input_list[[1]]$data_peptides %>%
                         filter(Protein.Group == protein_default) %>%
                         arrange(from, to) %>%
                         mutate(
                           peptide_display = paste0(
                             Peptide.Id, " (", from, "-", to, ")"
                           )
                         ) %>%
                         pull(Peptide.Id, name = peptide_display), 
                       selected = input_list[[1]]$data_peptides %>%
                         filter(Protein.Group == protein_default) %>%
                         arrange(from, to) %>%
                         mutate(
                           peptide_display = paste0(
                             Peptide.Id, " (", from, "-", to, ")"
                           )
                         ) %>%
                         pull(Peptide.Id, name = peptide_display) %>% 
                         first(), 
                       server = TRUE)
  
  
  
  
  
  # ---- Observe input changes ---- 
  
  # By gene
  observeEvent(input$select_gene, {
    if (!is.null(input$select_gene) && 
        input$select_gene != "" &&
        input$select_gene != protein_group_selected()) {
      
      protein_group_selected(input$select_gene)
      updateSelectizeInput(session,
                           "select_protein",
                           selected = input$select_gene)
      
      updateSelectizeInput(session,
                           "select_protein_name",
                           selected = input$select_gene)
      
      # Update selected peptide
      updateSelectizeInput(
        session,
        "select_peptide",
        choices = NULL,
        selected = NULL,
        server = TRUE
      )
      
      peptide_choices <- data_list()$data_peptides %>%
        filter(Protein.Group == input$select_gene) %>%
        arrange(from, to) %>%
        mutate(
          peptide_display = paste0(
            Peptide.Id, " (", from, "-", to, ")"
          )
        ) %>%
        pull(Peptide.Id, name = peptide_display)
      
      selected_peptide <- if (length(peptide_choices) > 0) {
        peptide_choices[[1]]  # Select the first peptide
      } else {
        NULL
      }
      peptide_selected(selected_peptide)
      updateSelectizeInput(
        session,
        "select_peptide",
        choices = peptide_choices,
        selected = selected_peptide,
        server = TRUE)
    }
    
  }, ignoreInit = TRUE)
  
  
  # By protein id
  observeEvent(input$select_protein, {
    if (!is.null(input$select_protein) && 
        input$select_protein != "" &&
        input$select_protein != protein_group_selected()) {
      protein_group_selected(input$select_protein)
      updateSelectizeInput(session,
                           "select_gene",
                           selected = input$select_protein)
      
      updateSelectizeInput(session,
                           "select_protein_name",
                           selected = input$select_protein)
      
      # Update selected peptide
      updateSelectizeInput(
        session,
        "select_peptide",
        choices = NULL,
        selected = NULL,
        server = TRUE
      )
      
      peptide_choices <- data_list()$data_peptides %>%
        filter(Protein.Group == input$select_protein) %>%
        arrange(from, to) %>%
        mutate(
          peptide_display = paste0(
            Peptide.Id, " (", from, "-", to, ")"
          )
        ) %>%
        pull(Peptide.Id, name = peptide_display)
      
      selected_peptide <- if (length(peptide_choices) > 0) {
        peptide_choices[[1]]  # Select the first peptide
      } else {
        NULL
      }
      peptide_selected(selected_peptide)
      updateSelectizeInput(
        session,
        "select_peptide",
        choices = peptide_choices,
        selected = selected_peptide,
        server = TRUE)
    }
    
  }, ignoreInit = TRUE)
  
  
  # By protein name
  observeEvent(input$select_protein_name, {
    if (!is.null(input$select_protein_name) && 
        input$select_protein_name != "" &&
        input$select_protein_name != protein_group_selected()) {
      protein_group_selected(input$select_protein_name)
      updateSelectizeInput(session,
                           "select_gene",
                           selected = input$select_protein_name)
      
      updateSelectizeInput(session,
                           "select_protein",
                           selected = input$select_protein_name)
      
      # Update selected peptide
      updateSelectizeInput(
        session,
        "select_peptide",
        choices = NULL,
        selected = NULL,
        server = TRUE)
      
      peptide_choices <- data_list()$data_peptides %>%
        filter(Protein.Group == input$select_protein_name) %>%
        arrange(from, to) %>%
        mutate(
          peptide_display = paste0(
            Peptide.Id, " (", from, "-", to, ")")) %>%
        pull(Peptide.Id, name = peptide_display)
      
      selected_peptide <- if (length(peptide_choices) > 0) {
        peptide_choices[[1]]  # Select the first peptide
      } else {
        NULL
      }
      peptide_selected(selected_peptide)
      updateSelectizeInput(
        session,
        "select_peptide",
        choices = peptide_choices,
        selected = selected_peptide,
        server = TRUE)
    }
    
  }, ignoreInit = TRUE)
  
  
  
  # By keystroke
  observeEvent(input$key_pressed_arrows, {
    
    if (input$dynamic_mode) {
      
      key <- input$key_pressed_arrows$key
      
      all_proteins <- data_list()$data_peptides %>% 
        arrange(p.value) %>% 
        filter(!duplicated(Protein.Group)) %>%  
        filter(if (input$filter_protein_selection) 
          regulation %in% c("stabilized", "destabilized") 
          else 
            T) %>% 
        arrange(
          (if ("log2.fc" %in% names(.)) 1 else -1)
          *
            (if (input$protein_selection_stabilized) 1 else -1)
          *
            (if("log2.fc" %in% names(.)) log2.fc else pEC50_signed)
        ) %>% 
        pull(Protein.Group)
      
      current_protein_pos <- match(protein_group_selected(), all_proteins)
      
      
      all_peptides <- data_list()$data_peptides %>%
        filter(Protein.Group == protein_group_selected()) %>%
        arrange(from, to) %>%
        mutate(
          peptide_display = paste0(
            Peptide.Id, " (", from, "-", to, ")"
          )
        ) %>%
        pull(Peptide.Id, name = peptide_display)
      
      current_peptide_pos <- match(peptide_selected(), all_peptides)
      
      
      # ArrowLeft: Previous protein
      if (key == "ArrowLeft") {
        
        if (current_protein_pos != 1) 
          new_protein <- all_proteins[current_protein_pos - 1]
        else 
          new_protein <- all_proteins[length(all_proteins)]
        # Update protein group and selections
        protein_group_selected(new_protein)
        updateSelectizeInput(session,
                             "select_gene",
                             selected = new_protein)
        
        updateSelectizeInput(session,
                             "select_protein",
                             selected = new_protein)
        
        updateSelectizeInput(session,
                             "select_protein_name",
                             selected = new_protein)
        
        # Update selected peptide
        updateSelectizeInput(
          session,
          "select_peptide",
          choices = NULL,
          selected = NULL,
          server = TRUE
        )
        peptide_choices <- data_list()$data_peptides %>%
          filter(Protein.Group == new_protein) %>%
          arrange(from, to) %>%
          mutate(
            peptide_display = paste0(
              Peptide.Id, " (", from, "-", to, ")"
            )
          ) %>%
          pull(Peptide.Id, name = peptide_display)
        
        selected_peptide <- if (length(peptide_choices) > 0) {
          peptide_choices[[1]]  # Select the first peptide
        } else {
          NULL
        }
        peptide_selected(selected_peptide)
        updateSelectizeInput(
          session,
          "select_peptide",
          choices = peptide_choices,
          selected = selected_peptide,
          server = TRUE)
        
        
        # ArrowRight: Next protein
      } else if (key == "ArrowRight") {
        
        if (current_protein_pos != length(all_proteins)) 
          new_protein <- all_proteins[current_protein_pos + 1]
        else 
          new_protein <- all_proteins[1]
        
        # Update protein group and selections
        protein_group_selected(new_protein)
        updateSelectizeInput(session,
                             "select_gene",
                             selected = new_protein)
        
        updateSelectizeInput(session,
                             "select_protein",
                             selected = new_protein)
        
        updateSelectizeInput(session,
                             "select_protein_name",
                             selected = new_protein)
        
        # Update selected peptide
        updateSelectizeInput(
          session,
          "select_peptide",
          choices = NULL,
          selected = NULL,
          server = TRUE
        )
        peptide_choices <- data_list()$data_peptides %>%
          filter(Protein.Group == new_protein) %>%
          arrange(from, to) %>%
          mutate(
            peptide_display = paste0(
              Peptide.Id, " (", from, "-", to, ")"
            )
          ) %>%
          pull(Peptide.Id, name = peptide_display)
        
        selected_peptide <- if (length(peptide_choices) > 0) {
          peptide_choices[[1]]  # Select the first peptide
        } else {
          NULL
        }
        peptide_selected(selected_peptide)
        updateSelectizeInput(
          session,
          "select_peptide",
          choices = peptide_choices,
          selected = selected_peptide,
          server = TRUE)
        
        
        # ArrowUp: Previous peptide
      } else if (key == "ArrowUp") {
        
        if (current_peptide_pos != 1) 
          new_peptide <- all_peptides[current_peptide_pos - 1]
        else 
          new_peptide <- all_peptides[length(all_peptides)]
        
        # updateSelectizeInput(
        #   session,
        #   "select_peptide",
        #   choices = NULL,
        #   selected = NULL,
        #   server = TRUE
        # )
        
        peptide_selected(new_peptide)
        
        updateSelectizeInput(
          session,
          "select_peptide",
          choices = all_peptides,
          selected = new_peptide,
          server = TRUE)
        
        
        # ArrowDown: Next peptide
      } else if (key == "ArrowDown") {
        
        if (current_peptide_pos != length(all_peptides)) 
          new_peptide <- all_peptides[current_peptide_pos + 1]
        else 
          new_peptide <- all_peptides[1]
        # updateSelectizeInput(
        #   session,
        #   "select_peptide",
        #   choices = NULL,
        #   selected = NULL,
        #   server = TRUE
        # )
        
        peptide_selected(new_peptide)
        
        updateSelectizeInput(
          session,
          "select_peptide",
          choices = all_peptides,
          selected = new_peptide,
          server = TRUE)
        
      }
    }
  })
  
  # ---- Output selected protein group (for control) ---- 
  output$protein_group <- renderText(
    paste0(protein_group_selected(), 
           " (# peptides = ", 
           nrow(filter(data_list()$data_peptides, 
                       Protein.Group == protein_group_selected())), 
           ")"))
  
  
  # ---- Select peptide ---- 
  observeEvent(input$select_peptide, {
    req(protein_group_selected(), peptide_selected())
    
    if (input$select_peptide == "") {
      
      peptide_choices <- data_list()$data_peptides %>%
        filter(Protein.Group == input$select_protein_name) %>%
        arrange(from, to) %>%
        mutate(
          peptide_display = paste0(
            Peptide.Id, " (", from, "-", to, ")")) %>%
        pull(Peptide.Id, name = peptide_display)
      
      selected_peptide <- if (length(peptide_choices) > 0) {
        peptide_choices[[1]]  # Select the first peptide
      } else {
        NULL
      }
      
      peptide_selected(selected_peptide)
      updateSelectizeInput(
        session,
        "select_peptide",
        choices = peptide_choices,
        selected = selected_peptide,
        server = TRUE)
      
    }
    
    if (!is.null(input$select_peptide) &&
        input$select_peptide != peptide_selected() && 
        input$select_peptide != "" && 
        input$select_peptide %in% 
        filter(data_list()$data_peptides, Protein.Group == 
               protein_group_selected())$Peptide.Id) {
      peptide_selected(input$select_peptide)
    }})
  
  peptide_selected_debounced <- debounce(reactive(peptide_selected()), 
                                         millis = 700)
  
  ## By clicking on overview plot 
  observeEvent(input$overview_plot_click, {
    
    if ("pEC50" %in% names(table_peptides_plot())) {
      np <- nearPoints(table_peptides_plot(), 
                       input$overview_plot_click, 
                       xvar = "pEC50", 
                       yvar = "CV_pEC50")
    } else {
      np <- nearPoints(table_peptides_plot(), 
                       input$overview_plot_click, 
                       xvar = "log2.fc", 
                       yvar = "log10.p.value")
    }
    
    if(nrow(np) > 0) {
      new_protein <- 
        filter(data_list()$data_peptides, 
               Peptide.Id == np[["Peptide.Id"]][1])$Protein.Group
      
      
      
      if (protein_group_selected() != new_protein) {
        
        protein_group_selected(new_protein)
        
        updateSelectizeInput(session,
                             "select_gene",
                             selected = new_protein)
        
        updateSelectizeInput(session,
                             "select_protein",
                             selected = new_protein)
        
        updateSelectizeInput(session,
                             "select_protein_name",
                             selected = new_protein)
      } 
      
      
      # Change peptide choices
      peptide_choices <- data_list()$data_peptides %>%
        filter(Protein.Group == new_protein) %>%
        arrange(from, to) %>%
        mutate(
          peptide_display = paste0(
            Peptide.Id, " (", from, "-", to, ")"
          )
        ) %>%
        pull(Peptide.Id, name = peptide_display)
      
      peptide_selected(np[["Peptide.Id"]][1])
      
      updateSelectizeInput(
        session,
        "select_peptide",
        choices = NULL,
        selected = NULL,
        server = TRUE)
      
      updateSelectizeInput(
        session,
        "select_peptide",
        choices = peptide_choices,
        selected = np[["Peptide.Id"]][1],
        server = TRUE)
    }
  })
  
  
  
  
  # ---- Plot overview EC50/volcano ---- 
  
  # Define limits 
  overview_plot_limits <- reactiveValues()
  
  observeEvent(input$overview_plot_brush, {
    overview_plot_limits$xmax <- input$overview_plot_brush$xmax
    overview_plot_limits$ymax <- input$overview_plot_brush$ymax
    
    overview_plot_limits$xmin <- input$overview_plot_brush$xmin
    overview_plot_limits$ymin <- input$overview_plot_brush$ymin
  })
  
  observeEvent(c(input$overview_plot_dbclick), {
    overview_plot_limits$xmin <- NULL
    overview_plot_limits$ymin <- NULL
    overview_plot_limits$xmax <- NULL
    overview_plot_limits$ymax <- NULL
  })
  
  # Data for plot 
  table_peptides_plot <- 
    reactive({
      if ("pEC50_signed" %in% names(data_list()$data_peptides)) {
        data_list()$data_peptides %>% 
          mutate(peptide_type = case_when(
            Protein.Group == protein_group_selected() ~ "highlight", 
            .default = "regular")) %>% #.default = regulation)) %>% 
          arrange(desc(pEC50)) %>% 
          filter(if (input$collapse_proteins) !duplicated(Protein.Group) else T) %>% 
          arrange(rank(match(peptide_type, names(peptide_colors)), 
                       na.last = F, ties.method = "first"))
      } else {
        data_list()$data_peptides %>% 
          mutate(peptide_type = case_when(
            Protein.Group == protein_group_selected() ~ "highlight", 
            !regulation_peptide %in% c("stabilized", "destabilized") ~ "regular", 
            .default = "regular"), 
            log10.p.value = -log10(p.value)) %>% 
          arrange(p.value) %>% 
          filter(if (input$collapse_proteins) !duplicated(Protein.Group) else T) %>% 
          arrange(rank(match(peptide_type, names(peptide_colors)), 
                       na.last = F, ties.method = "first"))
      }
    })
  
  # Render plot
  output$overview_plot <- renderPlot({
    
    if ("pEC50_signed" %in% names(table_peptides_plot())) {
      
      table_peptides_plot() %>% 
        filter(regulation %in% c("stabilized", "destabilized")) %>% 
        ggplot(aes(x = pEC50, 
                   y = CV_pEC50, 
                   color = regulation,  
                   alpha = peptide_type, 
                   size = peptide_type)) + 
        geom_point(shape = 16) + 
        theme_classic(15) + 
        theme(
          plot.title = element_textbox_simple()) + 
        scale_color_manual(values = peptide_colors) +
        scale_alpha_manual(values = c(regular = 0.3, 
                                      highlight = 1, 
                                      stabilized = 0.2, 
                                      destabilized = 0.2)) + 
        scale_size_manual(values = c(regular = 2, 
                                     highlight = 3, 
                                     stabilized = 2, 
                                     destabilized = 2)) + 
        coord_cartesian(xlim = c(overview_plot_limits$xmin, 
                                 overview_plot_limits$xmax), 
                        ylim = c(overview_plot_limits$ymin, 
                                 overview_plot_limits$ymax)) + 
        labs(title = if (!input$collapse_proteins) 
          paste0("<span style='color:#000000'>Peptides highlighted for ", 
                 filter(data_list()$data_peptides, Protein.Group == protein_group_selected())$Genes[[1]], 
                 " (<span style='color:#FF0000'>de-/<span style='color:#0000FF'>stabilized<span style='color:#000000'>)</span>")
          else 
            paste0("<span style='color:#000000'>", 
                   filter(data_list()$data_peptides, Protein.Group == protein_group_selected())$Genes[[1]], 
                   " (<span style='color:#FF0000'>de-/<span style='color:#0000FF'>stabilized<span style='color:#000000'>)</span>"), 
          x = expression(paste("mean ", pEC[50])), 
          y = expression(paste("CV ", pEC[50]))) + 
        guides(alpha = "none", 
               size = "none")
      
      # Volcano plot 
    } else {
      table_peptides_plot() %>% 
        ggplot(aes(x = log2.fc, 
                   y = log10.p.value, 
                   color = regulation_peptide,  
                   alpha = peptide_type, 
                   size = peptide_type)) + 
        geom_point(shape = 16) + 
        theme_classic(15) + 
        theme(
          plot.title = element_textbox_simple()) + 
        scale_color_manual(values = peptide_colors) +
        scale_alpha_manual(values = c(regular = 0.3, 
                                      highlight = 1, 
                                      stabilized = 0.2, 
                                      destabilized = 0.2)) + 
        scale_size_manual(values = c(regular = 2, 
                                     highlight = 3, 
                                     stabilized = 2, 
                                     destabilized = 2)) + 
        coord_cartesian(xlim = c(overview_plot_limits$xmin, 
                                 overview_plot_limits$xmax), 
                        ylim = c(overview_plot_limits$ymin, 
                                 overview_plot_limits$ymax)) + 
        labs(title = if (!input$collapse_proteins) 
          paste0("<span style='color:#000000'>Peptides highlighted for ", 
                 filter(data_list()$data_peptides, Protein.Group == protein_group_selected())$Genes[[1]], 
                 " (<span style='color:#FF0000'>de-/<span style='color:#0000FF'>stabilized<span style='color:#000000'>)</span>")
          else 
            paste0("<span style='color:#000000'>", 
                   filter(data_list()$data_peptides, Protein.Group == protein_group_selected())$Genes[[1]], 
                   " (<span style='color:#FF0000'>de-/<span style='color:#0000FF'>stabilized<span style='color:#000000'>)</span>"), 
          x = expression(paste(log[2], " fold-change")), 
          y = expression(paste(-log[10], " p-value"))) + 
        guides(alpha = "none", 
               size = "none")
      
      
    }
  })
  
  
  
  
  
  # ---- Protein table ---- 
  output$results_table <- renderDT({
    if ("pEC50_signed" %in% names(data_list()$data_peptides)) {
      data_list()$data_peptides %>% 
        mutate(pEC50 = signif(pEC50, 3), 
               pEC50_signed = signif(pEC50_signed, 3), 
               CV_pEC50 = signif(CV_pEC50, 3)) %>% 
        arrange(desc(pEC50))
    } else {
      data_list()$data_peptides %>% 
        mutate(log2.fc = signif(log2.fc, 3), 
               p.value = signif(p.value, 3)) %>% 
        arrange(p.value)
    }
  }, 
  rownames = FALSE)
  
  
  
  # ---- Plot protein sequence ---- 
  
  # Manage plot range 
  protein_sequence_plot_range <- reactiveValues()
  
  observeEvent(input$protein_sequence_plot_brush, {
    protein_sequence_plot_range$xmin <- floor(input$protein_sequence_plot_brush$xmin)
    protein_sequence_plot_range$xmax <- ceiling(input$protein_sequence_plot_brush$xmax)
  })
  
  observeEvent(c(input$protein_sequence_plot_dbclick, 
                 protein_group_selected()), {
                   protein_sequence_plot_range$xmin <- NULL
                   protein_sequence_plot_range$xmax <- NULL
                 })
  
  
  # Plot 
  output$protein_sequence <- renderPlot({
    
    # Define plot x-range
    if (is.null(protein_sequence_plot_range$xmin)) 
      protein_range <- c(1, 
                         filter(data_list()$data_protein_annotation, 
                                Protein.Group == 
                                  protein_group_selected())$Length[[1]])
    else 
      protein_range <- c(protein_sequence_plot_range$xmin, 
                         protein_sequence_plot_range$xmax)
    
    
    # Decide EC50 or log2 FC
    data_plot <- data_list()$data_peptides %>% 
      filter(Protein.Group == protein_group_selected()) %>% 
      mutate(label = paste0(Peptide.Id, " [", 
                            from, "-", to, "]"), 
             .by = "Peptide.Id")
    
    
    # Initialize protein sequence plot 
    p_seq <- plot_protein_seq(protein_range = protein_range, 
                              protein_width = 0.6, 
                              x_lab = paste0(
                                data_list()$data_peptides %>% 
                                  filter(Protein.Group == protein_group_selected()) %>% 
                                  pull(Genes) %>% 
                                  first(), 
                                " (", protein_range[1], "-", protein_range[2], ")"), 
                              y_lab = if("pEC50_signed" %in% names(data_plot)) 
                                expression(paste("signed local ", pEC[50])) 
                              else expression(paste(log[2], " fold change")), 
                              color = "#b0b0b0", 
                              alpha = 0.5, 
                              base_size = 16)
    
    
    # Add protein sites 
    if (protein_group_selected() %in% 
        data_list()$data_protein_annotation$Protein.Group &&
        !is.na(filter(data_list()$data_protein_annotation, 
                      Protein.Group == protein_group_selected())$`Binding site`)) {
      p_seq <- add_plot_protein_sites(
        p = p_seq, 
        data_sites = data_list()$data_protein_annotation %>% 
          filter(Protein.Group == protein_group_selected()) %>% 
          extract_UniProt_seqinfo("Binding site"), 
        start = "from", 
        end = "to", 
        name = "ligand", 
        add.length = 0.5, 
        merge_name_position = T, 
        legend.position = "right")
    }
    
    
    # Add peptides 
    if ("pEC50_signed" %in% names(data_plot)) {
      
      p_seq <- data_plot %>% 
        mutate(
          # Assign peptide groups to highlight
          peptide_group = case_when(
            Peptide.Id == peptide_selected() ~ "highlight", 
            .default = "regular"), 
          pEC50 = replace_na(pEC50_signed, 0)) %>% 
        plot_peptides_seq2(
          p = p_seq, 
          yvalue = "pEC50_signed", 
          Start = "from", 
          End = "to", 
          color = "regulation", 
          color.scale = peptide_colors, 
          protein_range = protein_range, 
          add.labels = protein_range[1] != 1, 
          label = "label", 
          nudge_x = 0.1, 
          nudge_y = 0, 
          hjust = 0,
          vjust = 0, 
          direction = "both", 
          min.segment.length = 2)
      
    } else {
      
      p_seq <- data_plot %>% 
        mutate(
          # Assign peptide groups to highlight
          peptide_group = case_when(
            Peptide.Id == peptide_selected() ~ "highlight", 
            .default = "regular"), 
          log2.fc = replace_na(log2.fc, 0)) %>% 
        plot_peptides_seq2(
          p = p_seq, 
          yvalue = "log2.fc", 
          Start = "from", 
          End = "to", 
          color = "regulation_peptide", 
          color.scale = peptide_colors, 
          protein_range = protein_range, 
          add.labels = protein_range[1] != 1, 
          label = "label", 
          nudge_x = 0.1, 
          nudge_y = 0, 
          hjust = 0,
          vjust = 0, 
          direction = "both", 
          min.segment.length = 2)
      
    }
    
    p_seq
    
  })
  
  
  # ---- Visualize protein structure ---- 
  protein_group_structure <- reactiveVal(protein_default)
  
  # observeEvent(protein_group_selected)
  
  
  output$protein_structure <- renderNGLVieweR({
    pdb <- NGLVieweR_AFdb(protein_group_selected()) 
    
    if (is.list(pdb)) {
      
      pdb <- pdb %>% 
        stageParameters(backgroundColor = "white", 
                        zoomSpeed = 2) %>%
        addRepresentation("cartoon", param = list(color = "#b0b0b0", 
                                                  opacity = 1))
      
      peptides_data <- data_list()$data_peptides %>% 
        filter(Protein.Group == protein_group_selected()) %>% 
        arrange(rank(match(regulation, names(peptide_colors)), 
                     na.last = F, 
                     ties.method = "first"))
      
      for (i in seq_len(nrow(peptides_data))) {
        
        selection <- paste0(peptides_data$from[i], "-", peptides_data$to[i])
        
        pdb <- pdb %>% 
          addRepresentation("cartoon", 
                            param = list(name = peptides_data$Peptide.Id[i], 
                                         sele = selection, 
                                         colorValue = peptide_colors[[peptides_data$regulation[i]]], 
                                         colorScheme = "element"))
      }
      
      # Add ATP binding site 
      if (protein_group_selected() %in% data_list()$data_protein_annotation$Protein.Group &&
          !is.na(filter(data_list()$data_protein_annotation, Protein.Group == protein_group_selected())$`Binding site`) && 
          str_detect(filter(data_list()$data_protein_annotation, Protein.Group == 
                            protein_group_selected())$`Binding site`, "ATP")) {
        
        atp_sites <- data_list()$data_protein_annotation %>% 
          filter(Protein.Group == protein_group_selected()) %>% 
          extract_UniProt_seqinfo("Binding site") %>% 
          filter(ligand == "ATP")
        
        for (i in seq_len(nrow(atp_sites))) {
          
          selection <- paste0(atp_sites$from[i], "-", atp_sites$to[i])
          
          pdb <- pdb %>% 
            addRepresentation("cartoon", 
                              param = list(name = atp_sites$`Binding site`[i], 
                                           sele = selection, 
                                           colorValue = "#18974C", 
                                           colorScheme = "element"))
        }
      }
    } else {
      pdb <- NGLVieweR_pdb("1XYZ") %>% 
        stageParameters(backgroundColor = "grey", 
                        zoomSpeed = 2)
    }
    
    pdb %>% 
      setQuality(quality = "high")
  }
  )
  
  
  # ---- Snapshot of structure ---- 
  # observeEvent(input$snapshot, {
  #   NGLVieweR_proxy("protein_structure") %>%
  #     snapShot("Snapshot",
  #              param = list(
  #                antialias = TRUE,
  #                trim = TRUE,
  #                transparent = TRUE,
  #                scale = 5))})
  # 
  
  
  # ---- Explorer mode ---- 
  
  ## Peptides ---- 
  
  ## Plot peptide curves 
  output$peptide_ds_curves <- renderPlot({
    if (peptide_selected_debounced() %in% data_list()$data_raw$Peptide.Id) {
      data_list()$data_raw %>% 
        filter(Peptide.Id == peptide_selected_debounced()) %>% 
        ggplot(aes(log10(concentration), Ratio, color = rep)) + 
        geom_smooth(method = "drm", 
                    method.args = list(fct = LL.4(), logDose = 10), 
                    se = FALSE, 
                    linewidth=0.5) + 
        geom_point(shape = 16, alpha = 0.8, size = 2) + 
        theme_classic(15) + 
        theme(legend.position = "bottom") + 
        labs(title = paste0("Peptide ", 
                            peptide_selected_debounced(), 
                            " (", 
                            filter(data_list()$data_peptides, Peptide.Id == 
                                     peptide_selected_debounced())$Genes[[1]], 
                            " ", 
                            filter(data_list()$data_peptides, Peptide.Id == 
                                     peptide_selected_debounced())$from[[1]], 
                            "-", 
                            filter(data_list()$data_peptides, Peptide.Id == 
                                     peptide_selected_debounced())$to[[1]], 
                            ")"), 
             color = "Replicate")
    } else {
      ggplot() + 
        theme_classic(15) + 
        labs(title = "No data available.")
    }
  })
  
  ## Generate peptide table 
  output$peptide_EC50s <- 
    renderTable(data_list()$data_peptides %>% 
                  filter(Peptide.Id == peptide_selected_debounced()) %>% 
                  select(Peptide.Id, 
                         EC50_rep1, EC50_rep2, EC50_rep3, EC50_rep4, EC50, 
                         R2_rep1, R2_rep2, R2_rep3, R2_rep4, CV_pEC50) %>% 
                  rename(pEC50_mean = EC50, 
                         R2_mean = CV_pEC50) %>% 
                  pivot_longer(
                    -Peptide.Id, 
                    names_to = c("set", ".value"),
                    names_pattern = "(EC50|R2)_(rep1|rep2|rep3|rep4|mean)") %>% 
                  select(-Peptide.Id) %>% 
                  mutate(across(where(is.numeric), \(x) round(x, 3))) %>% 
                  mutate(mean = paste(c("mean =", "CV ="), mean)) %>% 
                  rename(pEC50 = mean))
  
  
  
  
  
  ## Structures ----
  
  # Select PDB structure 
  observeEvent(c(protein_group_selected(), input$expmode_navset), {
    if (input$explorer_mode && input$expmode_navset == "pdb") {
      if (is.na(filter(data_list()$data_protein_annotation, 
                       Protein.Group == protein_group_selected())$PDB)) 
        pdb_ids <- c("AFdb" = protein_group_selected())
      else
        pdb_ids <- data_list()$data_protein_annotation %>% 
          filter(Protein.Group == protein_group_selected()) %>% 
          #mutate(PDB = str_replace_all(PDB, "\\.;", ".here;")) %>% 
          separate_longer_delim(cols = "PDB", delim = ".;") %>% 
          filter(PDB != "") %>% 
          separate_wider_delim(cols = "PDB", delim = "; ", 
                               names = c("pdb_id", 
                                         "method", 
                                         "resolution", 
                                         "subunits")) %>% 
          # mutate(title = get_pdb_title(pdb_id), .after = "pdb_id") %>% 
          # with(setNames(pdb_id, 
          #               paste0(pdb_id, " - ", 
          #                      title)))
          pull(pdb_id, pdb_id) %>% 
          c("AFdb" = protein_group_selected(), .)
      
      
      updateSelectizeInput(
        session,
        "expmode_pdb_id",
        choices = NULL,
        selected = NULL,
        server = TRUE
      )
      
      updateSelectizeInput(
        session, 
        "expmode_pdb_id", 
        paste0("Select PDB entry (", length(pdb_ids) - 1, " found):"), 
        choices = pdb_ids, 
        selected = pdb_ids[1], 
        server = T)
    }
  })
  
  output$expmode_pdb_title <- renderText({
    if (input$expmode_pdb_id != isolate(protein_group_selected()))
      get_pdb_title(input$expmode_pdb_id)
    else
      "AlphaFold pLDDT"
  })
  
  output$expmode_pdb_structure <- renderNGLVieweR({
    if (input$expmode_pdb_id != isolate(protein_group_selected()))
      pdb <- NGLVieweR_pdb(input$expmode_pdb_id) %>% 
        stageParameters(backgroundColor = "white", 
                        zoomSpeed = 2) %>%
        addRepresentation("cartoon")
    else {
      pdb <- NGLVieweR_AFdb(input$expmode_pdb_id) 
      
      if (is.list(pdb)) 
        pdb <- pdb %>% 
          stageParameters(backgroundColor = "white", 
                          zoomSpeed = 2) %>% 
          addRepresentation("cartoon", 
                            param = list(colorScheme = "bfactor"))
      else 
        pdb <- NGLVieweR_pdb("1XYZ") %>% 
          stageParameters(backgroundColor = "grey", 
                          zoomSpeed = 2)
    }
    
    pdb %>% 
      setQuality(quality = "high")
  })
  
  
}

# Start app
shinyApp(ui, server)
