# R Shiny script to plot fungal ITS2 metabarcoding data
# Craig F. Barrett, Department of Biology, West Virginia University

# Install packages

# Base CRAN packages
# install.packages(c(
#   "shiny", "phyloseq", "ggplot2", "plotly", 
#   "dplyr", "tidyr", "ggh4x"
# ))

# You may need BiocManager to install phyloseq
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("phyloseq")

# Citations

# Core packages
# citation("shiny")     # RStudio (2021)
# citation("phyloseq")  # McMurdie & Holmes (2013)
# citation("ggplot2")   # Wickham (2016)
# citation("plotly")    # Sievert (2020)
# citation("dplyr")     # Wickham et al. (2023)
# citation("tidyr")     # Wickham & Girlich (2022)

# Script development and debugging was assisted by R Wizard, an R-specialized GPT-4 assistant developed using OpenAI’s ChatGPT platform.
# OpenAI (2025). ChatGPT (GPT-4o). Retrieved from https://chat.openai.com/


library(shiny)
library(phyloseq)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)

ui <- fluidPage(
  titlePanel("Interactive Phyloseq Barplot Viewer"),
  sidebarLayout(
    sidebarPanel(
      fileInput("phyloseq_file", "Upload phyloseq .rds file", accept = ".rds"),
      selectInput("tax_rank", "Select taxonomic rank", choices = NULL),
      actionButton("uncheck_all", "Uncheck all taxa"),
      uiOutput("taxa_ui"),
      uiOutput("subtaxa_ui"),
      selectInput("facet_var", "Facet by (optional)", choices = NULL, selected = NULL),
      numericInput("facet_ncol", "Facets per row (ncol)", value = NA, min = 1),
      hr(),
      numericInput("pdf_width", "PDF width (in)", value = 11, min = 4),
      numericInput("pdf_height", "PDF height (in)", value = 8.5, min = 4),
      downloadButton("download_plot", "Download PDF")
    ),
    mainPanel(
      plotlyOutput("barplot", height = "1000px")
    )
  )
)

server <- function(input, output, session) {
  ps_data <- reactive({
    req(input$phyloseq_file)
    readRDS(input$phyloseq_file$datapath)
  })

  observeEvent(ps_data(), {
    ranks <- rank_names(ps_data())
    updateSelectInput(session, "tax_rank", choices = ranks, selected = tail(ranks, 1))

    sample_vars <- sample_variables(ps_data())
    updateSelectInput(session, "facet_var", choices = c("None" = "", sample_vars))
  })

  all_taxa <- reactiveVal()

  output$taxa_ui <- renderUI({
    req(input$tax_rank)
    taxa <- unique(tax_table(ps_data())[, input$tax_rank])
    taxa <- as.character(taxa[!is.na(taxa)])
    all_taxa(taxa)
    checkboxGroupInput("taxa", "Select taxa to display", choices = taxa, selected = taxa)
  })

  output$subtaxa_ui <- renderUI({
    req(input$tax_rank, input$taxa)

    ranks <- rank_names(ps_data())
    current_rank_index <- which(ranks == input$tax_rank)
    if (current_rank_index >= length(ranks)) return(NULL)

    child_rank <- ranks[current_rank_index + 1]
    tax_table_df <- as.data.frame(tax_table(ps_data()))
    out <- tagList()

    for (parent in input$taxa) {
      child_taxa <- tax_table_df %>%
        filter(.data[[input$tax_rank]] == parent) %>%
        pull(.data[[child_rank]]) %>%
        unique() %>%
        na.omit() %>%
        as.character()

      if (length(child_taxa) > 0) {
        out <- tagList(
          out,
          checkboxGroupInput(
            inputId = paste0("subtaxa_", parent),
            label = paste(child_rank, "under", parent),
            choices = child_taxa,
            selected = child_taxa
          )
        )
      }
    }

    out
  })

  observeEvent(input$uncheck_all, {
    updateCheckboxGroupInput(session, "taxa", selected = character(0))
  })

  plot_data <- reactive({
    req(input$tax_rank, input$taxa)

    ps <- ps_data()
    ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
    df <- psmelt(ps_rel)

    df <- df[df[[input$tax_rank]] %in% input$taxa, ]

    ranks <- rank_names(ps_data())
    current_rank_index <- which(ranks == input$tax_rank)
    if (current_rank_index < length(ranks)) {
      child_rank <- ranks[current_rank_index + 1]
      selected_subtaxa <- unlist(
        lapply(input$taxa, function(taxon) input[[paste0("subtaxa_", taxon)]])
      )
      if (!is.null(selected_subtaxa) && length(selected_subtaxa) > 0) {
        df <- df[df[[child_rank]] %in% selected_subtaxa, ]
      }
    }

    df <- df %>%
      group_by(Sample) %>%
      filter(sum(Abundance) > 0) %>%
      ungroup()

    tax <- as(tax_table(ps), "matrix")
    tax_strings <- apply(tax, 1, function(x) paste(names(x), na.omit(x), sep = ": ", collapse = "; "))
    df$tooltip <- paste("Sample:", df$Sample, "<br>", tax_strings[match(df$OTU, rownames(tax))])

    df$rank_color <- df[[input$tax_rank]]

    df
  })

  static_plot <- reactive({
    df <- plot_data()
    req(nrow(df) > 0)

    if (input$facet_var != "") {
      df <- df[!is.na(df[[input$facet_var]]), ]
      df[[input$facet_var]] <- factor(df[[input$facet_var]])
      df[[input$facet_var]] <- droplevels(df[[input$facet_var]])
    }

    p <- ggplot(
      df,
      aes(
        x = Sample,
        y = Abundance,
        fill = rank_color,
        subgroup = OTU,
        text = tooltip
      )
    ) +
      geom_bar(
        stat = "identity",
        position = "stack",
        colour = "black",
        linewidth = 0.1
      ) +
      labs(
        x = "Sample",
        y = "Relative Abundance",
        fill = input$tax_rank
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10),
        strip.background = element_rect(fill = "grey90")
      )

    if (input$facet_var != "") {
      if (!is.na(input$facet_ncol)) {
        p <- p + facet_wrap(as.formula(paste("~", input$facet_var)), scales = "free_x", ncol = input$facet_ncol)
      } else {
        p <- p + facet_wrap(as.formula(paste("~", input$facet_var)), scales = "free_x")
      }
    }

    p
  })

  output$barplot <- renderPlotly({
    ggplotly(static_plot(), tooltip = "text") %>% layout(barmode = "stack")
  })

  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("phyloseq_barplot_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      ggsave(file, plot = static_plot(), width = input$pdf_width, height = input$pdf_height)
    }
  )
}

shinyApp(ui = ui, server = server)

