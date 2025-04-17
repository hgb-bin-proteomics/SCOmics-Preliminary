# scomics - a proof of concept implementation of automated visualization of multiomics results
# Copyright (C) 2025  Micha J. Birklbauer <micha.birklbauer@fh-hagenberg.at>

library(shiny)
library(bslib)

library(tidyverse)
library(proDA)
library(corrplot)

ui <- page_sidebar(
  theme = bs_theme(bootswatch = "lux"),
  title = "Automated Single Cell Proteomics and Transcriptomics Analysis",
  
  sidebar = sidebar(
    title = "Data Input",
    width = 500,
    open = "open",
  
    tags$b("Proteomics"),
    
    fileInput(
      inputId = "data_prot",
      label = "MS1 Protein Quantification from Spectronaut"
    ),
    
    tags$b("Transcriptomics"),
    
    fileInput(
      inputId = "data_trans",
      label = "Gene Counts from Salmon"
    ),
    
    tags$b("Parameters"),
    
    input_switch(id = "rri",
                 label = "RRI Sample",
                 value = T
    ),
    
    input_switch(id = "pn",
                 label = "Normalize Proteomics Data",
                 value = T
    ),
    
    input_switch(id = "tn",
                 label = "Normalize Transcriptomics Data",
                 value = T
    ),
    
    actionButton(inputId = "do",
                 label = "Run Analysis!"
    ),
  ),
  
  tags$h3("Results"),
  
  tags$hr(),

  card(
    min_height = 1300,
    card_header("Visualizations"),
    layout_columns(
      min_height = 600,
      card(
        min_height = 600,
        card_header(
          "Protein Group Correlation"
        ),
        card_body(
          plotOutput("plot_prot", height = 500, fill = F),
        )
      ),
      card(
        min_height = 600,
        card_header(
          "Gene Correlation"
        ),
        card_body(
          plotOutput("plot_trans", height = 500, fill = F),
        )
      ),
    ),
    layout_columns(
      min_height = 600,
      card(
        min_height = 600,
        card_header(
          "Protein-Gene Correlation"
        ),
        card_body(
          plotOutput("plot_multi", height = 500, fill = F),
        )
      ),
      card(
        min_height = 600,
        card_header(
          "Cell Correlation"
        ),
        card_body(
          plotOutput("plot_cell", height = 500, fill = F),
        )
      ),
    ),
  ),
)

server <- function(input, output){
  
  #### RNAseq Sample Annotation ####
  # Prot: B20-B24
  SC_RRI <- c(345279:345275)
  # Prot: A20-A24
  SC_noRRI <- c(345287:345283)
  
  data <- reactive({
    
    if(is.null(input$data_prot)){
      return("")
    }
    
    if(is.null(input$data_trans)){
      return("")
    }
    
    return(0)
  })
  
  output$plot_prot <- renderPlot({
    req(data())
    
    # data input
    data <- read.table(input$data_prot$datapath,
                       header = T,
                       sep = "\t")
    keep_cols <- c(grep(pattern = "PG.ProteinGroups", colnames(data)),
                   grep(pattern = "PG.Genes", colnames(data)),
                   grep(pattern = "*MS1Quantity", colnames(data)))
    data <- data[,keep_cols]
    
    # remove missing values
    data <- data[complete.cases(data),]
    
    # remove contaminants
    data <- data |> filter(!grepl("cont_*", PG.ProteinGroups))
    
    # remove ambigous PGs
    data <- data |> filter(!grepl(";", PG.Genes))
    
    # get genes
    protein_genes <- data$PG.Genes
    
    # median normalization
    data_m <- data.matrix(data[,3:7])
    normalized <- if (input$pn) proDA::median_normalization(data_m) else data_m
    
    cor_data <- data.frame(t(normalized))
    colnames(cor_data) <- data$PG.ProteinGroups
    
    pg_cor <- cor(cor_data, method = "pearson", use = "pairwise.complete.obs")
    corrplot(pg_cor, method = "color", order = "hclust", hclust.method = "centroid", tl.pos = "n")
  }) |>
    bindEvent(input$do)
  
  output$plot_trans <- renderPlot({
    req(data())
    
    # data input
    data_ngs <- read.table(input$data_trans$datapath,
                           header = T,
                           sep = "\t")
    use = if (input$rri) SC_RRI else SC_noRRI
    
    keep_cols <- c("gene_name",
                   paste("X", use, sep = ""))
    data_ngs <- data_ngs[,keep_cols]
    
    # remove missing values
    # real zero or technical zeros?
    data_ngs[data_ngs==0] = NA
    data_ngs <- data_ngs[complete.cases(data_ngs),]
    
    # get genes
    rna_genes <- data_ngs$gene_name
    
    # median normalization
    data_m_ngs <- data.matrix(data_ngs[,2:6])
    normalized_ngs <- if (input$tn) proDA::median_normalization(data_m_ngs) else data_m_ngs
    
    cor_data_ngs <- data.frame(t(normalized_ngs))
    colnames(cor_data_ngs) <- data_ngs$gene_name
    
    g_cor <- cor(cor_data_ngs, method = "pearson", use = "pairwise.complete.obs")
    corrplot(g_cor, method = "color", order = "hclust", hclust.method = "centroid", tl.pos = "n")
  }) |>
    bindEvent(input$do)
  
  output$plot_multi <- renderPlot({
    req(data())
    
    #### PROTEOMICS ####
    
    data <- read.table(input$data_prot$datapath,
                       header = T,
                       sep = "\t")
    keep_cols <- c(grep(pattern = "PG.ProteinGroups", colnames(data)),
                   grep(pattern = "PG.Genes", colnames(data)),
                   grep(pattern = "*MS1Quantity", colnames(data)))
    data <- data[,keep_cols]
    summary(data)
    
    # remove missing values
    data <- data[complete.cases(data),]
    
    # remove contaminants
    data <- data |> filter(!grepl("cont_*", PG.ProteinGroups))
    
    # remove ambigous PGs
    data <- data |> filter(!grepl(";", PG.Genes))
    
    # get genes
    protein_genes <- data$PG.Genes
    
    # median normalization
    data_m <- data.matrix(data[,3:7])
    normalized <- if (input$pn) proDA::median_normalization(data_m) else data_m
    
    cor_data <- data.frame(t(normalized))
    colnames(cor_data) <- data$PG.ProteinGroups
    
    #### RNA-Seq ####
    
    data_ngs <- read.table(input$data_trans$datapath,
                           header = T,
                           sep = "\t")
    use = if (input$rri) SC_RRI else SC_noRRI
    
    keep_cols <- c("gene_name",
                   paste("X", use, sep = ""))
    data_ngs <- data_ngs[,keep_cols]
    
    # remove missing values
    # real zero or technical zeros?
    data_ngs[data_ngs==0] = NA
    data_ngs <- data_ngs[complete.cases(data_ngs),]
    
    # get genes
    rna_genes <- data_ngs$gene_name
    
    # median normalization
    data_m_ngs <- data.matrix(data_ngs[,2:6])
    normalized_ngs <- if (input$tn) proDA::median_normalization(data_m_ngs) else data_m_ngs
    
    cor_data_ngs <- data.frame(t(normalized_ngs))
    colnames(cor_data_ngs) <- data_ngs$gene_name
    
    #### INTEGRATION ####
    
    row.names(cor_data_ngs) <- row.names(cor_data)
    colnames(cor_data) <- data$PG.Genes
    
    genes <- base::intersect(protein_genes, rna_genes)
    
    cor_data <- cor_data[,genes]
    cor_data_ngs <- cor_data_ngs[,genes]
    
    colnames(cor_data) <- paste("Protein_", colnames(cor_data), sep = "")
    colnames(cor_data_ngs) <- paste("RNA_", colnames(cor_data_ngs), sep = "")
    
    cor_data_multi <- cbind(cor_data, cor_data_ngs)
    cor_multi <- cor(cor_data_multi, method = "pearson", use = "pairwise.complete.obs")
    cor_multi <- cor_multi[1:length(genes),(length(genes)+1):(length(genes)*2)]
    corrplot(cor_multi, method = "color", order = "hclust", hclust.method = "centroid", tl.pos = "n")
  }) |>
    bindEvent(input$do)
  
  output$plot_cell <- renderPlot({
    req(data())
    
    use = if (input$rri) SC_RRI else SC_noRRI
    rownames_rri = c("cell e", "cell d", "cell c", "cell b", "cell a")
    rownames_norri = c("cell m", "cell l", "cell k", "cell j", "cell i")
    
    # set row names for plot
    rownames_prot <- if (input$rri) paste(rownames_rri, "PROT", sep = " ") else paste(rownames_norri, "PROT", sep = " ")
    rownames_rna <- if (input$rri) paste(rownames_rri, "RNA", sep = " ") else paste(rownames_norri, "RNA", sep = " ")
    
    #### PROTEOMICS ####
    
    data <- read.table(input$data_prot$datapath,
                       header = T,
                       sep = "\t")
    keep_cols <- c(grep(pattern = "PG.ProteinGroups", colnames(data)),
                   grep(pattern = "PG.Genes", colnames(data)),
                   grep(pattern = "*MS1Quantity", colnames(data)))
    data <- data[,keep_cols]
    
    # remove missing values
    data <- data[complete.cases(data),]
    
    # remove contaminants
    data <- data |> filter(!grepl("cont_*", PG.ProteinGroups))
    
    # remove ambigous PGs
    data <- data |> filter(!grepl(";", PG.Genes))
    
    # get genes
    protein_genes <- data$PG.Genes
    
    # median normalization
    data_m <- data.matrix(data[,3:7])
    normalized <- if (input$pn) proDA::median_normalization(data_m) else data_m
    
    cor_data <- data.frame(t(normalized))
    colnames(cor_data) <- data$PG.ProteinGroups
    
    #### RNA-Seq ####
    
    data_ngs <- read.table(input$data_trans$datapath,
                           header = T,
                           sep = "\t")
    
    keep_cols <- c("gene_name",
                   paste("X", use, sep = ""))
    data_ngs <- data_ngs[,keep_cols]
    
    # remove missing values
    # real zero or technical zeros?
    data_ngs[data_ngs==0] = NA
    data_ngs <- data_ngs[complete.cases(data_ngs),]
    
    # get genes
    rna_genes <- data_ngs$gene_name
    
    # median normalization
    data_m_ngs <- data.matrix(data_ngs[,2:6])
    normalized_ngs <- if (input$tn) proDA::median_normalization(data_m_ngs) else data_m_ngs
    
    cor_data_ngs <- data.frame(t(normalized_ngs))
    colnames(cor_data_ngs) <- data_ngs$gene_name
    
    #### INTEGRATION ####
    
    colnames(cor_data) <- data$PG.Genes
    row.names(cor_data_ngs) <- rownames_rna
    row.names(cor_data) <- rownames_prot
    
    cor_data <- cor_data[sort(rownames_prot),]
    cor_data_ngs <- cor_data_ngs[sort(rownames_rna),]
    
    genes <- base::intersect(protein_genes, rna_genes)
    
    cor_data <- cor_data[,genes]
    cor_data_ngs <- cor_data_ngs[,genes]
    
    cor_data_multi <- t(rbind(cor_data, cor_data_ngs))
    cor_multi <- cor(cor_data_multi, method = "pearson", use = "pairwise.complete.obs")
    corrplot(cor_multi, method = "color", tl.col = "black", addCoef.col = "white")
  }) |>
    bindEvent(input$do)
  
}

shinyApp(ui = ui, server = server)