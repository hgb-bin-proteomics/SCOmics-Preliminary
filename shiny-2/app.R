# scomics - a proof of concept implementation of automated visualization of multiomics results
# Copyright (C) 2025  Micha J. Birklbauer <micha.birklbauer@fh-hagenberg.at>

library(shiny)
library(bslib)

library(tidyverse)
library(proDA)
library(corrplot)
library(httr)
library(jsonlite)
library(ggsci)

get_location <- function(data, colname){
  
  gene <- data[colname]
  url1 <- paste("https://rest.uniprot.org/uniprotkb/search?query=gene:", gene, "+AND+organism_id:9606&format=tsv", sep = "")
  response <- GET(url1)
  if (status_code(response) != 200) 
    return(NA)
  url1_content <- content(response, "text", encoding = "utf-8")
  tryCatch(
    {
      url1_table = read.csv(text=url1_content, header = T, sep = "\t")
    },
    error = function(cond) {
      message("Could not get reponse for gene: ", gene)
      return(NA)
    }
  )
  if (nrow(url1_table) == 0)
    return(NA)
  accession = url1_table[1,"Entry"]
  
  url2 <- paste("https://rest.uniprot.org/uniprotkb/", accession, ".json", sep = "")
  response <- GET(url2)
  if (status_code(response) != 200)
    return(NA)
  url2_content <- content(response, "text", encoding = "utf-8")
  url2_json <- fromJSON(url2_content)
  if (is.null(url2_json$comments))
    return(NA)
  if(is.null(url2_json$comments$subcellularLocations))
    return(NA)
  idx <- match("SUBCELLULAR LOCATION", url2_json$comments$commentType)
  if (is.na(idx))
    return(NA)
  loc <- url2_json$comments$subcellularLocations[[idx]]$location$value[[1]][1]
  return(loc)
}

get_cv <- function(x){
  return(sd(x, na.rm = T) / mean(x, na.rm = T))
}

pp_str <- function(data, colname) {
  val <- data[colname]
  return(str_trim(str_split(val, ",")[[1]][1]))
}

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
          "Coefficients of Variation"
        ),
        card_body(
          plotOutput("plot_cv", height = 500, fill = F),
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
    layout_columns(
      min_height = 600,
      card(
        min_height = 600,
        card_header(
          "CCO Distribution Proteomics"
        ),
        card_body(
          plotOutput("plot_cco_prot", height = 500, fill = F),
        )
      ),
      card(
        min_height = 600,
        card_header(
          "CCO Distribution Transcriptomics"
        ),
        card_body(
          plotOutput("plot_cco_trans", height = 500, fill = F),
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
  
  #### GLOBAL ####
  top_cco <- 9
  
  data <- reactive({
    
    if(is.null(input$data_prot)){
      return("")
    }
    
    if(is.null(input$data_trans)){
      return("")
    }
    
    return(0)
  })
  
  output$plot_cv <- renderPlot({
    req(data())
    
    na_drop_threshold = 2
    normalize_proteomics = F
    normalize_transcriptomics = F
    use = if (input$rri) SC_RRI else SC_noRRI
    colnames_rri = c("cell e", "cell d", "cell c", "cell b", "cell a")
    colnames_norri = c("cell m", "cell l", "cell k", "cell j", "cell i")
    
    colnames_use <- if (input$rri) colnames_rri else colnames_norri
    
    #### PROTEOMICS ####
    
    data <- read.table(input$data_prot$datapath,
                       header = T,
                       sep = "\t")
    keep_cols <- c(grep(pattern = "PG.ProteinGroups", colnames(data)),
                   grep(pattern = "PG.Genes", colnames(data)),
                   grep(pattern = "*MS1Quantity", colnames(data)))
    data <- data[,keep_cols]
    
    # remove missing values
    #data <- data[complete.cases(data),]
    data <- data[!rowSums(is.na(data)) > na_drop_threshold,]
    
    # remove contaminants
    #data <- data |> filter(!grepl("cont_*", PG.ProteinGroups))
    
    # median normalization
    data[,3:7] <- if (normalize_proteomics) proDA::median_normalization(data.matrix(data[,3:7])) else data[,3:7]
    data <- data[,3:7]
    colnames(data) <- colnames_use
    
    data$cv <- apply(data, 1, get_cv)
    data$mean_cv <- mean(data$cv)
    data$sd_cv <- sd(data$cv)
    sample_size_prot = nrow(data)
    data$group = rep("Proteomics", sample_size_prot)
    
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
    #data_ngs <- data_ngs[complete.cases(data_ngs),]
    data_ngs <- data_ngs[!rowSums(is.na(data_ngs)) > na_drop_threshold,]
    
    # median normalization
    data_ngs[,2:6] <- if (normalize_transcriptomics) proDA::median_normalization(data.matrix(data_ngs[,2:6])) else data_ngs[,2:6]
    data_ngs <- data_ngs[,2:6]
    colnames(data_ngs) <- colnames_use
    
    data_ngs$cv <- apply(data_ngs, 1, get_cv)
    data_ngs$mean_cv <- mean(data_ngs$cv)
    data_ngs$sd_cv <- sd(data_ngs$cv)
    sample_size_ngs = nrow(data_ngs)
    data_ngs$group = rep("Transcriptomics", sample_size_ngs)
    
    #### CV PLOT ####
    
    caption = paste0("Coefficients of Variation for Proteomics (n = ", sample_size_prot, ") and Transcriptomics (n = ", sample_size_ngs, ").")
    
    df <- rbind(data, data_ngs)
    ggplot(df, aes(x=group, y=cv, fill=group)) + 
      geom_violin(trim = F, colour = NA) + 
      geom_text(aes(x=group, y=mean_cv, label=round(mean_cv, 2)), hjust = -1) +
      geom_point(aes(x=group, y=mean_cv), shape = 18, size = 5, fill = "black") +
      geom_errorbar(aes(x=group, ymin=mean_cv-sd_cv, ymax=mean_cv+sd_cv), width = 0.1) +
      scale_fill_manual(values = alpha(c("#6d4aff", "#c999ff"), 0.8),
                        labels = c("Proteomics", "Transcriptomics")) +
      xlab("Analysis Type") +
      ylab("Coefficient of Variation (CV)") +
      labs(fill="Analysis Type", caption = caption) +
      theme_minimal(base_size = 18) + 
      theme(plot.caption.position = "plot",
            plot.caption = element_text(color = "grey", face = "italic", hjust = 0)) +
      guides(fill = "none") #guide_legend(ncol = 1, label.position = "right"))
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
  
  output$plot_cco_prot <- renderPlot({
    req(data())
    
    #### PROTEOMICS ####
    
    data <- read.table(input$data_prot$datapath,
                       header = T,
                       sep = "\t")
    keep_cols <- c(grep(pattern = "PG.ProteinGroups", colnames(data)),
                   grep(pattern = "PG.Genes", colnames(data)),
                   grep(pattern = "*MS1Quantity", colnames(data)))
    data <- data[,keep_cols]
    
    # remove data with no gene annotation
    data[data==""] = NA
    data <- data[complete.cases(data$PG.Genes),]
    
    # remove ambigous PGs
    data <- data |> filter(!grepl(";", PG.Genes))
    
    data$CCO <- apply(data, 1, get_location, colname = "PG.Genes")
    data <- data[complete.cases(data$CCO),]

    data$CCO <- apply(data, 1, pp_str, colname = "CCO")
    data$CCO <- as.factor(data$CCO)
    
    cco <- as.data.frame(table(data$CCO))
    cco <- cco[order(-cco$Freq),]
    cco <- cco[1:top_cco,]
    cco <- rbind(cco, data.frame(Var1="Other",Freq=nrow(data)-sum(cco$Freq)))
    
    ggplot(cco, aes(x=reorder(Var1, -Freq), y=Freq, fill=Var1)) +
      geom_bar(stat="identity", color="black", width = 0.7) +
      scale_fill_cosmic() +
      xlab("Cellular Component Ontology (CCO)") +
      ylab("Number of Proteins") +
      theme_minimal(base_size = 18) +
      theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
      guides(fill = "none") #guide_legend(ncol = 1, label.position = "right"))
  }) |>
    bindEvent(input$do)
  
  output$plot_cco_trans <- renderPlot({
    req(data())
    
    use = if (input$rri) SC_RRI else SC_noRRI
    #### RNA-Seq ####
    
    data_ngs <- read.table(input$data_trans$datapath,
                           header = T,
                           sep = "\t")
    
    keep_cols <- c("gene_name",
                   paste("X", use, sep = ""))
    data_ngs <- data_ngs[,keep_cols]
    
    # remove data with no gene annotation and too many missing values
    data_ngs[data_ngs==0] = NA
    data_ngs[data_ngs==""] = NA
    data_ngs <- data_ngs[complete.cases(data_ngs$gene_name),]
    data_ngs <- data_ngs[!rowSums(is.na(data_ngs)) > 1,]
    
    data_ngs$CCO <- apply(data_ngs, 1, get_location, colname = "gene_name")
    data <- data_ngs[complete.cases(data_ngs$CCO),]
    
    data$CCO <- apply(data, 1, pp_str, colname = "CCO")
    data$CCO <- as.factor(data$CCO)
    
    cco <- as.data.frame(table(data$CCO))
    cco <- cco[order(-cco$Freq),]
    cco <- cco[1:top_cco,]
    cco <- rbind(cco, data.frame(Var1="Other",Freq=nrow(data)-sum(cco$Freq)))
    
    ggplot(cco, aes(x=reorder(Var1, -Freq), y=Freq, fill=Var1)) +
      geom_bar(stat="identity", color="black", width = 0.7) +
      scale_fill_cosmic() +
      xlab("Cellular Component Ontology (CCO)") +
      ylab("Number of Genes") +
      theme_minimal(base_size = 18) +
      theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
      guides(fill = "none") #guide_legend(ncol = 1, label.position = "right"))
  }) |>
    bindEvent(input$do)
  
}

shinyApp(ui = ui, server = server)