library(tidyverse)
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

pp_str <- function(data, colname) {
  val <- data[colname]
  return(str_trim(str_split(val, ",")[[1]][1]))
}

#### RNAseq Sample Annotation ####
# Prot: B20-B24
SC_RRI <- c(345279:345275)
# Prot: A20-A24
SC_noRRI <- c(345287:345283)

#### PARAMETERS ####

top_cco = 9
na_drop_threshold = 2

RRI = "data/spectronaut/20250325_111740_MM_A1_SCMO001_1cell_split_+RRI_directDIA_default/MM_A1_SCMO001_1cell_split_+RRI_directDIA_default_Report_Protein_JB (Pivot).tsv"
noRRI = "data/spectronaut/20250325_111639_MM_A1_SCMO001_1cell_split_directDIA_default/MM_A1_SCMO001_1cell_split_directDIA_default_Report_Protein_JB (Pivot).tsv"
proteomics_data = RRI

counts <- "data/ngs/salmon.merged.gene_counts.tsv"
counts_scaled <- "data/ngs/salmon.merged.gene_counts_scaled.tsv"
counts_scaled_len <- "data/ngs/salmon.merged.gene_counts_length_scaled.tsv"
transcriptomics_data <- counts
use <- SC_RRI

#### PROTEOMICS ####

data <- read.table(proteomics_data,
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

#data$CCO <- apply(data, 1, get_location, colname = "PG.Genes")
#data <- data[complete.cases(data$CCO),]
#saveRDS(data, "data/cco_prot.rds")

data <- readRDS("data/cco_prot.rds")
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
  ylab("CCO Distribution in Number of Proteins") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
  guides(fill = "none") #guide_legend(ncol = 1, label.position = "right"))

#### RNA-Seq ####

data_ngs <- read.table(transcriptomics_data,
                       header = T,
                       sep = "\t")

keep_cols <- c("gene_name",
               paste("X", use, sep = ""))
data_ngs <- data_ngs[,keep_cols]

# remove data with no gene annotation and too many missing values
data_ngs[data_ngs==0] = NA
data_ngs[data_ngs==""] = NA
data_ngs <- data_ngs[complete.cases(data_ngs$gene_name),]
data_ngs <- data_ngs[!rowSums(is.na(data_ngs)) > na_drop_threshold,]

#data_ngs$CCO <- apply(data_ngs, 1, get_location, colname = "gene_name")
#data_ngs <- data_ngs[complete.cases(data_ngs$CCO),]
#saveRDS(data_ngs, "data/cco_trans.rds")

data <- readRDS("data/cco_trans.rds")
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
  ylab("CCO Distribution in Number of Genes") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
  guides(fill = "none") #guide_legend(ncol = 1, label.position = "right"))