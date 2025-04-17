library(tidyverse)
library(proDA)

get_cv <- function(x){
  return(sd(x)/mean(x))
}

#### RNAseq Sample Annotation ####
# Prot: B20-B24
SC_RRI <- c(345279:345275)
# Prot: A20-A24
SC_noRRI <- c(345287:345283)

#### PARAMETERS ####

RRI = "data/spectronaut/20250325_111740_MM_A1_SCMO001_1cell_split_+RRI_directDIA_default/MM_A1_SCMO001_1cell_split_+RRI_directDIA_default_Report_Protein_JB (Pivot).tsv"
noRRI = "data/spectronaut/20250325_111639_MM_A1_SCMO001_1cell_split_directDIA_default/MM_A1_SCMO001_1cell_split_directDIA_default_Report_Protein_JB (Pivot).tsv"
proteomics_data = RRI

counts <- "data/ngs/salmon.merged.gene_counts.tsv"
counts_scaled <- "data/ngs/salmon.merged.gene_counts_scaled.tsv"
counts_scaled_len <- "data/ngs/salmon.merged.gene_counts_length_scaled.tsv"
transcriptomics_data <- counts
use <- SC_RRI

normalize_proteomics = F
normalize_transcriptomics = F

#### PROTEOMICS ####

data <- read.table(proteomics_data,
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

# median normalization
data_m <- data.matrix(data[,3:7])
normalized <- if (normalize_proteomics) proDA::median_normalization(data_m) else data_m

cor_data <- data.frame(t(normalized))
colnames(cor_data) <- data$PG.ProteinGroups

cv <- apply(cor_data, 1, get_cv)

#### RNA-Seq ####

data_ngs <- read.table(transcriptomics_data,
                       header = T,
                       sep = "\t")

keep_cols <- c("gene_name",
               paste("X", use, sep = ""))
data_ngs <- data_ngs[,keep_cols]
summary(data_ngs)

# remove missing values
# real zero or technical zeros?
data_ngs[data_ngs==0] = NA
data_ngs <- data_ngs[complete.cases(data_ngs),]


# median normalization
data_m_ngs <- data.matrix(data_ngs[,2:6])
normalized_ngs <- if (normalize_transcriptomics) proDA::median_normalization(data_m_ngs) else data_m_ngs

cor_data_ngs <- data.frame(t(normalized_ngs))
colnames(cor_data_ngs) <- data_ngs$gene_name

cv_ngs <- apply(cor_data_ngs, 1, get_cv)

# df <- data.frame(cv = c(cv, cv_ngs), group = c(rep("Proteomics", 5), rep("Transcriptomics", 5)))
df <- data.frame(cv = cv, group = rep("Proteomics", 5))
ggplot(df, aes(x=group, y=cv, fill=group)) + geom_violin(trim = F) + geom_jitter() + scale_fill_manual(values = alpha(c("#6d4aff"), 0.5),labels = c("Proteomics"))
