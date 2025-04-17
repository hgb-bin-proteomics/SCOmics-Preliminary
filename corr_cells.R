library(tidyverse)
library(proDA)
library(corrplot)

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

normalize_proteomics = T
normalize_transcriptomics = T

# set row names for plot
rownames_prot <- c("cell e PROT", "cell d PROT", "cell c PROT", "cell b PROT", "cell a PROT")
rownames_rna <- c("cell e RNA", "cell d RNA", "cell c RNA", "cell b RNA", "cell a RNA")

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

# remove ambigous PGs
data <- data |> filter(!grepl(";", PG.Genes))

# get genes
protein_genes <- data$PG.Genes

# median normalization
# described in proDA package
# https://bioconductor.org/packages/release/bioc/vignettes/proDA/inst/doc/Introduction.html
data_m <- data.matrix(data[,3:7])
normalized <- if (normalize_proteomics) proDA::median_normalization(data_m) else data_m

cor_data <- data.frame(t(normalized))
colnames(cor_data) <- data$PG.ProteinGroups

#### RNA-Seq ####

data_ngs <- read.table(transcriptomics_data,
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
normalized_ngs <- if (normalize_transcriptomics) proDA::median_normalization(data_m_ngs) else data_m_ngs

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
#corrplot(cor_multi, method = "color", tl.col = "black")
corrplot(cor_multi, method = "color", tl.col = "black", addCoef.col = "white")
#corrplot(cor_multi, method = "color", order = "hclust", hclust.method = "centroid", tl.col = "black")
#corrplot(cor_multi, method = "color", order = "hclust", hclust.method = "centroid", tl.col = "black", addCoef.col = "white")
