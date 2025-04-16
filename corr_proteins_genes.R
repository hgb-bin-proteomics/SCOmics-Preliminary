library(tidyverse)
library(proDA)
library(corrplot)

#### RNAseq Sample Annotation ####
# A1-E1
SC_RRI <- c(345275:345279)
# A2-E2
SC_noRRI <- c(345283:345287)

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

#pg_cor <- cor(cor_data, method = "pearson", use = "pairwise.complete.obs")
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
#corrplot(pg_cor, method = "color", type = "lower", tl.pos = "n")
#corrplot(pg_cor, method = "color", order = "hclust", hclust.method = "centroid", tl.pos = "n")
#heatmap(pg_cor, Rowv = NA, Colv = NA)

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

#g_cor <- cor(cor_data_ngs, method = "pearson", use = "pairwise.complete.obs")
#corrplot(g_cor, method = "color", order = "hclust", hclust.method = "centroid", tl.pos = "n")

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
#corrplot(cor_multi, method = "color", type = "lower", tl.pos = "n")
#corrplot(cor_multi, method = "color", order = "hclust", hclust.method = "centroid", tl.pos = "n")

cor_multi <- cor_multi[1:length(genes),(length(genes)+1):(length(genes)*2)]
#corrplot(cor_multi, method = "color", type = "lower", tl.pos = "n")
corrplot(cor_multi, method = "color", order = "hclust", hclust.method = "centroid", tl.pos = "n")
