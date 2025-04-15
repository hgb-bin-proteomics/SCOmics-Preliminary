library(tidyverse)

counts <- read.table("data/ngs/salmon.merged.gene_counts.tsv",
                     header = T,
                     sep = "\t")
counts_cols <- colnames(counts)
summary(counts)

counts_norm <- read.table("data/ngs/salmon.merged.gene_counts_scaled.tsv",
                          header = T,
                          sep = "\t")
counts_norm_cols <- colnames(counts_norm)
summary(counts_norm)

counts_len_norm <- read.table("data/ngs/salmon.merged.gene_counts_length_scaled.tsv",
                              header = T,
                              sep = "\t")
counts_len_norm_cols <- colnames(counts_len_norm)
summary(counts_len_norm)

# A1-E1
SC_RRI <- c(345275:345279)
# A2-E2
SC_noRRI <- c(345283:345287)
