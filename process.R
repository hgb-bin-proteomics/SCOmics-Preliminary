library(tidyverse)
library(proDA)
library(corrplot)

# 20250325_111639_MM_A1_SCMO001_1cell_split_directDIA_default
data <- read.table("data/spectronaut/20250325_111639_MM_A1_SCMO001_1cell_split_directDIA_default/MM_A1_SCMO001_1cell_split_directDIA_default_Report_Protein_JB (Pivot).tsv",
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
# described in proDA package
# https://bioconductor.org/packages/release/bioc/vignettes/proDA/inst/doc/Introduction.html
data_m <- data.matrix(data[,3:7])
normalized <- proDA::median_normalization(data_m)

cor_data <- data.frame(t(normalized))
colnames(cor_data) <- data$PG.ProteinGroups

pg_cor <- cor(cor_data, method = "pearson", use = "pairwise.complete.obs")
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
corrplot(pg_cor, method = "color", type = "lower", tl.pos = "n")
#heatmap(pg_cor, Rowv = NA, Colv = NA)

