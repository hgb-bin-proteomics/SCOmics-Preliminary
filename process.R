library(tidyverse)

# 20250324_170201_MM_A1_SCMO001_directDIA_matching_default_set
data1 <- read.table("data/spectronaut/20250324_170201_MM_A1_SCMO001_directDIA_matching_default_set/MM_A1_SCMO001_directDIA_matching_default_set_Report_Protein_JB (Pivot).tsv",
                    header = T,
                    sep = "\t")
data1_cols <- colnames(data1)
data1_cols
summary(data1[,1:8])

# 20250325_111639_MM_A1_SCMO001_1cell_split_directDIA_default
data2 <- read.table("data/spectronaut/20250325_111639_MM_A1_SCMO001_1cell_split_directDIA_default/MM_A1_SCMO001_1cell_split_directDIA_default_Report_Protein_JB (Pivot).tsv",
                    header = T,
                    sep = "\t")
data2_cols <- colnames(data2)
data2_cols
summary(data2[,1:8])

# 20250325_111657_MM_A1_SCMO001_0cell_ctrls_directDIA_default
data3 <- read.table("data/spectronaut/20250325_111657_MM_A1_SCMO001_0cell_ctrls_directDIA_default/MM_A1_SCMO001_0cell_ctrls_directDIA_default_Report_Protein_JB (Pivot).tsv",
                    header = T,
                    sep = "\t")
data3_cols <- colnames(data3)
data3_cols
summary(data3[,1:8])

# 20250325_111740_MM_A1_SCMO001_1cell_split_+RRI_directDIA_default
data4 <- read.table("data/spectronaut/20250325_111740_MM_A1_SCMO001_1cell_split_+RRI_directDIA_default/MM_A1_SCMO001_1cell_split_+RRI_directDIA_default_Report_Protein_JB (Pivot).tsv",
                    header = T,
                    sep = "\t")
data4_cols <- colnames(data4)
data4_cols
summary(data4[,1:8])

# 20250325_111826_MM_A1_SCMO001_1cell_directDIA_default
data5 <- read.table("data/spectronaut/20250325_111826_MM_A1_SCMO001_1cell_directDIA_default/MM_A1_SCMO001_1cell_directDIA_default_Report_Protein_JB (Pivot).tsv",
                    header = T,
                    sep = "\t")
data5_cols <- colnames(data5)
data5_cols
summary(data5[,1:8])

# 20250325_111947_MM_A1_SCMO001_1cell_+RRI_directDIA_default
data6 <- read.table("data/spectronaut/20250325_111947_MM_A1_SCMO001_1cell_+RRI_directDIA_default/MM_A1_SCMO001_1cell_+RRI_directDIA_default_Report_Protein_JB (Pivot).tsv",
                    header = T,
                    sep = "\t")
data6_cols <- colnames(data6)
data6_cols
summary(data6[,1:8])

# 20250325_112124_MM_A1_SCMO001_10_100cells_split_noMBR_directDIA_default
data7 <- read.table("data/spectronaut/20250325_112124_MM_A1_SCMO001_10_100cells_split_noMBR_directDIA_default/MM_A1_SCMO001_10_100cells_split_noMBR_directDIA_default_Report_Protein_JB (Pivot).tsv",
                    header = T,
                    sep = "\t")
data7_cols <- colnames(data7)
data7_cols
summary(data7[,1:8])

# condition setup
overview <- read.table("data/spectronaut/MM_A1_SCMO001_directDIA_matching_default_set_ConditionSetup.tsv",
                       header = T,
                       sep = "\t",
                       comment.char = ">")

