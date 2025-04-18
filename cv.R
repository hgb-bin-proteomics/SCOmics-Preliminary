library(tidyverse)
library(proDA)

get_cv <- function(x){
  return(sd(x, na.rm = T) / mean(x, na.rm = T))
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

na_drop_threshold = 2
normalize_proteomics = F
normalize_transcriptomics = F

colnames_rri <- c("cell e", "cell d", "cell c", "cell b", "cell a")

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
#data <- data[complete.cases(data),]
data <- data[!rowSums(is.na(data)) > na_drop_threshold,]

# remove contaminants
#data <- data |> filter(!grepl("cont_*", PG.ProteinGroups))

# median normalization
data[,3:7] <- if (normalize_proteomics) proDA::median_normalization(data.matrix(data[,3:7])) else data[,3:7]
data <- data[,3:7]
colnames(data) <- colnames_rri

data$cv <- apply(data, 1, get_cv)
data$mean_cv <- mean(data$cv)
data$sd_cv <- sd(data$cv)
sample_size_prot = nrow(data)
data$group = rep("Proteomics", sample_size_prot)

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
#data_ngs <- data_ngs[complete.cases(data_ngs),]
data_ngs <- data_ngs[!rowSums(is.na(data_ngs)) > na_drop_threshold,]

# median normalization
data_ngs[,2:6] <- if (normalize_transcriptomics) proDA::median_normalization(data.matrix(data_ngs[,2:6])) else data_ngs[,2:6]
data_ngs <- data_ngs[,2:6]
colnames(data_ngs) <- colnames_rri

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
