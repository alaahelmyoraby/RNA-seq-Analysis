Here’s how you can convert your code into a **README.md** file for GitHub documentation. This file will provide an overview of the project, explain the steps, and include the code in a readable format.

---

# RNA-Seq Differential Expression Analysis Pipeline

This repository contains an RNA-Seq differential expression analysis pipeline using DESeq2. The pipeline processes raw count data, performs exploratory data analysis, normalizes the data, identifies differentially expressed genes (DEGs), and visualizes the results.

---

## Table of Contents
1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Steps](#steps)
   - [Load Libraries](#load-libraries)
   - [Load Data](#load-data)
   - [Exploratory Data Analysis](#exploratory-data-analysis)
   - [Gene ID Mapping](#gene-id-mapping)
   - [Metadata Preprocessing](#metadata-preprocessing)
   - [Expression Matrix Preprocessing](#expression-matrix-preprocessing)
   - [DESeq2 Analysis](#deseq2-analysis)
   - [Visualization](#visualization)
4. [Results](#results)
5. [License](#license)

---

## Overview
This pipeline processes RNA-Seq data to identify differentially expressed genes (DEGs) between two conditions (e.g., Control vs. Treatment). It includes:
- Exploratory data analysis (boxplots, histograms, PCA).
- Gene ID mapping (ENTREZ to gene symbols).
- Normalization using DESeq2.
- DEG identification and visualization (heatmaps, volcano plots, bar plots).

---

## Requirements
To run this pipeline, you need the following R packages:
```R
install.packages(c("readr", "org.Hs.eg.db", "vsn", "DESeq2", "dplyr", "ggplot2", "ComplexHeatmap", "circlize", "rgl"))
```

---

## Steps

### Load Libraries
```R
library(readr)
library(org.Hs.eg.db)
library("vsn")
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(rgl)
```

---

### Load Data
Load the raw count data and metadata:
```R
# Load raw counts
read_counts <- read.delim("C:/My pc/Egcombio/MODA/assignment/GSE275290_raw_counts.tsv")

# Load metadata
metadata <- read.csv("C:/My pc/Egcombio/MODA/assignment/Phenotable.csv")
```

---

### Exploratory Data Analysis
Perform initial exploratory analysis:
```R
# Boxplot
boxplot(log2(read_counts[,-1] + 1), main = "Exploratory Box plot", ylab = "log2_count", las = 2)

# Histogram
hist(as.matrix(log2(read_counts[,-1] + 1)), main = "Exploratory Histogram", xlab = "Sample", ylab = "log2_count", breaks = 50)
```

---

### Gene ID Mapping
Convert ENTREZ IDs to gene symbols:
```R
gene_ids <- as.character(read_counts$GeneID)
gene_symbol <- mapIds(
  org.Hs.eg.db, 
  keys = gene_ids, 
  keytype = "ENTREZID",
  column = "SYMBOL", 
  multiVals = "first"
)

# Create gene symbol data frame
gene_df <- data.frame(
  GeneID = names(gene_symbol),
  Gene_symbol = as.vector(gene_symbol),
  stringsAsFactors = FALSE
)

# Merge with raw counts
data <- merge(read_counts, gene_df, by = "GeneID", all.x = TRUE)
data <- data %>%
  dplyr::select(Gene_symbol, everything(), -GeneID)
```

---

### Metadata Preprocessing
Preprocess the metadata:
```R
meta <- metadata %>%
  select(Sample.Name, treatment) %>%
  rename(sampleid = Sample.Name, Condition = treatment)

meta <- meta %>%
  mutate(Condition = case_when(
    Condition == "FBZ for 48 hours" ~ "FBZ",
    Condition == "DMSO; Vehicle" ~ "Control",
    TRUE ~ Condition
  ))
```

---

### Expression Matrix Preprocessing
Preprocess the expression matrix:
```R
# Remove NA values
data <- na.omit(data)

# Aggregate by gene symbol
exp_data <- data %>% 
  dplyr::select(-Gene_symbol)

exp_data_agg <- aggregate(exp_data, by = list(data$Gene_symbol), FUN = mean)

# Set row names and remove Group.1 column
row.names(exp_data_agg) <- exp_data_agg$Group.1
exp_data_agg <- exp_data_agg[,-1]

# Remove genes with zero variance
varrow <- apply(exp_data_agg, 1, var, na.rm = TRUE)
cons_var <- (varrow == 0 | is.na(varrow))
exp_data_agg <- exp_data_agg[!cons_var,]

# Perform PCA
pca <- prcomp(t(log2(exp_data_agg + 1)), scale. = TRUE)
plot(pca$x[,1], pca$x[,2], main = "PCA of Filtered Raw Counts", xlab = "PC1", ylab = "PC2")
```

---

### DESeq2 Analysis
Run DESeq2 for differential expression analysis:
```R
# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = exp_data_agg,
  colData = meta,
  design = ~ Condition
)

# Filter low-count genes
dds <- dds[rowSums(counts(dds)) >= 10,]

# Run DESeq2
dds_run <- DESeq(dds)

# Extract results
res <- results(dds_run, contrast = c("Condition", "FBZ", "Control"), alpha = 0.05)
res <- res[complete.cases(res),]

# Add significance column
res_df <- as.data.frame(res)
res_df <- res_df %>%
  mutate(significance = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "Not Significant"))
```

---

### Visualization
Visualize the results:
```R
# Volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = -log(padj), color = significance)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  labs(
    title = "Volcano Plot",
    x = "log Fold change",
    y = "-log10 adjusted Pvalue"
  ) +
  geom_vline(xintercept = c(1, -1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")

# Heatmap of top 100 DEGs
deg_100 <- res_df[order(res_df$padj, -abs(res_df$log2FoldChange)),]
deg_100 <- deg_100[1:100,]
deg_100_exp <- vsn_norm[rownames(deg_100),]

Heatmap(
  matrix = deg_100_exp,
  top_annotation = column_annot,
  row_title = "Top 100 DEGs",
  column_title = "Samples",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 4)
)
```

---

## Results
The pipeline generates the following outputs:
- **Exploratory plots**: Boxplots, histograms, and PCA plots.
- **DEGs**: Lists of upregulated and downregulated genes.
- **Visualizations**: Volcano plots, heatmaps, and bar plots.

---

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

This **README.md** file provides a clear and structured overview of your RNA-Seq analysis pipeline. You can upload it to GitHub along with your code for documentation. Let me know if you need further assistance! 😊
