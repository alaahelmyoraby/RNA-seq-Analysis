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
4. [Further Analysis](#further-analysis)
5. [Results](#results)
6. [License](#license)

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
read_counts <- read.delim("C:/Egcombio/MODA/GSE275290_raw_counts.tsv")

# Load metadata
metadata <- read.csv("C:/Egcombio/MODA/Phenotable.csv")
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
sum(is.na(data))  # Check for NA values
data <- na.omit(data)  # Remove rows with NA values

dim(data)  # Check dimensions of the data after removing NAs

# Check for duplicated gene symbols
sum(duplicated(data$Gene_symbol))

# Aggregate by gene symbol
exp_data <- data %>% 
  dplyr::select(-Gene_symbol)

is.numeric(exp_data)  # Check if the data is numeric

exp_data_agg <- aggregate(exp_data, by = list(data$Gene_symbol), FUN = mean)

sum(duplicated(exp_data_agg$Group.1))  # Check for duplicated gene symbols after aggregation

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

dim(exp_data_agg)  # Check dimensions of the filtered data
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

# Normalize data using VST
ntd <- vst(dds)
exp.norm <- assay(ntd)

# Alternatively, use normTransform for normalization
ntd <- normTransform(dds)
exp.norm <- assay(ntd)

# Extract normalized expression levels of DEGs
degs <- res[res$padj < 0.05 & res$log2FoldChange > 1,]
degs.genes <- rownames(degs)
degs.exp <- vsn_norm[degs.genes,]
write.csv(degs.exp, "degs.exp.csv")
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

# Create column annotations for the heatmap
sam_condition <- meta$Condition
column_annot <- HeatmapAnnotation(
  Condition = sam_condition,
  col = list(Condition = c("Control" = "blue", "FBZ" = "red"))
)

# Generate the heatmap
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

# Perform PCA on DEGs
degs_normalized <- vsn_norm[rownames(vsn_norm) %in% rownames(degs),]
pca_result <- prcomp(degs_normalized, scale. = TRUE)
pca_data <- pca_result$x

# Plot 2D PCA
ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point() +
  theme_classic()

# Plot 3D PCA
mycolors <- ifelse(meta$Condition == "FBZ", "red", "lightgreen")
plot3d(pca_result$x[, 1:3], col = mycolors, size = 12, type = "s", main = "3D PCA Plot")

# Save normalized expression data for Enrichr
exp_norm <- read.csv("vsn_norm.csv", header = TRUE, row.names = 1)
write.table(exp_norm, file = "exp_norm.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# Convert DESeqResults to a data frame
res_df <- as.data.frame(res)

# Add a classification column
res_df <- res_df %>%
  mutate(classification = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "significantly upregulated",
    padj < 0.05 & log2FoldChange < -1 ~ "significantly downregulated",
    TRUE ~ "not significant"
  ))

# Extract top 10 upregulated and downregulated genes
top_genes <- res_df %>%
  filter(classification %in% c("significantly upregulated", "significantly downregulated")) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 10)

# Ensure the Gene column exists
top_genes <- top_genes %>%
  mutate(Gene = rownames(top_genes))  # Add Gene column if not already present

# Create plot for up/down regulated
ggplot(top_genes, aes(x = reorder(Gene, log2FoldChange), y = log2FoldChange, fill = classification)) +
  geom_bar(stat = "identity") +  # Bar plot
  coord_flip() +  # Flip coordinates for horizontal bars
  scale_fill_manual(values = c(
    "significantly upregulated" = "steelblue",  # Color for upregulated genes
    "significantly downregulated" = "coral"     # Color for downregulated genes
  )) +
  labs(
    title = "Top 10 Upregulated and Downregulated Genes",
    x = "Genes",
    y = "Log2 Fold Change",
    fill = "Classification"
  ) +
  theme_minimal() +  # Minimal theme
  theme(
    axis.text.y = element_text(size = 10),  # Adjust gene name text size
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")  # Center and style the title
  )
```
**VST (DESeq2)** and **VSN**:

| Feature                  | VST (DESeq2)                          | VSN                                |
|--------------------------|---------------------------------------|------------------------------------|
| **Designed for**         | RNA-seq count data                   | General high-throughput data      |
| **Distributional Assumption** | Negative binomial distribution     | Log-normal distribution           |
| **Variance Stabilization** | Parametric, based on mean-variance relationship | Non-parametric, glog transformation |
| **Normalization**        | Integrated into DESeq2 pipeline      | Standalone normalization method   |
| **Speed**                | Fast and efficient                   | Can be slower for large datasets  |
| **Use Case**             | RNA-seq analysis, visualization      | Multi-platform data integration   | 
---

## Further Analysis

### 5- DO Enrichment Analysis Using DAVID, GeneTrail, Enrichr
To perform Disease Ontology (DO) enrichment analysis, you can use the following platforms:
- **DAVID**: Upload your DEG list to [DAVID](https://david.ncifcrf.gov/) and select "Disease" or "Functional Annotation" for enrichment analysis.
- **GeneTrail**: Use [GeneTrail](https://genetrail.bioinf.uni-sb.de/) to analyze DEGs for disease associations.
- **Enrichr**: Upload your DEG list to [Enrichr](https://maayanlab.cloud/Enrichr/) and select "Disease" or "Drug" libraries for enrichment analysis.

---

### 6- Run L1000CDS2 Tool to Characterize Small Molecules
To identify small molecules that could reverse the disease signature, use the **L1000CDS2** tool:
1. Visit the [L1000CDS2](https://maayanlab.cloud/l1000cds2/) website.
2. Upload your DEG list (upregulated and downregulated genes).
3. Run the analysis to identify small molecules that reverse the disease signature.

---

## Results
The pipeline generates the following outputs:
- **Exploratory plots**: Boxplots, histograms, and PCA plots.
- **DEGs**: Lists of upregulated and downregulated genes.
- **Visualizations**: Volcano plots, heatmaps, and bar plots.

---
