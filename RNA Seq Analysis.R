# Load necessary libraries
library(readr)
library(org.Hs.eg.db)
library("vsn")
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(rgl)


# DESeq2 must deal with raw data
read_counts <- read.delim("C:/My pc/Egcombio/MODA/assignment/GSE275290_raw_counts.tsv")

# Exploratory Data Analysis: Box plot and Histogram
boxplot(log2(read_counts[,-1]+1), main="Exploratory Box plot", ylab="log2_count", las=2)
hist(as.matrix(log2(read_counts[,-1]+1)), main="Exploratory Histogram", xlab="Sample", ylab="log2_count", breaks=50)

# Mapping gene IDs to gene symbols
gene_ids <- as.character(read_counts$GeneID)  # Ensure gene IDs are characters

# Convert ENTREZ IDs to gene symbols using org.Hs.eg.db
gene_symbol <- mapIds(
  org.Hs.eg.db, keys = gene_ids, keytype = "ENTREZID",
  column = "SYMBOL", multiVals = "first"
)

# Create a data frame with gene IDs and corresponding gene symbols
gene_df <- data.frame(
  GeneID = names(gene_symbol),
  Gene_symbol = as.vector(gene_symbol),
  stringsAsFactors = FALSE
)

# Merge the gene symbols with the raw counts data
data <- merge(read_counts, gene_df, by="GeneID", all.x=TRUE)
data <- data %>%
  dplyr::select(Gene_symbol, everything(), -GeneID)

# Read metadata
metadata <- read.csv("C:/My pc/Egcombio/MODA/assignment/Phenotable.csv")

# Preprocess metadata
meta <- metadata %>%
  select(Sample.Name, treatment) %>%
  rename(sampleid = Sample.Name, Condition = treatment)

meta <- meta %>%
  mutate(Condition = case_when(
    Condition == "FBZ for 48 hours" ~ "FBZ",
    Condition == "DMSO; Vehicle" ~ "Control",
    TRUE ~ Condition
  ))

# Preprocessing expression matrix
sum(is.na(data))  # Check for NA values
data <- na.omit(data)  # Remove rows with NA values

dim(data)  # Check dimensions of the data after removing NAs

# Check for duplicated gene symbols
sum(duplicated(data$Gene_symbol))

# Aggregate expression data by gene symbol
exp_data <- data %>% 
  dplyr::select(-Gene_symbol)

is.numeric(exp_data)  # Check if the data is numeric

exp_data_agg <- aggregate(exp_data, by=list(data$Gene_symbol), FUN=mean)

sum(duplicated(exp_data_agg$Group.1))  # Check for duplicated gene symbols after aggregation

# Set row names and remove the Group.1 column
row.names(exp_data_agg) <- exp_data_agg$Group.1
exp_data_agg <- exp_data_agg[,-1]

# Remove genes with zero variance
varrow <- apply(exp_data_agg, 1, var, na.rm=TRUE)
cons_var <- (varrow == 0 | is.na(varrow))
exp_data_agg <- exp_data_agg[!cons_var,]

# Perform PCA on the filtered data
pca <- prcomp(t(log2(exp_data_agg + 1)), scale. = TRUE)
plot(pca$x[,1], pca$x[,2], main = "PCA of Filtered Raw Counts", xlab = "PC1", ylab = "PC2")

dim(exp_data_agg)  # Check dimensions of the filtered data

# DESeq2 analysis
exp <- exp_data_agg

# Ensure column names in expression data match metadata
all(colnames(exp) %in% meta$sampleid)
all(colnames(exp) == meta$sampleid)

# Reorder expression data to match metadata
exp <- exp[, meta$sampleid]

# Round expression values and factorize conditions
exp <- round(exp)
meta$Condition <- factor(meta$Condition, levels = c("Control", "FBZ"))

class(meta$Condition)  # Check the class of the Condition column

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = exp,
  colData = meta,
  design = ~ Condition
)

# Filter out low-count genes
dds <- dds[rowSums(counts(dds)) >= 10,]

# Run DESeq2 analysis
dds_run <- DESeq(dds)

# Extract results for the contrast FBZ vs Control
res <- results(dds_run, contrast=c("Condition", "FBZ", "Control"), alpha=0.05)
res <- res[complete.cases(res),]  # Remove rows with NA values

summary(res)  # Summarize the results

# Convert results to a data frame and add significance column
res_df <- as.data.frame(res)
res_df <- res_df %>%
  mutate(significance = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "Not Significant"))

# Normalize data using VST
row_count <- counts(dds, normalized=FALSE)
vsn_data <- vsn2(row_count)
vsn_norm <- exprs(vsn_data)
write.csv(vsn_norm, "vsn_norm.csv")

# Alternatively, use normTransform for normalization
ntd <- normTransform(dds)
exp.norm <- assay(ntd)

# Extract normalized expression levels of DEGs
degs <- res[res$padj < 0.05 & res$log2FoldChange > 1,]
degs.genes <- rownames(degs)
degs.exp <- vsn_norm[degs.genes,]
write.csv(degs.exp, "degs.exp.csv")

# Create a heatmap for the top 100 DEGs
deg_100 <- degs[order(degs$padj, -abs(degs$log2FoldChange)),]
deg_100 <- degs[1:100,]
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