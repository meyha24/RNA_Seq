# RNA-seq downstream analysis
# Author: Meyha Bishnoi



# Install required packages 
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("DESeq2", "ggplot2", "ggpubr", "pheatmap", "RColorBrewer", "clusterProfiler", "enrichplot", "org.Mm.eg.db"))

# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

# Step 1: Load and Clean Data
# Load the counts table
data_count <- read.delim("/Users/meyha/Desktop/counts_new.txt", header = TRUE, row.names = 1)

# Clean up column names to retain only sample IDs
colnames(data_count) <- sub(".*SRR", "SRR", colnames(data_count))

# Remove unnecessary columns like "Chr", "Start", "End", "Strand", and "Length"
data_count <- data_count[, grep("SRR", colnames(data_count))]

# Load the sample metadata
data_col <- read.csv("/Users/meyha/Desktop/Samples.csv", header = TRUE, row.names = 1)

# Ensure the order of samples matches between counts and metadata
# Identify mismatches
missing_in_counts <- setdiff(rownames(data_col), colnames(data_count))
missing_in_metadata <- setdiff(colnames(data_count), rownames(data_col))

if (length(missing_in_counts) > 0 || length(missing_in_metadata) > 0) {
  stop("Mismatch between data_count and data_col: Check sample names.")
}

# Remove the suffix "_20250124_sorted.bam" from column names
colnames(data_count) <- gsub("_20250124_sorted.bam", "", colnames(data_count))

# Re-check for mismatches
print("Re-checking mismatches after cleaning column names...")
print(setdiff(rownames(data_col), colnames(data_count)))  # Should return character(0)
print(setdiff(colnames(data_count), rownames(data_col)))  # Should return character(0)

# Get common sample names
common_samples <- intersect(rownames(data_col), colnames(data_count))

# Subset data_count and data_col
data_count <- data_count[, common_samples, drop = FALSE]
data_col <- data_col[common_samples, , drop = FALSE]

# Check for NA values in the count matrix
if (any(is.na(data_count))) {
  print("NA values found in the count matrix. Replacing NA values with 0...")
  data_count[is.na(data_count)] <- 0
}

# Remove rows with all zero counts
data_count <- data_count[rowSums(data_count) > 0, ]

# Step 2: Create DESeqDataSet
data_col$Group <- as.factor(data_col$Group)
dds <- DESeqDataSetFromMatrix(countData = data_count, colData = data_col, design = ~ Group)

# Step 3: Run DESeq
dds <- DESeq(dds)

# Step 4: Variance Stabilizing Transformation (VST)
vst_data <- vst(dds, blind = TRUE)

# Step 5: PCA Plot
pca_plot <- plotPCA(vst_data, intgroup = "Group") +
  theme_minimal(base_size = 14) +
  ggtitle("PCA Plot of Gene Expression") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "right",
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold")
  ) +
  scale_color_manual(
    name = "Group",
    values = c(
      "Blood_WT_Case" = "#D72638",
      "Blood_WT_Control" = "#1B998B",
      "Lung_WT_Case" = "#F46036",
      "Lung_WT_Control" = "#1E3888"
    )
  )
ggsave("Improved_PCA_Plot_Custom_Colors.png", plot = pca_plot, width = 8, height = 6, dpi = 300)

# Step 6: Heatmap of Sample-to-Sample Distances
sample_dist <- dist(t(assay(vst_data)))
sample_dist_matrix <- as.matrix(sample_dist)
rownames(sample_dist_matrix) <- colnames(vst_data)
colnames(sample_dist_matrix) <- colnames(vst_data)
heatmap_colors <- colorRampPalette(rev(brewer.pal(9, "Reds")))(255)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         col = heatmap_colors,
         main = "Sample-to-Sample Distance Heatmap")

# Step 7: Differential Expression Analysis
res <- results(dds, contrast = c("Group", "Lung_WT_Control", "Blood_WT_Case"))
res <- res[order(res$padj, na.last = NA), ]
sig_genes <- res[res$padj < 0.05 & !is.na(res$padj), ]
write.csv(as.data.frame(sig_genes), "/Users/meyha/Desktop/significant_genes.csv")

# Step 8: Volcano Plot
volcano_data <- as.data.frame(res)
volcano_data$Significance <- ifelse(volcano_data$padj < 0.05, "Significant", "Not Significant")

ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  ggtitle("Volcano Plot of DE Genes") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Center and adjust title size
    legend.position = "right"  # Keep the legend to the right
  )

# Step 9: GO Enrichment Analysis
sig_gene_ids <- rownames(sig_genes)
ego <- enrichGO(
  gene = sig_gene_ids,
  universe = rownames(res),
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
write.csv(as.data.frame(ego), "/Users/meyha/Desktop/GO_enrichment_results.csv")

# Simplify enrichment results
simplified_ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
write.csv(as.data.frame(simplified_ego), "/Users/meyha/Desktop/simplified_GO_results.csv")
