# RNA-seq downstream analysis
# Author: Meyha Bishnoi

# Install required packages 
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("DESeq2", "ggplot2", "ggpubr", "pheatmap", "RColorBrewer", "clusterProfiler", "enrichplot", "org.Mm.eg.db", "BiocParallel"))

# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(BiocParallel)

# Register parallel processing
register(MulticoreParam(4))  # Adjust the number of cores as needed

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
dds <- DESeq(dds, parallel = TRUE)  # Enable parallel processing

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
         main = "Sample-to-Sample Distance Heatmap",
         annotation_col = data_col["Group"])  # Add group annotations

# Step 7: Differential Expression Analysis (Lung & Blood Case vs Control)

# 1. Lung Case vs Control
comparison_name_lung <- "Lung case v control"
res_lung <- results(dds, contrast = c("Group", "Lung_WT_Case", "Lung_WT_Control"), parallel = TRUE)
res_lung <- res_lung[order(res_lung$padj, na.last = NA), ]

# 2. Blood Case vs Control
comparison_name_blood <- "Blood case v control"
res_blood <- results(dds, contrast = c("Group", "Blood_WT_Case", "Blood_WT_Control"), parallel = TRUE)
res_blood <- res_blood[order(res_blood$padj, na.last = NA), ]

# Save significant genes
sig_genes_lung <- subset(res_lung, padj < 0.05)
sig_genes_blood <- subset(res_blood, padj < 0.05)

# Add Regulation Column (Upregulated or Downregulated based on log2FoldChange)
sig_genes_lung$Regulation <- ifelse(sig_genes_lung$log2FoldChange > 0, "Upregulated", "Downregulated")
sig_genes_blood$Regulation <- ifelse(sig_genes_blood$log2FoldChange > 0, "Upregulated", "Downregulated")

# Save significant genes with regulation information
write.csv(as.data.frame(sig_genes_lung), "/Users/meyha/Desktop/Lung_Significant_Genes_With_Regulation.csv")
write.csv(as.data.frame(sig_genes_blood), "/Users/meyha/Desktop/Blood_Significant_Genes_With_Regulation.csv")

# Count the number of Upregulated and Downregulated Genes
lung_regulation_summary <- table(sig_genes_lung$Regulation)
blood_regulation_summary <- table(sig_genes_blood$Regulation)

# Print summary to console
print("Lung Regulation Summary:")
print(lung_regulation_summary)

print("Blood Regulation Summary:")
print(blood_regulation_summary)

# Step 7.1: Filter for Genes of Interest
genes_of_interest <- c("Ifit1", "Ifit3", "Oas1a", "Oas2", "Oas3", "Mx1", "Irf7", "Irf9", "Stat2", "Rsad2", "Isg15",
                       "Ifng", "Irf1", "Irf8", "Gbp2", "Gbp5", "Gbp7", "Tap1", "Tap2", "Nos2",
                       "III7a", "III7f", "Cxcl9", "Cxcl10", "Pf4", "Mpo",
                       "Unc93b1", "Slamf7", "Trim30a", "Trim30d", "Pml", "Sp100", "Dhx58", "Phf11d", "Usp18", "Fcgr1", "Siglec1", "Aif1")

# Map gene symbols to ENSEMBL IDs
genes_of_interest_ensembl <- mapIds(
  org.Mm.eg.db,
  keys = genes_of_interest,
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first"
)

# Remove any NA values (if some symbols don't have ENSEMBL mappings)
genes_of_interest_ensembl <- na.omit(genes_of_interest_ensembl)

# Filter Lung and Blood DE genes
filtered_lung_genes <- sig_genes_lung[rownames(sig_genes_lung) %in% genes_of_interest_ensembl, ]
filtered_blood_genes <- sig_genes_blood[rownames(sig_genes_blood) %in% genes_of_interest_ensembl, ]

# Save filtered genes
write.csv(filtered_lung_genes, "/Users/meyha/Desktop/Lung_Filtered_Genes_of_Interest.csv")
write.csv(filtered_blood_genes, "/Users/meyha/Desktop/Blood_Filtered_Genes_of_Interest.csv")

# Debugging: Check overlaps
overlap_lung <- intersect(rownames(sig_genes_lung), genes_of_interest_ensembl)
overlap_blood <- intersect(rownames(sig_genes_blood), genes_of_interest_ensembl)

print("Genes of interest found in lung DE results:")
print(overlap_lung)

print("Genes of interest found in blood DE results:")
print(overlap_blood)

# Step 8: Volcano Plot for Lung and Blood
generate_volcano_plot <- function(res, comparison_name, filename) {
  volcano_data <- as.data.frame(res)
  volcano_data$Significance <- ifelse(volcano_data$padj < 0.05, "Significant", "Not Significant")
  volcano_data$GeneName <- rownames(volcano_data)  # Extract gene names
  
  # Highlight genes of interest
  volcano_data$Highlight <- ifelse(volcano_data$GeneName %in% genes_of_interest_ensembl, "Yes", "No")
  
  ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance, shape = Highlight)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("grey", "red")) +
    scale_shape_manual(values = c(16, 17)) +  # Different shapes for highlighted genes
    theme_minimal() +
    ggtitle(paste("Volcano Plot -", comparison_name)) +
    theme(plot.title = element_text(hjust = 0.5, size = 16), legend.position = "right")
  
  ggsave(filename, width = 8, height = 6, dpi = 300)
}

# Generate & Save Volcano Plots
generate_volcano_plot(res_lung, comparison_name_lung, "Volcano_Lung.png")
generate_volcano_plot(res_blood, comparison_name_blood, "Volcano_Blood.png")

# Step 9: GO Enrichment Analysis (BP, MF, CC)
run_go_analysis <- function(sig_gene_ids, ont) {
  enrichGO(gene = sig_gene_ids,
           universe = rownames(dds),
           OrgDb = org.Mm.eg.db,
           keyType = "ENSEMBL",
           ont = ont,
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05)
}

# Get Significant Gene IDs
sig_gene_ids_lung <- rownames(sig_genes_lung)
sig_gene_ids_blood <- rownames(sig_genes_blood)

# Run GO Analysis for Lung
go_bp_lung <- run_go_analysis(sig_gene_ids_lung, "BP")
go_mf_lung <- run_go_analysis(sig_gene_ids_lung, "MF")
go_cc_lung <- run_go_analysis(sig_gene_ids_lung, "CC")

# Run GO Analysis for Blood
go_bp_blood <- run_go_analysis(sig_gene_ids_blood, "BP")
go_mf_blood <- run_go_analysis(sig_gene_ids_blood, "MF")
go_cc_blood <- run_go_analysis(sig_gene_ids_blood, "CC")

# Save GO Results
write.csv(go_bp_lung@result, "/Users/meyha/Desktop/GO_BP_Lung.csv")
write.csv(go_mf_lung@result, "/Users/meyha/Desktop/GO_MF_Lung.csv")
write.csv(go_cc_lung@result, "/Users/meyha/Desktop/GO_CC_Lung.csv")

write.csv(go_bp_blood@result, "/Users/meyha/Desktop/GO_BP_Blood.csv")
write.csv(go_mf_blood@result, "/Users/meyha/Desktop/GO_MF_Blood.csv")
write.csv(go_cc_blood@result, "/Users/meyha/Desktop/GO_CC_Blood.csv")

# Step 10: GO Dot Plot
create_go_dotplot <- function(go_bp, go_mf, go_cc, filename) {
  go_bp@result$Category <- "Biological Process"
  go_mf@result$Category <- "Molecular Function"
  go_cc@result$Category <- "Cellular Component"
  
  # Merge results
  go_combined <- dplyr::bind_rows(go_bp@result, go_mf@result, go_cc@result)
  
  # Select top 10 terms per category
  go_combined_top <- go_combined %>%
    dplyr::group_by(Category) %>%
    dplyr::top_n(-10, p.adjust)
  
  # Create Dot Plot
  ggplot(go_combined_top, aes(x = Category, y = Description, size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    theme_minimal() +
    labs(title = "GO Enrichment Analysis (BP, MF, CC)",
         x = "GO Category", y = "GO Term", size = "Gene Count", color = "Adjusted p-value") +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave(filename, width = 10, height = 6, dpi = 300)
}

# Generate Dot Plots
create_go_dotplot(go_bp_lung, go_mf_lung, go_cc_lung, "GO_DotPlot_Lung.png")
create_go_dotplot(go_bp_blood, go_mf_blood, go_cc_blood, "GO_DotPlot_Blood.png")

# Step 11: Simplify GO Enrichment Results
simplified_ego_lung <- simplify(go_bp_lung, cutoff = 0.7, by = "p.adjust", select_fun = min)
simplified_ego_blood <- simplify(go_bp_blood, cutoff = 0.7, by = "p.adjust", select_fun = min)

# Save Simplified Results
write.csv(as.data.frame(simplified_ego_lung), "/Users/meyha/Desktop/Simplified_GO_Lung.csv")
write.csv(as.data.frame(simplified_ego_blood), "/Users/meyha/Desktop/Simplified_GO_Blood.csv")

# List of all ENSMUSG gene IDs
all_ensmusg_ids <- c("ENSMUSG00000105504", "ENSMUSG00000028270", "ENSMUSG00000027514", 
                     "ENSMUSG00000027639", "ENSMUSG00000028268", "ENSMUSG00000059089", 
                     "ENSMUSG00000028793", "ENSMUSG00000079363", "ENSMUSG00000015947", 
                     "ENSMUSG00000059498", "ENSMUSG00000028494", "ENSMUSG00000038507", 
                     "ENSMUSG00000023206", "ENSMUSG00000039753", "ENSMUSG00000002103", 
                     "ENSMUSG00000028792")

# Mapping ENSMUSG IDs to gene symbols
gene_symbols <- mapIds(org.Mm.eg.db, keys = all_ensmusg_ids, 
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Mappping ENSMUSG IDs to gene descriptions
gene_descriptions <- mapIds(org.Mm.eg.db, keys = all_ensmusg_ids, 
                            column = "GENENAME", keytype = "ENSEMBL", multiVals = "first")

# Create a dataframe
all_gene_info <- data.frame(ENSEMBL = all_ensmusg_ids, 
                            Gene_Symbol = gene_symbols, 
                            Description = gene_descriptions)

# Save results and view output
write.csv(all_gene_info, "All_Gene_Annotations.csv", row.names = FALSE)

print(all_gene_info)




