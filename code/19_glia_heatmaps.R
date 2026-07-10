#' ======================================================================
#' 19_glia_heatmaps.R
#'
#' Condition-based gene expression heatmaps for glial subpopulations.
#'
#' Inputs:  ../r_objects/DRGSNMI_glia.RDS
#' Outputs: ../plots/glia_heatmap_*.pdf
#'
#' Author: Michael B. Thomsen, Ph.D.
#' Created: 2024-07-01
#' ======================================================================

# Load necessary libraries
library(Seurat)
library(pheatmap)

# 1. Define the desired conditions
desired_conditions <- c(
  "cavalli_aih",
  "cavalli_contra",
  "cavalli_drc",
  "cavalli_naive",
  "cavalli_sci",
  "cavalli_snc"
)

# 2. Subset the Seurat object to include only the desired conditions
glia_subset <- subset(glia, subset = Study %in% desired_conditions)

glia_subset@meta.data$celltype <- Idents(glia_subset)

# 3. Set the active identity to "Study" for grouping
Idents(glia_subset) <- "Study"

# 4. Specify your list of genes for differential expression
# Replace the example genes below with your actual gene list
gene_list <- c("Selenow", "Selenof", "Selenok","Selenom","Selenop","Selenot","Selenos","Selenoh", "Gfap", "Cd74")

# Optional: Verify that all genes in gene_list are present in the dataset
genes_present <- gene_list %in% rownames(glia_subset)
if(!all(genes_present)) {
  warning("Some genes in gene_list are not present in the dataset and will be excluded.")
  gene_list <- gene_list[genes_present]
}

# 5. Calculate the average expression of the specified genes across the conditions
average_expression <- AverageExpression(glia_subset, features = gene_list, group.by = "Study")$RNA

# 6. Scale the average expression for better visualization
scaled_avg_exp <- t(scale(t(average_expression)))

# 7. Generate the heatmap using pheatmap
pheatmap(
  scaled_avg_exp,
  cluster_rows = TRUE,       # Cluster genes
  cluster_cols = TRUE,       # Cluster conditions
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 10,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  border_color = NA,
  main = "Differential Gene Expression Heatmap",
  display_numbers = FALSE
)
