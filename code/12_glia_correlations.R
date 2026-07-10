#' ======================================================================
#' 12_glia_correlations.R
#'
#' Pearson correlation of average variable-gene expression across glial
#' subpopulations (Minor-level idents). Produces an unordered heatmap and
#' a hierarchically-clustered heatmap.
#'
#' Inputs:  ../r_objects/DRGSNMI_glia.RDS
#' Outputs: ../plots/correlations/glia_named_rna_correlation_matrix.pdf
#'          ../plots/correlations/glia_named_rna_correlation_matrix_hclust.pdf
#'
#' Author: Michael B. Thomsen, Ph.D.
#' Created: 2024-07-01
#' ======================================================================


# Set working directory to script location (portable: Rscript or RStudio)
local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) setwd(dirname(normalizePath(sub("^--file=", "", f))))
  else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))
})

# Imports
source('utils/libraries.R')
source('utils/helper_functions.R')
source('utils/variables.R')


# ===== LOAD OBJECT =====
glia <- readRDS("../r_objects/DRGSNMI_glia.RDS")
DefaultAssay(glia) <- "RNA"


# ===== OUTPUT DIRECTORY =====
output_directory <- "../plots/correlations"
dir.create(output_directory, recursive = TRUE)


# ===== VARIABLE FEATURES + AVERAGE EXPRESSION =====
glia <- FindVariableFeatures(glia)
variable_gene_names <- VariableFeatures(glia, assay = "RNA")

# Average expression per ident, restricted to variable genes
avg_exp <- AverageExpression(glia, assays = "RNA", features = variable_gene_names)
avg_exp_var_dense <- as.matrix(avg_exp$RNA)

# Correlation matrix (cluster x cluster)
cor_matrix <- cor(avg_exp_var_dense)


# ===== UNORDERED HEATMAP =====
melted_cor <- reshape2::melt(cor_matrix)

p <- ggplot(melted_cor, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(-1, 1)) +
  geom_text(aes(label = sprintf("%.2f", value)), vjust = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  labs(fill = "Correlation", x = "Cluster", y = "Cluster") +
  coord_fixed()

print(p)

ggsave(
  file.path(output_directory, "glia_named_rna_correlation_matrix.pdf"),
  plot = p, device = "pdf", width = 5, height = 5, dpi = 300
)


# ===== HIERARCHICALLY CLUSTERED HEATMAP =====
hc <- hclust(as.dist(1 - cor_matrix), method = "average")
cluster_order <- order.dendrogram(as.dendrogram(hc))

cor_matrix_ordered <- cor_matrix[cluster_order, cluster_order]

melted_cor <- reshape2::melt(cor_matrix_ordered)
melted_cor$Var1 <- factor(melted_cor$Var1, levels = unique(melted_cor$Var1))
melted_cor$Var2 <- factor(melted_cor$Var2, levels = unique(melted_cor$Var2))

p <- ggplot(melted_cor, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(-1, 1)) +
  geom_text(aes(label = sprintf("%.2f", value)), vjust = 1, size = 7) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16)
  ) +
  labs(fill = "Correlation", x = "Cluster", y = "Cluster") +
  coord_fixed()

print(p)

ggsave(
  file.path(output_directory, "glia_named_rna_correlation_matrix_hclust.pdf"),
  plot = p, device = "pdf", width = 12, height = 12, dpi = 300
)

message("12_glia_correlations complete.")
