#' ======================================================================
#' 03_init_correlations.R
#'
#' Compute Pearson correlation matrix across initial clusters on variable genes.
#'
#' Inputs:  ../r_objects/combined_init_object.RDS
#' Outputs: ../plots/correlations/combined_init_rna_correlation_hclust.pdf
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


###########################################
###########################################
###                                     ###
###        COMBINED CORRELATIONS        ###
###                                     ###
###########################################
###########################################

combined_object <- readRDS("../r_objects/combined_init_object.RDS")
dir.create("../plots/correlations", recursive = TRUE, showWarnings = FALSE)

# Resolution 0.9 selected in 02_init_clustree.R
combined_object <- FindClusters(combined_object, resolution = 0.9)
DefaultAssay(combined_object) <- "RNA"
combined_object <- NormalizeData(combined_object)
combined_object <- JoinLayers(combined_object)

# Extract the average gene expression for each cluster
avg_exp <- AverageExpression(combined_object)

# Restrict analysis to variable genes
combined_object <- FindVariableFeatures(object = combined_object)
variable_gene_names <- VariableFeatures(combined_object, assay = "RNA")

avg_exp_var <- avg_exp$RNA[variable_gene_names, ]

# Convert sparse matrix to regular matrix
avg_exp_var_dense <- as.matrix(avg_exp_var)

# Compute the correlation matrix
cor_matrix <- cor(avg_exp_var_dense)

# Melt the matrix for ggplot2 plotting
# Use the melt function from reshape2 explicitly
melted_cor <- reshape2::melt(cor_matrix)

# Hierarchical clustering
hc <- hclust(as.dist(1 - cor_matrix), method="average")  # 1 minus to transform correlation to distance
cluster_order <- as.dendrogram(hc) %>% order.dendrogram()

# Reorder the correlation matrix based on the hierarchical clustering
cor_matrix_ordered <- cor_matrix[cluster_order, cluster_order]

# Melt the matrix for ggplot2 plotting
melted_cor <- reshape2::melt(cor_matrix_ordered)
melted_cor$Var1 <- factor(melted_cor$Var1, levels = unique(melted_cor$Var1))
melted_cor$Var2 <- factor(melted_cor$Var2, levels = unique(melted_cor$Var2))

# Plot
p <- ggplot(melted_cor, aes(Var1, Var2, fill=value)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(-1, 1)) +
  geom_text(aes(label=sprintf("%.2f", value)), vjust=1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust=1, size=12, hjust=1),
        axis.text.y = element_text(size=12)) +
  labs(fill="Correlation", x="Cluster", y="Cluster") +
  coord_fixed()

print(p)

ggsave('../plots/correlations/combined_init_rna_correlation_hclust.pdf', plot=p, device='pdf', width=20, height=20, dpi=300)

message("03_init_correlations complete.")
