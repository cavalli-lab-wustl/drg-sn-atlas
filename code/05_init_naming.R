#' ======================================================================
#' 05_init_naming.R
#'
#' Assign 51 initial clusters to 7 major cell classes; remove batch-driven clusters.
#'
#' Inputs:  ../r_objects/combined_init_object.RDS
#' Outputs: ../r_objects/combined_object_init_named.RDS
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
###    CLUSTER NAMING & UMAPS EXPORT    ###
###                                     ###
###########################################
###########################################

combined_object <- readRDS("../r_objects/combined_init_object.RDS")
# Resolution 0.9 selected in 02_init_clustree.R
combined_object <- FindClusters(combined_object, resolution = 0.9)
DefaultAssay(combined_object) <- "RNA"

# Export highlight plots for each initial cluster (visual reference)
dir.create("../plots/umaps/init_cluster_highlights", recursive = TRUE, showWarnings = FALSE)
highlightAndSaveAllClusters(combined_object, save_directory = "../plots/umaps/init_cluster_highlights/")

# Create a copy for naming clusters
merged_named <- combined_object

# NOTE: decisions were made based on clear QC metrics, all contained in "combined_initial_cellcounts.xlsx"

# Map 51 initial clusters (res 0.9) to 7 major cell classes + BATCH
# Assignments based on marker gene expression in combined_initial_cellcounts.xlsx
cluster_to_celltype <- c(
  "0" = "glia",
  "1" = "glia",
  "2" = "BATCH",
  "3" = "glia",
  "4" = "glia",
  "5" = "glia",
  "6" = "fibro",
  "7" = "immune",
  "8" = "fibro",
  "9" = "vec",
  "10" = "glia",
  "11" = "immune",
  "12" = "glia",
  "13" = "fibro",
  "14" = "neuron",
  "15" = "mural",
  "16" = "vec",
  "17" = "mural",
  "18" = "fibro",
  "19" = "fibro",
  "20" = "immune",
  "21" = "glia",
  "22" = "neuron",
  "23" = "fibro",
  "24" = "fibro",
  "25" = "BATCH",
  "26" = "glia",
  "27" = "neuron",
  "28" = "immune",
  "29" = "fibro",
  "30" = "glia",
  "31" = "neuron",
  "32" = "neuron",
  "33" = "mitotic",
  "34" = "neuron",
  "35" = "immune",
  "36" = "neuron",
  "37" = "neuron",
  "38" = "fibro",
  "39" = "glia",
  "40" = "ery",
  "41" = "neuron",
  "42" = "neuron",
  "43" = "immune",
  "44" = "BATCH",
  "45" = "fibro",
  "46" = "vec",
  "47" = "mitotic",
  "48" = "BATCH",
  "49" = "fibro",
  "50" = "BATCH"
)

merged_named <- RenameIdents(merged_named, cluster_to_celltype)

# Remove batch-driven clusters (2, 25, 44, 48, 50) — driven by technical artifacts
merged_named <- subset(merged_named, idents = c("BATCH"), invert = TRUE)

# Re-run UMAP with background removed
merged_named <- RunUMAP(merged_named, reduction = "integrated.dr", dims = 1:50)

major_levels <- c(
  "glia",
  "neuron",
  "fibro",
  "immune",
  "vec",
  "mural",
  "ery",
  "mitotic"
)

Idents(merged_named) <- factor(Idents(merged_named), levels = major_levels)

# Define major cell types and their colors
major_cell_type_colors <- c(
  "neuron" = "#762a83",
  "glia" = "#0072B5",
  "fibro" = "#C1DB1B",
  "immune" = "#CC79A7",
  "mitotic" = "#660033",
  "vec" = "#D7191C",
  "mural" = "#EA6337",
  "ery" = "#FFCC00"
)

dir.create("../plots/umaps/combined_init", recursive = TRUE, showWarnings = FALSE)
p <- DimPlot(object = merged_named, reduction = "umap", cols = major_cell_type_colors, raster = FALSE) + NoLegend()
print(p)
ggsave(filename = "../plots/umaps/combined_init/combined_umap_init_major_celltypes.png", plot = p, device = "png", width = 10, height = 10, dpi = 300)


# Before major cell type annotation, mitotic cells need to be split (they are mixed)

# Subcluster mitotic cells separately to assign them to their true major types
mitotic_cells <- subset(merged_named, idents = "mitotic")
DefaultAssay(mitotic_cells) <- "SCT"
# resolution = 0.5: lower resolution sufficient for small mitotic subset
mitotic_cells <- FindClusters(mitotic_cells, resolution = 0.5, reduction = "integrated.dr")
mitotic_cells <- RunUMAP(mitotic_cells, reduction = "integrated.dr", dims = 1:50)
DefaultAssay(mitotic_cells) <- "RNA"
mitotic_cells <- NormalizeData(mitotic_cells)
mitotic_cells <- JoinLayers(mitotic_cells)

# Create a directory for saving output files
output_directory <- "../markers/combined_init_mitotic_cells"
dir.create(output_directory, showWarnings = FALSE)

# logfc.threshold = 0.25: default threshold (less strict than main clustering)
find_and_save_all_markers(mitotic_cells, output_directory, min.pct = 0.10, logfc.threshold = 0.25, summary = TRUE)

# Extract the average gene expression for each cluster
avg_exp <- AverageExpression(mitotic_cells)

# Restrict analysis to variable genes
mitotic_cells <- FindVariableFeatures(object = mitotic_cells)
variable_gene_names <- VariableFeatures(mitotic_cells, assay = "RNA")

avg_exp_var <- avg_exp$RNA[variable_gene_names, ]

# Convert sparse matrix to regular matrix
avg_exp_var_dense <- as.matrix(avg_exp_var)

# Compute the correlation matrix
cor_matrix <- cor(avg_exp_var_dense)

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
dir.create("../plots/correlations", recursive = TRUE, showWarnings = FALSE)
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
ggsave('../plots/correlations/mitotic_cells_seuratclusters_rna_correlation_matrix_hclust.pdf', plot=p, device='pdf', width=7, height=7, dpi=300)

init_abundances <- table(mitotic_cells@meta.data$study, Idents(mitotic_cells))
write.csv(init_abundances, "../metadata/mitotic_cells_initial_cellcounts_study.csv")

init_abundances <- table(mitotic_cells@meta.data$orig.ident, Idents(mitotic_cells))
write.csv(init_abundances, "../metadata/mitotic_cells_initial_cellcounts_batch.csv")

# Extract data
data <- FetchData(mitotic_cells, vars = c("nFeature_RNA", "nCount_RNA"))

# Add cluster information
data$cluster <- Idents(mitotic_cells)

# Calculate the average 'nFeature_RNA' and 'nCount_RNA' per cluster
average_values <- data %>%
  group_by(cluster) %>%
  summarise(
    Average_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
    Average_nCount_RNA = mean(nCount_RNA, na.rm = TRUE)
  )

# Export the result to a CSV file
write.csv(average_values, "../metadata/mitotic_cells_avg_nFeaturenCount.csv", row.names = FALSE)

mitotic_named <- mitotic_cells

# Assign mitotic subclusters to their true major cell types
cluster_to_celltype <- c(
  "0" = "glia",
  "1" = "glia",
  "2" = "immune",
  "3" = "glia",
  "4" = "fibro",
  "5" = "fibro",
  "6" = "glia",
  "7" = "glia",
  "8" = "immune",
  "9" = "BATCH"
)

mitotic_named <- RenameIdents(mitotic_named, cluster_to_celltype)

# Transfer mitotic cell labels back into the main object
mitotic_labels <- Idents(mitotic_named)
mitotic_cells <- Cells(mitotic_named)

Idents(merged_named, cells = mitotic_cells) <- mitotic_labels

# Remove batch-driven mitotic subcluster
merged_named <- subset(merged_named, idents = c("BATCH"), invert = TRUE)

saveRDS(merged_named, "../r_objects/combined_object_init_named.RDS")

message("05_init_naming complete.")
