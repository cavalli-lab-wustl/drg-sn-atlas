#' ======================================================================
#' 06_main_vec_clustering.R
#'
#' Subset, re-integrate, and cluster vascular endothelial populations at res=0.8.
#'
#' Inputs:  ../r_objects/combined_object_init_named.RDS
#' Outputs: ../r_objects/vecs_object_init.RDS, ../metadata/vecs_*_celltypes.rds
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


merged_named <- readRDS("../r_objects/combined_object_init_named.RDS")

vecs <- subset(merged_named, idents = c("vec"))

options(future.globals.maxSize = 3e+09)
vecs <- SCTransform(vecs)
vecs <- RunPCA(vecs, npcs = 50)
vecs <- IntegrateLayers(
  object = vecs,
  method = RPCAIntegration,
  normalization.method = "SCT",
  k.weight = 50  # reduced from default 100; smallest VEC batch has ~50 cells
)
vecs <- FindNeighbors(vecs, dims = 1:50, reduction = "integrated.dr")
vecs <- FindClusters(vecs, resolution = 0.5)
vecs <- RunUMAP(vecs, reduction = "integrated.dr", dims = 1:50)

# Plot UMAP, group by 'seurat_clusters'
dir.create("../plots/umaps/vecs", recursive = TRUE, showWarnings = FALSE)
p1 <- DimPlot(vecs, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, raster = F) + NoLegend()
print(p1)
ggsave("../plots/umaps/vecs/vecs_initial_umap.png", p1, width = 10, height = 10, units = "in", dpi = 300)

saveRDS(vecs, "../r_objects/vecs_object_init.RDS")

###########################################
###########################################
###                                     ###
###    CLUSTER RESOLUTION + CLUSTREE    ###
###                                     ###
###########################################
###########################################

vecs <- readRDS("../r_objects/vecs_object_init.RDS")

# Do not set RNA assay or normalize - FindClusters requires SCT assay

# Array of resolutions
dir.create("../plots/clustree/vecs", recursive = TRUE, showWarnings = FALSE)
resolutions <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

for (res in resolutions) {
  # Find clusters and create plots
  vecs <- FindClusters(vecs, resolution = res)
  p <- DimPlot(vecs, label = T, raster = F)

  # Save plots with and without legend
  ggplot2::ggsave(paste0("../plots/clustree/vecs/clustree_vecs_", res, "res.png"), plot = p, device = "png", width = 10, height = 10, dpi = 300)
}

# Run clustree to examine clustering layout
p <- clustree(vecs, prefix = "SCT_snn_res.") + theme(text = element_text(size = 20))

# Saving the plot with specified parameters
ggsave(filename = "../plots/clustree/vecs/clustree_vecs_summary_plot.png", plot = p, device = "png", width = 18, height = 9, dpi = 300)

# Resolution 0.8 selected via clustree

vecs <- readRDS("../r_objects/vecs_object_init.RDS")
vecs <- FindClusters(vecs, resolution = 0.8)

DefaultAssay(vecs) <- "RNA"
vecs <- NormalizeData(vecs)
vecs <- JoinLayers(vecs)

# Create a directory for saving output files
output_directory <- "../markers/vecs_init_full"
dir.create(output_directory, showWarnings = FALSE)
find_and_save_all_markers(vecs, output_directory, min.pct = 0.1, logfc.threshold = 0.25, summary = T)

init_abundances <- table(vecs@meta.data$study, Idents(vecs))
write.csv(init_abundances, "../metadata/vecs_initial_cellcounts.csv")

# Extract data
data <- FetchData(vecs, vars = c("nFeature_RNA", "nCount_RNA"))
data$cluster <- Idents(vecs)
average_values <- data %>%
  group_by(cluster) %>%
  summarise(
    Average_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
    Average_nCount_RNA = mean(nCount_RNA, na.rm = TRUE)
  )
# Export the result to a CSV file
write.csv(average_values, "../metadata/vecs_avg_nFeaturenCount.csv", row.names = FALSE)

dir.create("../plots/umaps/vecs_highlights", recursive = TRUE, showWarnings = FALSE)
highlightAndSaveAllClusters(vecs, save_directory = "../plots/umaps/vecs_highlights/")

# Extract the average gene expression for each cluster
avg_exp <- AverageExpression(vecs)

# Restrict analysis to variable genes
vecs <- FindVariableFeatures(object = vecs)
variable_gene_names <- VariableFeatures(vecs, assay = "RNA")

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

ggsave('../plots/correlations/vecs_seuratclusters_correlation_matrix_hclust.pdf', plot=p, device='pdf', width=12, height=12, dpi=300)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6309638/

# Lymphatic = Sox18, Prox1, Flt4 (Vegfr3)
# Venous = Ephb4, Smarca4 (Brg1), Nr2f2 (CoupTFII)

### Interactive marker exploration (run manually in RStudio) ###
# FeaturePlot(vecs, features = c("Ephb4", "Smarca4", "Nr2f2", "Flt4"), order = T)
#
# # Remove batch-driven clusters for clean dotplots
# vecs_plotting <- subset(vecs, idents = c(2, 7, 14), invert = T)
#
# DotPlot(vecs, features = kidney_artery_genes)
# DotPlot(vecs, features = kidney_vein_genes)
# DotPlot(vecs, features = kidney_capillary_genes)
# DotPlot(vecs, features = brain_artery_genes)
# DotPlot(vecs, features = brain_vein_genes)
# DotPlot(vecs, features = brain_capillary_genes)
# DotPlot(vecs, features = heart_artery_genes)
# DotPlot(vecs, features = heart_vein_genes)
# DotPlot(vecs, features = heart_capillary_genes)
#
# tissue_specific_vec_genes <- c(
#   "Fabp4",
#   "Aqp7",
#   "Scgb3a1",
#   "Nkx2-3",
#   "Oit3",
#   "Star",
#   "Pglyrp1",
#   "Lcn2",
#   "Igfbp5",
#   "Tmem100"
# )
#
# DotPlot(vecs_plotting, features = rev(tissue_specific_vec_genes)) + coord_flip()
# DotPlot(vecs_plotting, features = tissue_specific_vec_genes, group.by="tissue")

# Major-level naming
vecs_named <- vecs

# renaming all clusters
cluster_to_celltype <- c(
  "0" = "VEC A",
  "1" = "VEC A",
  "2" = "BATCH",
  "3" = "VEC A",
  "4" = "VEC V",
  "5" = "VEC A",
  "6" = "VEC A",
  "7" = "BATCH",
  "8" = "VEC V",
  "9" = "LEC",
  "10" = "VEC V",
  "11" = "BATCH",
  "12" = "VEC A",
  "13" = "VEC V",
  "14" = "BATCH",
  "15" = "VEC V",
  "16" = "E_CYC",
  "17" = "BATCH",
  "18" = "VEC A"
)

vecs_named <- RenameIdents(vecs_named, cluster_to_celltype)

major_levels <- c(
  "VEC A",
  "VEC V",
  "LEC",
  "E_CYC",
  "BATCH"
)

Idents(vecs_named) <- factor(Idents(vecs_named), levels = major_levels)

# Export idents (to be used in main object)
vecs_labels <- Idents(vecs_named)
vecs_cells <- Cells(vecs_named)

# Export minor cell type idents and cells
saveRDS(vecs_labels, "../metadata/vecs_major_celltypes.rds")
saveRDS(vecs_cells, "../metadata/vecs_major_cells.rds")

# REMOVE background, technical artifacts
vecs_named <- subset(vecs_named, idents = c("BATCH"), invert = TRUE)

# Quick UMAP rerun after doublet removal
vecs_named <- RunUMAP(vecs_named, reduction = "integrated.dr", dims = 1:50)

# https://pmc.ncbi.nlm.nih.gov/articles/PMC7137757/

### Interactive marker exploration (run manually in RStudio) ###
# FeaturePlot(vecs_named, features = c("Ephb2", "Hey1", "Gja5", "Dll4"), order = T)
# FeaturePlot(vecs_named, features = c("Nrp2", "Ephb4", "Nr2f2", "Slc38a5"), order = T)
# FeaturePlot(vecs_named, features = c("Gja5", "Nrp1", "Sox17", "Alk1"), order = T)
# FeaturePlot(vecs_named, features = c("Kdr", "Mfsd2a", "Car4", "Slc38a5"), order = T)
# FeaturePlot(vecs_named, features = c("Bgn", "Ctla2a", "Apoe", "Il6st"), order = T)
# FeaturePlot(vecs_named, features = c("Clu", "Crip1", "Fbln2", "Mecom"), order = T)
# FeaturePlot(vecs_named, features = c("Mgp", "Fbln5", "Vwf", "Rbp7"), order = T)
# FeaturePlot(vecs_named, features = c("Hey1", "Dll4", "Nr2f2", "Nrp2"), order = T)
# FeaturePlot(vecs_named, features = c("Kdr", "Mfsd2a", "Car4", "Slc38a5"), order = T)

# Plot UMAP
p1 <- DimPlot(vecs_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F) + NoLegend()
print(p1)
ggsave(filename = "../plots/umaps/vecs/vecs_umap_major_celltypes.png", plot = p1, device = "png", width = 5, height = 5, dpi = 300)

p2 <- DimPlot(vecs_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F, split.by = "tissue") + NoLegend()
print(p2)
ggsave(filename = "../plots/umaps/vecs/vecs_umap_major_celltypes_split.png", plot = p2, device = "png", width = 10, height = 5, dpi = 300)

p3 <- DimPlot(vecs_named, reduction = "umap", label = F, repel = TRUE, raster = F, group.by = "seurat_clusters") + NoLegend()
print(p3)
ggsave(filename = "../plots/umaps/vecs/vecs_umap_major_celltypes_unlabeled.png", plot = p3, device = "png", width = 5, height = 5, dpi = 300)

p4 <- DimPlot(vecs_named, reduction = "umap", label = F, repel = TRUE, raster = F, split.by = "tissue", group.by = "seurat_clusters") + NoLegend()
print(p4)
ggsave(filename = "../plots/umaps/vecs/vecs_umap_major_celltypes_split_unlabeled.png", plot = p4, device = "png", width = 10, height = 5, dpi = 300)

#### FIGURE S1 EXPORTS

dir.create("../plots/s1_umaps", recursive = TRUE, showWarnings = FALSE)
p_supp_major <- DimPlot(vecs_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F, split.by = "tissue") +
  labs(x = "UMAP-1", y = "UMAP-2") +
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  NoLegend()

print(p_supp_major)
ggsave(filename = "../plots/s1_umaps/vecs_umap_major_celltypes_split.png", plot = p_supp_major, device = "png", width = 7, height = 4, dpi = 300)

ggsave(filename = "../plots/s1_umaps/vecs_umap_major_celltypes_split.pdf", plot = p_supp_major, device = "pdf", width = 7, height = 4, dpi = 300)

# Minor-level naming
vecs_named <- vecs

# renaming all clusters
cluster_to_celltype <- c(
  "0" = "VEC A i",
  "1" = "VEC A ii",
  "2" = "BATCH",
  "3" = "VEC A iii",
  "4" = "VEC V i",
  "5" = "VEC A iv",
  "6" = "VEC LA",
  "7" = "BATCH",
  "8" = "VEC V ii",
  "9" = "LEC",
  "10" = "VEC V iii",
  "11" = "BATCH",
  "12" = "VEC A v",
  "13" = "VEC LV",
  "14" = "BATCH",
  "15" = "VEC V iv",
  "16" = "E_CYC",
  "17" = "BATCH",
  "18" = "VEC IFN"
)

vecs_named <- RenameIdents(vecs_named, cluster_to_celltype)

major_levels <- c(
  "VEC A i",
  "VEC A ii",
  "VEC A iii",
  "VEC A iv",
  "VEC A v",
  "VEC V i",
  "VEC V ii",
  "VEC V iii",
  "VEC V iv",
  "VEC LA",
  "VEC LV",
  "VEC IFN",
  "LEC",
  "E_CYC",
  "BATCH"
)

Idents(vecs_named) <- factor(Idents(vecs_named), levels = major_levels)

# Export idents (to be used in main object)
vecs_labels <- Idents(vecs_named)
vecs_cells <- Cells(vecs_named)

# Export minor cell type idents and cells
saveRDS(vecs_labels, "../metadata/vecs_minor_celltypes.rds")
saveRDS(vecs_cells, "../metadata/vecs_minor_cells.rds")

# REMOVE background, technical artifacts
vecs_named <- subset(vecs_named, idents = c("BATCH"), invert = TRUE)

# Quick UMAP rerun after doublet removal
vecs_named <- RunUMAP(vecs_named, reduction = "integrated.dr", dims = 1:50)

# Plot UMAP
p1 <- DimPlot(vecs_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F) + NoLegend()
print(p1)
ggsave(filename = "../plots/umaps/vecs/vecs_umap_minor_celltypes.png", plot = p1, device = "png", width = 5, height = 5, dpi = 300)

p2 <- DimPlot(vecs_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F, split.by = "tissue") + NoLegend()
print(p2)
ggsave(filename = "../plots/umaps/vecs/vecs_umap_minor_celltypes_split.png", plot = p2, device = "png", width = 10, height = 5, dpi = 300)

p3 <- DimPlot(vecs_named, reduction = "umap", label = F, repel = TRUE, raster = F, group.by = "seurat_clusters") + NoLegend()
print(p3)
ggsave(filename = "../plots/umaps/vecs/vecs_umap_minor_celltypes_unlabeled.png", plot = p3, device = "png", width = 5, height = 5, dpi = 300)

p4 <- DimPlot(vecs_named, reduction = "umap", label = F, repel = TRUE, raster = F, split.by = "tissue", group.by = "seurat_clusters") + NoLegend()
print(p4)
ggsave(filename = "../plots/umaps/vecs/vecs_umap_minor_celltypes_split_unlabeled.png", plot = p4, device = "png", width = 10, height = 5, dpi = 300)

#### FIGURE S1 EXPORTS

p_supp_minor <- DimPlot(vecs_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F, split.by = "tissue") +
  labs(x = "UMAP-1", y = "UMAP-2") +
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  NoLegend()

print(p_supp_minor)
ggsave(filename = "../plots/s1_umaps/vecs_umap_minor_celltypes_split.png", plot = p_supp_minor, device = "png", width = 7, height = 4, dpi = 300)

ggsave(filename = "../plots/s1_umaps/vecs_umap_minor_celltypes_split.pdf", plot = p_supp_minor, device = "pdf", width = 7, height = 4, dpi = 300)

message("06_main_vec_clustering complete.")
