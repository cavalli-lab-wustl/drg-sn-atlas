#' ======================================================================
#' 06_main_mural_clustering.R
#'
#' Subset, re-integrate, and cluster mural cell populations at res=0.8.
#'
#' Inputs:  ../r_objects/combined_object_init_named.RDS
#' Outputs: ../r_objects/mural_object_init.RDS, ../metadata/mural_*_celltypes.rds
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

mural <- subset(merged_named, idents = c("mural"))

options(future.globals.maxSize = 3e+09)
mural <- SCTransform(mural)
mural <- RunPCA(mural, npcs = 50)
mural <- IntegrateLayers(
  object = mural,
  method = RPCAIntegration,
  normalization.method = "SCT",
  k.weight = 50  # reduced from default 100; smallest mural batch has ~50 cells
)
mural <- FindNeighbors(mural, dims = 1:50, reduction = "integrated.dr")
mural <- FindClusters(mural, resolution = 0.5)
mural <- RunUMAP(mural, reduction = "integrated.dr", dims = 1:50)

# Plot UMAP, group by 'seurat_clusters'
dir.create("../plots/umaps/mural", recursive = TRUE, showWarnings = FALSE)
p1 <- DimPlot(mural, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, raster = F) + NoLegend()
print(p1)
ggsave("../plots/umaps/mural/mural_initial_umap.png", p1, width = 10, height = 10, units = "in", dpi = 300)

saveRDS(mural, "../r_objects/mural_object_init.RDS")


###########################################
###########################################
###                                     ###
###    CLUSTER RESOLUTION + CLUSTREE    ###
###                                     ###
###########################################
###########################################

mural <- readRDS("../r_objects/mural_object_init.RDS")

# Do not set RNA assay or normalize - FindClusters requires SCT assay

# Array of resolutions
dir.create("../plots/clustree/mural", recursive = TRUE, showWarnings = FALSE)
resolutions <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

for (res in resolutions) {
  # Find clusters and create plots
  mural <- FindClusters(mural, resolution = res)
  p <- DimPlot(mural, label = T, raster = F)

  # Save plots with and without legend
  ggplot2::ggsave(paste0("../plots/clustree/mural/clustree_mural_", res, "res.png"), plot = p, device = "png", width = 10, height = 10, dpi = 300)
}

# Run clustree to examine clustering layout
p <- clustree(mural, prefix = "SCT_snn_res.") + theme(text = element_text(size = 20))

# Saving the plot with specified parameters
ggsave(filename = "../plots/clustree/mural/clustree_mural_summary_plot.png", plot = p, device = "png", width = 18, height = 9, dpi = 300)

# Resolution 0.8 selected via clustree

mural <- readRDS("../r_objects/mural_object_init.RDS")
mural <- FindClusters(mural, resolution = 0.8)

DefaultAssay(mural) <- "RNA"
mural <- NormalizeData(mural)
mural <- JoinLayers(mural)

# Create a directory for saving output files
output_directory <- "../markers/mural_init_full"
dir.create(output_directory, showWarnings = FALSE)
find_and_save_all_markers(mural, output_directory, min.pct = 0.1, logfc.threshold = 0.25, summary = T)


init_abundances <- table(mural@meta.data$study, Idents(mural))
write.csv(init_abundances, "../metadata/mural_initial_cellcounts.csv")


# Extract data
data <- FetchData(mural, vars = c("nFeature_RNA", "nCount_RNA"))
data$cluster <- Idents(mural)
average_values <- data %>%
  group_by(cluster) %>%
  summarise(
    Average_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
    Average_nCount_RNA = mean(nCount_RNA, na.rm = TRUE)
  )
# Export the result to a CSV file
write.csv(average_values, "../metadata/mural_avg_nFeaturenCount.csv", row.names = FALSE)

dir.create("../plots/umaps/mural_highlights", recursive = TRUE, showWarnings = FALSE)
highlightAndSaveAllClusters(mural, save_directory = "../plots/umaps/mural_highlights/")


# Extract the average gene expression for each cluster
avg_exp <- AverageExpression(mural)

# Restrict analysis to variable genes
mural <- FindVariableFeatures(object = mural)
variable_gene_names <- VariableFeatures(mural, assay = "RNA")

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

ggsave('../plots/correlations/mural_seuratclusters_correlation_matrix_hclust.pdf', plot=p, device='pdf', width=12, height=12, dpi=300)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6309638/

# Lymphatic = Sox18, Prox1, Flt4 (Vegfr3)
# Venous = Ephb4, Smarca4 (Brg1), Nr2f2 (CoupTFII)

### Interactive marker exploration (run manually in RStudio) ###
# FeaturePlot(mural, features = c("Ephb4", "Smarca4", "Nr2f2", "Flt4"), order = T)
# FeaturePlot(mural, features = c("Rgs5", "Ace2", "Acta2", "Flt1"), order = T)
# FeaturePlot(mural, features = c("Mki67", "Ace2", "Acta2", "Top2a"), order = T)
# FeaturePlot(mural, features = c("Pdgfra", "Dcn", "Lum", "Smoc2"), order = T)
# FeaturePlot(mural, features = c("Igfbp2", "Lrrc15", "Wnt11", "Actg2"), order = T)
#
# # Remove low-quality/batch clusters for clean dotplots
# mural_plotting <- subset(mural, idents = c(1, 3, 8, 10, 13), invert = T)
#
# DotPlot(mural, features = kidney_artery_genes) #5
# DotPlot(mural, features = kidney_vein_genes) #4/7
# DotPlot(mural, features = kidney_capillary_genes) #6/7
# DotPlot(mural, features = brain_artery_genes)
# DotPlot(mural, features = brain_vein_genes)
# DotPlot(mural, features = brain_capillary_genes)
# DotPlot(mural, features = heart_artery_genes)
# DotPlot(mural, features = heart_vein_genes)
# DotPlot(mural, features = heart_capillary_genes)
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
# DotPlot(mural_plotting, features = rev(tissue_specific_vec_genes)) + coord_flip()
# DotPlot(mural_plotting, features = tissue_specific_vec_genes, group.by="Tissue")

mural_named <- mural

# Major-level naming
# renaming all clusters
cluster_to_celltype <- c(
  "0" = "PERI",
  "1" = "PERI",
  "2" = "PERI",
  "3" = "VSMC",
  "4" = "VSMC",
  "5" = "VSMC",
  "6" = "VSMC",
  "7" = "VSMC",
  "8" = "VSMC",
  "9" = "BATCH",
  "10" = "PERI",
  "11" = "M_CYC",
  "12" = "BATCH",
  "13" = "BATCH",
  "14" = "PERI",
  "15" = "PERI",
  "16" = "PERI"
)

mural_named <- RenameIdents(mural_named, cluster_to_celltype)

major_levels <- c(
  "PERI",
  "VSMC",
  "M_CYC",
  "BATCH"
)

Idents(mural_named) <- factor(Idents(mural_named), levels = major_levels)


# Export idents (to be used in main object)
mural_labels <- Idents(mural_named)
mural_cells <- Cells(mural_named)


# REMOVE background, technical artifacts
mural_named <- subset(mural_named, idents = c("BATCH"), invert = TRUE)

# Quick UMAP rerun after doublet removal
mural_named <- RunUMAP(mural_named, reduction = "integrated.dr", dims = 1:50)

# Export major cell type idents and cells
saveRDS(mural_labels, "../metadata/mural_major_celltypes.rds")
saveRDS(mural_cells, "../metadata/mural_major_cells.rds")


#### FIGURE S1 EXPORTS

dir.create("../plots/s1_umaps", recursive = TRUE, showWarnings = FALSE)
p_supp_major <- DimPlot(mural_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F, split.by = "tissue") +
  labs(x = "UMAP-1", y = "UMAP-2") +
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  NoLegend()


print(p_supp_major)
ggsave(filename = "../plots/s1_umaps/mural_umap_major_celltypes_split.png", plot = p_supp_major, device = "png", width = 7, height = 4, dpi = 300)

ggsave(filename = "../plots/s1_umaps/mural_umap_major_celltypes_split.pdf", plot = p_supp_major, device = "pdf", width = 7, height = 4, dpi = 300)

mural_named <- mural

# Minor-level naming
# renaming all clusters
cluster_to_celltype <- c(
  "0" = "PERI I",
  "1" = "PERI I",
  "2" = "PERI I",
  "3" = "VSMC I",
  "4" = "VSMC I",
  "5" = "VSMC II",
  "6" = "VSMC I",
  "7" = "VSMC I",
  "8" = "VSMC III",
  "9" = "BATCH",
  "10" = "PERI I",
  "11" = "M_CYC",
  "12" = "BATCH",
  "13" = "BATCH",
  "14" = "PERI II",
  "15" = "PERI IFN",
  "16" = "PERI I"
)

mural_named <- RenameIdents(mural_named, cluster_to_celltype)

major_levels <- c(
  "PERI I",
  "PERI II",
  "PERI IFN",
  "VSMC I",
  "VSMC II",
  "VSMC III",
  "M_CYC",
  "BATCH"
)

Idents(mural_named) <- factor(Idents(mural_named), levels = major_levels)


# Export idents (to be used in main object)
mural_labels <- Idents(mural_named)
mural_cells <- Cells(mural_named)


# Export minor cell type idents and cells
saveRDS(mural_labels, "../metadata/mural_minor_celltypes.rds")
saveRDS(mural_cells, "../metadata/mural_minor_cells.rds")

# REMOVE background, technical artifacts
mural_named <- subset(mural_named, idents = c("BATCH"), invert = TRUE)

# Quick UMAP rerun after doublet removal
mural_named <- RunUMAP(mural_named, reduction = "integrated.dr", dims = 1:50)

# Plot UMAP
p1 <- DimPlot(mural_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F) + NoLegend()
print(p1)
ggsave(filename = "../plots/umaps/mural/mural_umap_major_celltypes.png", plot = p1, device = "png", width = 5, height = 5, dpi = 300)

p2 <- DimPlot(mural_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F, split.by = "tissue") + NoLegend()
print(p2)
ggsave(filename = "../plots/umaps/mural/mural_umap_major_celltypes_split.png", plot = p2, device = "png", width = 10, height = 5, dpi = 300)


p3 <- DimPlot(mural_named, reduction = "umap", label = F, repel = TRUE, raster = F) + NoLegend()
print(p3)
ggsave(filename = "../plots/umaps/mural/mural_umap_minor_celltypes_unlabeled.png", plot = p3, device = "png", width = 5, height = 5, dpi = 300)

p4 <- DimPlot(mural_named, reduction = "umap", label = F, repel = TRUE, raster = F, split.by = "tissue") + NoLegend()
print(p4)
ggsave(filename = "../plots/umaps/mural/mural_umap_minor_celltypes_split_unlabeled.png", plot = p4, device = "png", width = 10, height = 5, dpi = 300)

# https://pmc.ncbi.nlm.nih.gov/articles/PMC7137757/
# mostly used for VEC reference, but wanted to check a couple markers

### Interactive marker exploration (run manually in RStudio) ###
# FeaturePlot(mural_named, features = c("Acta2", "Tagln", "Myf9", "Kcnj8"), order = T)
# FeaturePlot(mural_named, features = c("Pdgfrb", "Des", "Mcam", "Baiap3"), order = T)

#### FIGURE S1 EXPORTS

p_supp_minor <- DimPlot(mural_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F, split.by = "tissue") +
  labs(x = "UMAP-1", y = "UMAP-2") +
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  NoLegend()


print(p_supp_minor)
ggsave(filename = "../plots/s1_umaps/mural_umap_minor_celltypes_split.png", plot = p_supp_minor, device = "png", width = 7, height = 4, dpi = 300)

ggsave(filename = "../plots/s1_umaps/mural_umap_minor_celltypes_split.pdf", plot = p_supp_minor, device = "pdf", width = 7, height = 4, dpi = 300)

message("06_main_mural_clustering complete.")
