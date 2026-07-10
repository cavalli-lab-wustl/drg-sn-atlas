#' ======================================================================
#' 06_main_fibro_clustering.R
#'
#' Subset, re-integrate, and cluster fibroblast populations at res=0.9.
#'
#' Inputs:  ../r_objects/combined_object_init_named.RDS
#' Outputs: ../r_objects/fibroblasts_object_init.RDS, ../metadata/fibroblasts_*_celltypes.rds
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

fibroblasts <- subset(merged_named, idents = c("fibro"))

options(future.globals.maxSize = 3e+09)
fibroblasts <- SCTransform(fibroblasts)
fibroblasts <- RunPCA(fibroblasts, npcs = 50)
fibroblasts <- IntegrateLayers(
  object = fibroblasts,
  method = RPCAIntegration,
  normalization.method = "SCT",
  k.weight = 80  # reduced from default 100; smallest batch has ~80 cells
)
fibroblasts <- FindNeighbors(fibroblasts, dims = 1:50, reduction = "integrated.dr")
fibroblasts <- FindClusters(fibroblasts, resolution = 0.5)
fibroblasts <- RunUMAP(fibroblasts, reduction = "integrated.dr", dims = 1:50)

# Plot UMAP, group by 'seurat_clusters'
dir.create("../plots/umaps/fibroblasts", recursive = TRUE, showWarnings = FALSE)
p1 <- DimPlot(fibroblasts, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, raster = F) + NoLegend()
print(p1)
ggsave("../plots/umaps/fibroblasts/fibroblasts_initial_umap.png", p1, width = 10, height = 10, units = "in", dpi = 300)

saveRDS(fibroblasts, "../r_objects/fibroblasts_object_init.RDS")

###########################################
###########################################
###                                     ###
###    CLUSTER RESOLUTION + CLUSTREE    ###
###                                     ###
###########################################
###########################################

fibroblasts <- readRDS("../r_objects/fibroblasts_object_init.RDS")

# Do not set RNA assay or normalize - FindClusters requires SCT assay

# Array of resolutions
dir.create("../plots/clustree/fibroblasts", recursive = TRUE, showWarnings = FALSE)
resolutions <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

for (res in resolutions) {
  # Find clusters and create plots
  fibroblasts <- FindClusters(fibroblasts, resolution = res)
  p <- DimPlot(fibroblasts, label = T, raster = F)

  # Save plots with and without legend
  ggplot2::ggsave(paste0("../plots/clustree/fibroblasts/clustree_fibroblasts_", res, "res.png"), plot = p, device = "png", width = 10, height = 10, dpi = 300)
}

# Run clustree to examine clustering layout
p <- clustree(fibroblasts, prefix = "SCT_snn_res.") + theme(text = element_text(size = 20))

# Saving the plot with specified parameters
ggsave(filename = "../plots/clustree/fibroblasts/clustree_fibroblasts_summary_plot.png", plot = p, device = "png", width = 18, height = 9, dpi = 300)


# Resolution 0.6 selected initially; later changed to 0.9 for finer subtypes

###########################################
###########################################
###                                     ###
###    METADATA EXPORT                  ###
###                                     ###
###########################################
###########################################

fibroblasts <- readRDS("../r_objects/fibroblasts_object_init.RDS")
fibroblasts <- FindClusters(fibroblasts, resolution = 0.9)

DefaultAssay(fibroblasts) <- "RNA"
fibroblasts <- NormalizeData(fibroblasts)
fibroblasts <- JoinLayers(fibroblasts)

# Create a directory for saving output files
output_directory <- "../markers/fibroblasts_init_full"
dir.create(output_directory, showWarnings = FALSE)
find_and_save_all_markers(fibroblasts, output_directory, min.pct = 0.1, logfc.threshold = 0.25, summary = T)


init_abundances <- table(fibroblasts@meta.data$study, Idents(fibroblasts))
write.csv(init_abundances, "../metadata/fibroblasts_initial_cellcounts.csv")


# Extract data
data <- FetchData(fibroblasts, vars = c("nFeature_RNA", "nCount_RNA"))
data$cluster <- Idents(fibroblasts)
average_values <- data %>%
  group_by(cluster) %>%
  summarise(
    Average_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
    Average_nCount_RNA = mean(nCount_RNA, na.rm = TRUE)
  )
# Export the result to a CSV file
write.csv(average_values, "../metadata/fibroblasts_avg_nFeaturenCount.csv", row.names = FALSE)

dir.create("../plots/umaps/fibroblasts_highlights", recursive = TRUE, showWarnings = FALSE)
highlightAndSaveAllClusters(fibroblasts, save_directory = "../plots/umaps/fibroblasts_highlights/")

# Extract the average gene expression for each cluster
avg_exp <- AverageExpression(fibroblasts)

# Restrict analysis to variable genes
fibroblasts <- FindVariableFeatures(object = fibroblasts)
variable_gene_names <- VariableFeatures(fibroblasts, assay = "RNA")

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
ggsave('../plots/correlations/fibroblasts_seuratclusters_correlation_matrix_hclust.pdf', plot=p, device='pdf', width=12, height=12, dpi=300)

### Interactive marker exploration (run manually in RStudio) ###
# FeaturePlot(fibroblasts, features = c("Col1a1", "Pcolce2", "Cdkn2a", "Cldn1"), order = T)
# FeaturePlot(fibroblasts, features = c("Lum", "Defb36", "Angptl7", "Phactr3"), order = T)
# FeaturePlot(fibroblasts, features = c("Sfrp4", "Fbln7", "Cfh", "Cited1"), order = T)
# FeaturePlot(fibroblasts, features = c("Cldn10", "Npy1r", "Mmp3", "Gpnmb"), order = T)
# FeaturePlot(fibroblasts, features = c("Kera", "Spink1", "Gjb2", "Cldn1"), order = T)
# FeaturePlot(fibroblasts, features = c("Chchd10", "Cox7a1", "Ocln", "Ppargc1a"), order = T)
# FeaturePlot(fibroblasts, features = c("Rorb", "Gpr1", "Edn1", "Smpd3"), order = T)
# FeaturePlot(fibroblasts, features = c("Sema3e", "Ackr2", "Dmrt2", "Stmn4"), order = T)
# FeaturePlot(fibroblasts, features = c("Xist", "Tsix", "Kdm5d", "Ddx3y"), order = T)
# FeaturePlot(fibroblasts, features = c("Pla1a", "Pcdh20", "Aif1l", "Fam167a"), order = T)
# FeaturePlot(fibroblasts, features = c("Engase", "Cxcl13", "Pde8a", "Npr1"), order = T)
# FeaturePlot(fibroblasts, features = c("Lypd6b", "Slc6a12", "Kcnj13", "Spink1"), order = T)
# FeaturePlot(fibroblasts, features = c("Fam78b", "Specc1", "Gna14", "Itgb4"), order = T)
# FeaturePlot(fibroblasts, features = c("Ripor2", "Ocln", "Col4a3", "Cldn1"), order = T)
# FeaturePlot(fibroblasts, features = c("Cgn", "Ocln", "Jam2", "Cldn1"), order = T) # BARRIER GENES!
# FeaturePlot(fibroblasts, features = c("Tagln", "Ccdc42", "Acta2", "Npnt"), order = T)
# FeaturePlot(fibroblasts, features = c("Msln", "Hcn4", "Folr1", "Wnt10a"), order = T)
# FeaturePlot(fibroblasts, features = c("Scube2", "Itih2", "Cldn22", "Tnnt3"), order = T)
# FeaturePlot(fibroblasts, features = c("Myh11", "Tpm2", "Mylk", "Mustn1"), order = T)
# FeaturePlot(fibroblasts, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "Meg3"), order = T)
# FeaturePlot(fibroblasts, features = c("Anxa8", "Gm17455", "Ociad2", "Tnc"), order = T)
# FeaturePlot(fibroblasts, features = c("Dio2", "Wnt2", "Dmkn", "Lyve1"), order = T)
# FeaturePlot(fibroblasts, features = c("P4ha3", "Cd109", "Nog", "Krtdap"), order = T)
# FeaturePlot(fibroblasts, features = c("Comp", "Cit", "Wnt9a", "Pcolce2"), order = T)
# FeaturePlot(fibroblasts, features = c("Cxcl1", "Btg2", "Ier2", "Junb"), order = T)
# FeaturePlot(fibroblasts, features = c("Jun", "Socs3", "Fos", "Cdkn2a"), order = T)
# FeaturePlot(fibroblasts, features = c("Ace2", "Mettl7b", "Krt19", "Sfrp5"), order = T)
# FeaturePlot(fibroblasts, features = c("Pabpc1l", "Cfap57", "Gm13261", "Lypd2"), order = T)
# FeaturePlot(fibroblasts, features = c("Perp", "Lad1", "Cldn22", "Hsd17b14"), order = T)
# FeaturePlot(fibroblasts, features = c("Efnb3", "Slc2a1", "Msln", "Cldn1"), order = T)
# FeaturePlot(fibroblasts, features = c("Col25a1", "Slc26a7", "Gm26740", "Slc47a2"), order = T)
# FeaturePlot(fibroblasts, features = c("Slc47a1", "Slc4a10", "Dlg2.1", "Slc35g1"), order = T)
# FeaturePlot(fibroblasts, features = c("Slit2", "Slc12a7", "Gm29478", "Col1a1"), order = T)
# FeaturePlot(fibroblasts, features = c("Glrb", "Hhip", "Cit", "Ptgdr"), order = T)
# FeaturePlot(fibroblasts, features = c("Birc5", "Top2a", "Mki67", "Aurkb"), order = T)
# FeaturePlot(fibroblasts, features = c("Ccl2", "Ccl7", "Cxcl2", "Cxcl10"), order = T)
# FeaturePlot(fibroblasts, features = c("Ptx3", "Tnfaip3", "Il6", "Ccl9"), order = T)
# FeaturePlot(fibroblasts, features = c("Cytl1", "Mgp", "Gata6", "Cdh13"), order = T)
# FeaturePlot(fibroblasts, features = c("Slc38a1", "Slc26a7", "Gdf10", "Ecrg4"), order = T)
# FeaturePlot(fibroblasts, features = c("Mgp", "Cst3", "Bgn", "Gdf10"), order = T)
# FeaturePlot(fibroblasts, features = c("Slc26a7", "Slc47a2", "Slc47a1", "Slc4a10"), order = T)
# FeaturePlot(fibroblasts, features = c("Pcolce2", "Col1a1", "Cldn1", "Cdkn2a"), order = T)

###########################################
###########################################
###                                     ###
###    CLUSTER NAMING                   ###
###                                     ###
###########################################
###########################################

fibroblasts_named <- fibroblasts

# Minor-level naming (30 clusters -> 18 subtypes + BATCH)
cluster_to_celltype <- c(
  "0" = "EM I",
  "1" = "PM I",
  "2" = "EF I",
  "3" = "EM II",
  "4" = "EM I",
  "5" = "EM III",
  "6" = "EF II",
  "7" = "EF III",
  "8" = "EF IV",
  "9" = "PM II",
  "10" = "BATCH",
  "11" = "EM IV",
  "12" = "PM II",
  "13" = "EM IFN",
  "14" = "EF II",
  "15" = "EM II",
  "16" = "EF V",
  "17" = "EF VI",
  "18" = "EM I",
  "19" = "EF VII",
  "20" = "F_CYC",
  "21" = "BATCH",
  "22" = "BATCH",
  "23" = "BATCH",
  "24" = "PM III",
  "25" = "EF IV",
  "26" = "EF VIII",
  "27" = "BATCH",
  "28" = "EF IFN",
  "29" = "BATCH"
)


fibroblasts_named <- RenameIdents(fibroblasts_named, cluster_to_celltype)

major_levels <- c(
  "EM I",
  "EM II",
  "EM III",
  "EM IV",
  "EM IFN",
  "EF I",
  "EF II",
  "EF III",
  "EF IV",
  "EF V",
  "EF VI",
  "EF VII",
  "EF VIII",
  "EF IFN",
  "PM I",
  "PM II",
  "PM III",
  "F_CYC",
  "BATCH"
)

Idents(fibroblasts_named) <- factor(Idents(fibroblasts_named), levels = major_levels)


# Export minor cell type idents (to be used in main object)
fibroblasts_labels <- Idents(fibroblasts_named)
fibroblasts_cells <- Cells(fibroblasts_named)

# Export minor cell type idents and cells
saveRDS(fibroblasts_labels, "../metadata/fibroblasts_minor_celltypes.rds")
saveRDS(fibroblasts_cells, "../metadata/fibroblasts_minor_cells.rds")

# REMOVE background, technical artifacts
fibroblasts_named <- subset(fibroblasts_named, idents = c("BATCH"), invert = TRUE)

# Quick UMAP rerun after doublet removal
fibroblasts_named <- RunUMAP(fibroblasts_named, reduction = "integrated.dr", dims = 1:50)

# Plot UMAP
p1 <- DimPlot(fibroblasts_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F) + NoLegend()
print(p1)
ggsave(filename = "../plots/umaps/fibroblasts/fibroblasts_umap_minor_celltypes.png", plot = p1, device = "png", width = 5, height = 5, dpi = 300)

p2 <- DimPlot(fibroblasts_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F, split.by = "tissue") + NoLegend()
print(p2)
ggsave(filename = "../plots/umaps/fibroblasts/fibroblasts_umap_minor_celltypes_split.png", plot = p2, device = "png", width = 10, height = 5, dpi = 300)

p3 <- DimPlot(fibroblasts_named, reduction = "umap", label = F, repel = TRUE, raster = F) + NoLegend()
print(p3)
ggsave(filename = "../plots/umaps/fibroblasts/fibroblasts_umap_minor_celltypes_unlabeled.png", plot = p3, device = "png", width = 5, height = 5, dpi = 300)

p4 <- DimPlot(fibroblasts_named, reduction = "umap", label = F, repel = TRUE, raster = F, split.by = "tissue") + NoLegend()
print(p4)
ggsave(filename = "../plots/umaps/fibroblasts/fibroblasts_umap_minor_celltypes_split_unlabeled.png", plot = p4, device = "png", width = 10, height = 5, dpi = 300)

#### FIGURE S1 EXPORTS

dir.create("../plots/s1_umaps", recursive = TRUE, showWarnings = FALSE)
p_supp_minor <- DimPlot(fibroblasts_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F, split.by = "tissue") +
  labs(x = "UMAP-1", y = "UMAP-2") +
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  NoLegend()


print(p_supp_minor)
ggsave(filename = "../plots/s1_umaps/fibroblasts_umap_minor_celltypes_split.png", plot = p_supp_minor, device = "png", width = 7, height = 4, dpi = 300)

ggsave(filename = "../plots/s1_umaps/fibroblasts_umap_minor_celltypes_split.pdf", plot = p_supp_minor, device = "pdf", width = 7, height = 4, dpi = 300)



fibroblasts_named <- fibroblasts

# Major-level naming (30 clusters -> 4 major types + BATCH)
# EM = endoneurial/mesenchymal, EF = epineurial fibroblast, PM = perineurial
cluster_to_celltype <- c(
  "0" = "EM",
  "1" = "PM",
  "2" = "EF",
  "3" = "EM",
  "4" = "EM",
  "5" = "EM",
  "6" = "EF",
  "7" = "EF",
  "8" = "EF",
  "9" = "PM",
  "10" = "BATCH",
  "11" = "EM",
  "12" = "PM",
  "13" = "EM",
  "14" = "EF",
  "15" = "EM",
  "16" = "EF",
  "17" = "EF",
  "18" = "EM",
  "19" = "EF",
  "20" = "F_CYC",
  "21" = "BATCH",
  "22" = "BATCH",
  "23" = "BATCH",
  "24" = "PM",
  "25" = "EF",
  "26" = "EF",
  "27" = "BATCH",
  "28" = "EF",
  "29" = "BATCH"
)


fibroblasts_named <- RenameIdents(fibroblasts_named, cluster_to_celltype)

major_levels <- c(
  "EM",
  "EF",
  "PM",
  "F_CYC",
  "BATCH"
)

Idents(fibroblasts_named) <- factor(Idents(fibroblasts_named), levels = major_levels)


# Export idents (to be used in main object)
fibroblasts_labels <- Idents(fibroblasts_named)
fibroblasts_cells <- Cells(fibroblasts_named)

# Save major cell type idents to load later for transfer to main object
saveRDS(fibroblasts_labels, "../metadata/fibroblasts_major_celltypes.rds")
saveRDS(fibroblasts_cells, "../metadata/fibroblasts_major_cells.rds")

# REMOVE background, technical artifacts
fibroblasts_named <- subset(fibroblasts_named, idents = c("BATCH"), invert = TRUE)

# Quick UMAP rerun after doublet removal
fibroblasts_named <- RunUMAP(fibroblasts_named, reduction = "integrated.dr", dims = 1:50)

# Plot UMAP
p1 <- DimPlot(fibroblasts_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F) + NoLegend()
print(p1)
ggsave(filename = "../plots/umaps/fibroblasts/fibroblasts_umap_major_celltypes.png", plot = p1, device = "png", width = 5, height = 5, dpi = 300)

p2 <- DimPlot(fibroblasts_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F, split.by = "tissue") + NoLegend()
print(p2)
ggsave(filename = "../plots/umaps/fibroblasts/fibroblasts_umap_major_celltypes_split.png", plot = p2, device = "png", width = 10, height = 5, dpi = 300)



#### FIGURE S1 EXPORTS

p_supp_major <- DimPlot(fibroblasts_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F, split.by = "tissue") +
  labs(x = "UMAP-1", y = "UMAP-2") +
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  NoLegend()


print(p_supp_major)
ggsave(filename = "../plots/s1_umaps/fibroblasts_umap_major_celltypes_split.png", plot = p_supp_major, device = "png", width = 7, height = 4, dpi = 300)

ggsave(filename = "../plots/s1_umaps/fibroblasts_umap_major_celltypes_split.pdf", plot = p_supp_major, device = "pdf", width = 7, height = 4, dpi = 300)

message("06_main_fibro_clustering complete.")
