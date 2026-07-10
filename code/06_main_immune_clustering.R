#' ======================================================================
#' 06_main_immune_clustering.R
#'
#' Subset, re-integrate, and cluster immune populations at res=0.9 (k.weight=39).
#'
#' Inputs:  ../r_objects/combined_object_init_named.RDS
#' Outputs: ../r_objects/immune_object_init.RDS, ../metadata/immune_*_celltypes.rds
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

# Ear2 https://pubmed.ncbi.nlm.nih.gov/35484716/
immune <- subset(merged_named, idents = c("immune"))

options(future.globals.maxSize = 3e+09)
immune <- SCTransform(immune)
immune <- RunPCA(immune, npcs = 50)
immune <- IntegrateLayers(
  object = immune,
  method = RPCAIntegration,
  normalization.method = "SCT",
  k.weight = 39  # reduced from default 100; smallest immune batch has ~39 cells
)
immune <- FindNeighbors(immune, dims = 1:50, reduction = "integrated.dr")
immune <- FindClusters(immune, resolution = 0.5)
immune <- RunUMAP(immune, reduction = "integrated.dr", dims = 1:50)


# Plot UMAP, group by 'seurat_clusters'
p1 <- DimPlot(immune, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, raster = F) + NoLegend()

# Display plot
print(p1)

dir.create("../plots/umaps/immune", recursive = TRUE, showWarnings = FALSE)
ggsave("../plots/umaps/immune/immune_initial_umap.png", p1, width = 10, height = 10, units = "in", dpi = 300)

saveRDS(immune, "../r_objects/immune_object_init.RDS")


###########################################
###########################################
###                                     ###
###    CLUSTER RESOLUTION + CLUSTREE    ###
###                                     ###
###########################################
###########################################

immune <- readRDS("../r_objects/immune_object_init.RDS")

# Do not set RNA assay or normalize - FindClusters requires SCT assay

# Array of resolutions
resolutions <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

dir.create("../plots/clustree/immune", recursive = TRUE, showWarnings = FALSE)
for (res in resolutions) {
  # Find clusters and create plots
  immune <- FindClusters(immune, resolution = res)
  p <- DimPlot(immune, label = T, raster = F)

  # Save plots with and without legend
  ggplot2::ggsave(paste0("../plots/clustree/immune/clustree_immune_", res, "res.png"), plot = p, device = "png", width = 10, height = 10, dpi = 300)
}

# Run clustree to examine clustering layout
p <- clustree(immune, prefix = "SCT_snn_res.") + theme(text = element_text(size = 20))

# Saving the plot with specified parameters
ggsave(filename = "../plots/clustree/immune/clustree_immune_summary_plot.png", plot = p, device = "png", width = 18, height = 9, dpi = 300)


# Resolution 0.9 selected

immune <- readRDS("../r_objects/immune_object_init.RDS")
immune <- FindClusters(immune, resolution = 0.9)

DefaultAssay(immune) <- "RNA"
immune <- NormalizeData(immune)
immune <- JoinLayers(immune)

# Create a directory for saving output files
output_directory <- "../markers/immune_init_full"
dir.create(output_directory, showWarnings = FALSE)
find_and_save_all_markers(immune, output_directory, min.pct = 0.1, logfc.threshold = 0.25, summary = T)


init_abundances <- table(immune@meta.data$study, Idents(immune))
write.csv(init_abundances, "../metadata/immune_initial_cellcounts.csv")


# Extract data
data <- FetchData(immune, vars = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
data$cluster <- Idents(immune)
average_values <- data %>%
  group_by(cluster) %>%
  summarise(
    Average_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
    Average_nCount_RNA = mean(nCount_RNA, na.rm = TRUE),
    Average_percentmt = mean(percent.mt, na.rm = TRUE)
  )
# Export the result to a CSV file
write.csv(average_values, "../metadata/immune_avg_nFeaturenCountpercentmt.csv", row.names = FALSE)

# Extract the average gene expression for each cluster
avg_exp <- AverageExpression(immune)

# Restrict analysis to variable genes
immune <- FindVariableFeatures(object = immune)
variable_gene_names <- VariableFeatures(immune, assay = "RNA")

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

dir.create("../plots/correlations", recursive = TRUE, showWarnings = FALSE)
ggsave('../plots/correlations/immune_seuratclusters_correlation_matrix_hclust.pdf', plot=p, device='pdf', width=12, height=12, dpi=300)

dir.create("../plots/umaps/immune_highlights", recursive = TRUE, showWarnings = FALSE)
highlightAndSaveAllClusters(immune, save_directory = "../plots/umaps/immune_highlights/")

### Interactive marker exploration (run manually in RStudio) ###
# FeaturePlot(immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "Meg3"), order = T)
# FeaturePlot(immune, features = c("Birc5", "Top2a", "Mki67", "Aurkb"), order = T)
# FeaturePlot(immune, features = c("Cd68", "Cd63", "Fabp5", "Ptprc"), order = T)
# FeaturePlot(immune, features = c("Cxcr1", "Cxcr2", "Atf3", "Fos"), order = T)
# FeaturePlot(immune, features = c("Il10", "Tnf", "Atf3", "Fos"), order = T)
# FeaturePlot(immune, features = c("Fcna", "Cd209f", "Cd209g", "Cd209b"), order = T)
# FeaturePlot(immune, features = c("Ccl7", "Ccl2", "Ccl12", "Ch25h"), order = T)
# FeaturePlot(immune, features = c("Akr1c18", "Lyz1", "Ear2", "Aldh1a2"), order = T)
# FeaturePlot(immune, features = c("Cd8a", "Cd68", "Ms4a1", "Cd79b"), order = T)
# FeaturePlot(immune, features = c("Cd8a", "Trac", "Ms4a1", "Cd4"), order = T)
# FeaturePlot(immune, features = c("Flt1", "Prom1", "Ms4a1", "Cd4"), order = T)
# FeaturePlot(immune, features = c("Gpnmb", "Fabp5", "Spp1", "Arg1"), order = T) # MAC FIB
# FeaturePlot(immune, features = c("Cd209f", "Cd209g", "Cd209d", "Arg1"), order = T) # epiMAC (Giger)
# FeaturePlot(immune, features = c("Cx3cr1", "Trem2", "Lilra5", "Unc93b1"), order = T) # endoMAC (Giger)
# FeaturePlot(immune, features = c("Chil3"), order = T) # Monocytes
# FeaturePlot(immune, features = c("Cxcr2"), order = T) # Granulocytes
# FeaturePlot(immune, features = c("Cd163"), order = T) # PNS tissue-resident Mac (Giger)
# FeaturePlot(immune, features = c("Fcna", "Ccl8", "Cd209f", "Cd209g"), order = T) # epiMAC (Giger)
# FeaturePlot(immune, features = c("Fcer1a", "Cpa3", "Kit", "Unc93b1"), order = T) # epiMAC (Giger)
# FeaturePlot(immune, features = c("Cd74", "H2-Ab1", "H2-Aa", "Ifit3"), order = T)
# FeaturePlot(immune, features = c("Pla2g2d", "Ccl8", "Gpx3", "Tslp"), order = T) # Cluster 3
# FeaturePlot(immune, features = c("Cd209a", "H2-Aa", "Cd68", "H2-Eb1"), order = T) # moDC
# FeaturePlot(immune, features = c("Gal3st4", "Capn3", "Ermap", "F11r"), order = T) # Cluster 20
# FeaturePlot(immune, features = c("Il2rb", "Klrc1", "Trbc1", "Cd79a"), order = T) # T/NK/B
#
# FeaturePlot(immune, features = c("Ncam1", "Cd7", "Eomes", "Klra8"), order = T)
# FeaturePlot(immune, features = c("Ncr1", "Klrb1c", "Ccr7", "Gzmb"), order = T)
# FeaturePlot(immune, features = c("Cox6a2", "Cd7", "Ly6d", "Klra17"), order = T)
# FeaturePlot(immune, features = c("Fn1", "Vcan", "Tppp3", "Hopx"), order = T)
#
# FeaturePlot(immune_named, features = c("Birc5", "Top2a", "Mki67", "Aurkb"), order = T)

# 20250721 - iMoonglia subset update
imoonglia <- subset(immune, idents = c("11"))

DefaultAssay(imoonglia) <- "SCT"

imoonglia <- FindNeighbors(imoonglia, dims = 1:50)
imoonglia <- FindClusters(imoonglia, resolution = 0.2, reduction = "integrated.dr")  # low resolution — only 3 subclusters in small IMG population
imoonglia <- RunUMAP(imoonglia, reduction = "integrated.dr", dims = 1:50)

DefaultAssay(imoonglia) <- "RNA"

# Plot UMAP, group by 'seurat_clusters'
p1 <- DimPlot(imoonglia, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, raster = F)

# Display plot
print(p1)


# renaming all clusters
cluster_to_celltype <- c(
  "0" = "IMG",
  "1" = "BATCH",
  "2" = "BATCH"
)

imoonglia <- RenameIdents(imoonglia, cluster_to_celltype)

# Export idents (to be used in main object)
imoonglia_labels <- Idents(imoonglia)
imoonglia_cells <- Cells(imoonglia)

# Export minor cell type idents and cells
# Since imoonglia are tracked at both levels in the same manner these files are the same
# but used in different label transfer steps for naming consistency
saveRDS(imoonglia_labels, "../metadata/imoonglia_minor_celltypes.rds")
saveRDS(imoonglia_cells, "../metadata/imoonglia_minor_cells.rds")

saveRDS(imoonglia_labels, "../metadata/imoonglia_major_celltypes.rds")
saveRDS(imoonglia_cells, "../metadata/imoonglia_major_cells.rds")

### Interactive marker exploration: iMoonglia and cross-object (run manually in RStudio) ###
# FeaturePlot(imoonglia, features = c("Fabp7", "Plp1", "Ncmap", "Aif1")) &
#   theme(
#     plot.title  = element_text(size = 20, face = "italic", hjust = 0.5),
#     axis.text   = element_blank(),
#     axis.ticks  = element_blank(),
#     axis.title.x = element_text(size = 10, hjust = 0),
#     axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
#   )
#
# FeaturePlot(immune, features = c("Fabp7", "Plp1", "Ncmap", "Aif1")) &
#   theme(
#     plot.title  = element_text(size = 20, face = "italic", hjust = 0.5),
#     axis.text   = element_blank(),
#     axis.ticks  = element_blank(),
#     axis.title.x = element_text(size = 10, hjust = 0),
#     axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
#   )
#
# FeaturePlot(immune, features = c("Syt1", "Snap25", "Calca", "Prph")) &
#   theme(
#     plot.title  = element_text(size = 20, face = "italic", hjust = 0.5),
#     axis.text   = element_blank(),
#     axis.ticks  = element_blank(),
#     axis.title.x = element_text(size = 10, hjust = 0),
#     axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
#   )
#
# FeaturePlot(immune, features = c("Col1a1", "Lum", "Col6a1", "Pdgfra")) &
#   theme(
#     plot.title  = element_text(size = 20, face = "italic", hjust = 0.5),
#     axis.text   = element_blank(),
#     axis.ticks  = element_blank(),
#     axis.title.x = element_text(size = 10, hjust = 0),
#     axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
#   )
#
# FeaturePlot(immune, features = c("Plp1", "Fabp7", "Bcan", "Kcnj10")) &
#   theme(
#     plot.title  = element_text(size = 20, face = "italic", hjust = 0.5),
#     axis.text   = element_blank(),
#     axis.ticks  = element_blank(),
#     axis.title.x = element_text(size = 10, hjust = 0),
#     axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
#   )
#
# FeaturePlot(immune, features = c("Plp1", "Ncmap", "Gldn", "Pmp22")) &
#   theme(
#     plot.title  = element_text(size = 20, face = "italic", hjust = 0.5),
#     axis.text   = element_blank(),
#     axis.ticks  = element_blank(),
#     axis.title.x = element_text(size = 10, hjust = 0),
#     axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
#   )
#
# FeaturePlot(immune, features = c("Plp1", "Scn7a", "Ntrk3", "Pou3f1")) &
#   theme(
#     plot.title  = element_text(size = 20, face = "italic", hjust = 0.5),
#     axis.text   = element_blank(),
#     axis.ticks  = element_blank(),
#     axis.title.x = element_text(size = 10, hjust = 0),
#     axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
#   )
#
# FeaturePlot(immune, features = c("Tmem132c", "Anks1b", "Sorbs1", "Psd3")) &
#   theme(
#     plot.title  = element_text(size = 20, face = "italic", hjust = 0.5),
#     axis.text   = element_blank(),
#     axis.ticks  = element_blank(),
#     axis.title.x = element_text(size = 10, hjust = 0),
#     axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
#   )
#
# FeaturePlot(merged_named, features = c("Tmem132c", "Anks1b", "Sorbs1", "Psd3")) &
#   theme(
#     plot.title  = element_text(size = 20, face = "italic", hjust = 0.5),
#     axis.text   = element_blank(),
#     axis.ticks  = element_blank(),
#     axis.title.x = element_text(size = 10, hjust = 0),
#     axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
#   )
#
# FeaturePlot(merged_named, features = c("Plekhg5", "Taco1", "Fgfr2", "Fhit")) &
#   theme(
#     plot.title  = element_text(size = 20, face = "italic", hjust = 0.5),
#     axis.text   = element_blank(),
#     axis.ticks  = element_blank(),
#     axis.title.x = element_text(size = 10, hjust = 0),
#     axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
#   )
#
# FeaturePlot(merged_named, features = c("Hdac9", "Foxo3", "Mctp1", "P3h2")) &
#   theme(
#     plot.title  = element_text(size = 20, face = "italic", hjust = 0.5),
#     axis.text   = element_blank(),
#     axis.ticks  = element_blank(),
#     axis.title.x = element_text(size = 10, hjust = 0),
#     axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
#   )
#
# FeaturePlot(immune, features = c("Ngfr", "Nes", "Shh", "Cxcl10")) &
#   theme(
#     plot.title  = element_text(size = 20, face = "italic", hjust = 0.5),
#     axis.text   = element_blank(),
#     axis.ticks  = element_blank(),
#     axis.title.x = element_text(size = 10, hjust = 0),
#     axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
#   )

immune_named <- immune

# Major-level naming
cluster_to_celltype <- c(
  "0" = "MAC",
  "1" = "MAC",
  "2" = "DC",
  "3" = "MAC",
  "4" = "MAC",
  "5" = "MAC",
  "6" = "LEUK",
  "7" = "LEUK",
  "8" = "BATCH",
  "9" = "EOS",
  "10" = "MAC",
  "11" = "BATCH",
  "12" = "MAC",
  "13" = "MONO",
  "14" = "H_CYC",
  "15" = "MAC",
  "16" = "H_CYC",
  "17" = "NT",
  "18" = "DC",
  "19" = "MAST",
  "20" = "MAC",
  "21" = "BATCH",
  "22" = "NT",
  "23" = "MAC",
  "24" = "MAC",
  "25" = "DC",
  "26" = "DC",
  "27" = "BATCH"
)

immune_named <- RenameIdents(immune_named, cluster_to_celltype)

# ADDED 20251022 FOR IMG
Idents(immune_named, cells = imoonglia_cells) <- imoonglia_labels


major_levels <- c(
  "MAC",
  "MONO",
  "IMG",
  "DC",
  "NT",
  "EOS",
  "MAST",
  "LEUK",
  "H_CYC",
  "BATCH"
)

Idents(immune_named) <- factor(Idents(immune_named), levels = major_levels)

# Export idents (to be used in main object)
immune_labels <- Idents(immune_named)
immune_cells <- Cells(immune_named)

# Export major cell type idents and cells
saveRDS(immune_labels, "../metadata/immune_major_celltypes.rds")
saveRDS(immune_cells, "../metadata/immune_major_cells.rds")

# REMOVE background, technical artifacts
immune_named <- subset(immune_named, idents = c("BATCH"), invert = TRUE)

# Quick UMAP rerun after doublet removal
immune_named <- RunUMAP(immune_named, reduction = "integrated.dr", dims = 1:50)

print(DimPlot(immune_named, label = T))

#### FIGURE S1 EXPORTS

p_supp_major <- DimPlot(immune_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F, split.by = "tissue") +
  labs(x = "UMAP-1", y = "UMAP-2") +
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  NoLegend()


print(p_supp_major)
dir.create("../plots/s1_umaps", recursive = TRUE, showWarnings = FALSE)
ggsave(filename = "../plots/s1_umaps/immune_umap_major_celltypes_split.png", plot = p_supp_major, device = "png", width = 7, height = 4, dpi = 300)

ggsave(filename = "../plots/s1_umaps/immune_umap_major_celltypes_split.pdf", plot = p_supp_major, device = "pdf", width = 7, height = 4, dpi = 300)


immune_named <- immune

# Minor-level naming
cluster_to_celltype <- c(
  "0" = "MAC M1",
  "1" = "MAC M2",
  "2" = "MoDC",
  "3" = "MAC A",
  "4" = "MAC FIB",
  "5" = "MAC B",
  "6" = "T NK",
  "7" = "BC",
  "8" = "BATCH",
  "9" = "EOS",
  "10" = "MAC A",
  "11" = "BATCH",
  "12" = "MAC A",
  "13" = "MONO",
  "14" = "H_CYC",
  "15" = "epiMAC",
  "16" = "H_CYC",
  "17" = "NT A",
  "18" = "cDC",
  "19" = "MAST",
  "20" = "MLC", # RENAMED 20251107
  "21" = "BATCH",
  "22" = "NT B",
  "23" = "MAC M1",
  "24" = "MAC M1",
  "25" = "pDC",
  "26" = "DCx",
  "27" = "BATCH"
)

immune_named <- RenameIdents(immune_named, cluster_to_celltype)

# ADDED 20251110 FOR IMG
Idents(immune_named, cells = imoonglia_cells) <- imoonglia_labels


major_levels <- c(
  "MAC M1",
  "MAC M2",
  "MAC A",
  "MAC B",
  "MLC", # RENAMED 20251107
  "epiMAC",
  "MAC FIB",
  "DCx",
  "cDC",
  "pDC",
  "MoDC",
  "MAST",
  "MONO",
  "T NK",
  "BC",
  "H_CYC",
  "EOS",
  "NT A",
  "NT B",
  "IMG",
  "BATCH"
)

Idents(immune_named) <- factor(Idents(immune_named), levels = major_levels)

# Export idents (to be used in main object)
immune_labels <- Idents(immune_named)
immune_cells <- Cells(immune_named)

# Export minor cell type idents and cells
saveRDS(immune_labels, "../metadata/immune_minor_celltypes.rds")
saveRDS(immune_cells, "../metadata/immune_minor_cells.rds")

# REMOVE background, technical artifacts
immune_named <- subset(immune_named, idents = c("BATCH"), invert = TRUE)

# Quick UMAP rerun after doublet removal
immune_named <- RunUMAP(immune_named, reduction = "integrated.dr", dims = 1:50)

# Plot UMAP
p1 <- DimPlot(immune_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F) + NoLegend()
print(p1)
ggsave("../plots/umaps/immune/immune_umap_minor_celltypes.png", p1, width = 5, height = 5, units = "in", dpi = 300)

p2 <- DimPlot(immune_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F, split.by = "tissue") + NoLegend()
print(p2)
ggsave(filename = "../plots/umaps/immune/immune_umap_minor_celltypes_split.png", plot = p2, device = "png", width = 10, height = 5, dpi = 300)


p3 <- DimPlot(immune_named, reduction = "umap", label = F, repel = TRUE, raster = F) + NoLegend()
print(p3)
ggsave(filename = "../plots/umaps/immune/immune_umap_minor_celltypes_unlabeled.png", plot = p3, device = "png", width = 5, height = 5, dpi = 300)

p4 <- DimPlot(immune_named, reduction = "umap", label = F, repel = TRUE, raster = F, split.by = "tissue") + NoLegend()
print(p4)
ggsave(filename = "../plots/umaps/immune/immune_umap_minor_celltypes_split_unlabeled.png", plot = p4, device = "png", width = 10, height = 5, dpi = 300)

#### FIGURE S1 EXPORTS

p_supp_minor <- DimPlot(immune_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F, split.by = "tissue") +
  labs(x = "UMAP-1", y = "UMAP-2") +
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  NoLegend()


print(p_supp_minor)
ggsave(filename = "../plots/s1_umaps/immune_umap_minor_celltypes_split.png", plot = p_supp_minor, device = "png", width = 7, height = 4, dpi = 300)

ggsave(filename = "../plots/s1_umaps/immune_umap_minor_celltypes_split.pdf", plot = p_supp_minor, device = "pdf", width = 7, height = 4, dpi = 300)

###########################################
###########################################
###                                     ###
###   ABUNDANCE PLOTS                   ###
###                                     ###
###########################################
###########################################

immune_named$celltype <- Idents(immune_named)

# Extract metadata from the Seurat object
metadata <- immune_named@meta.data

# 1. Calculate counts of cells by treatment and celltype
cell_counts <- metadata %>%
  group_by(treatment_simple, celltype) %>%
  summarize(count = n(), .groups = "drop")

# 2. Calculate relative abundance for each treatment
cell_counts <- cell_counts %>%
  group_by(treatment_simple) %>%
  mutate(total = sum(count),
         rel_abundance = count / total) %>%
  ungroup()

# 3a. Plot: Stacked Bar Plot of Relative Abundance
p_stacked <- ggplot(cell_counts, aes(x = treatment_simple, y = rel_abundance, fill = celltype)) +
  geom_bar(stat = "identity") +
  labs(x = "Treatment", y = "Relative Abundance",
       title = "Stacked Bar Plot of Relative Abundance by Cell Type and Treatment") +
  theme_minimal()

# 3b. Plot: Grouped (Dodge) Bar Plot of Relative Abundance
p_grouped <- ggplot(cell_counts, aes(x = treatment_simple, y = rel_abundance, fill = celltype)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Treatment", y = "Relative Abundance",
       title = "Grouped Bar Plot of Relative Abundance by Cell Type and Treatment") +
  theme_minimal()

# 4. Calculate relative change in abundance for each non-WT condition compared to WT.
#    First, extract the WT relative abundance for each cell type.
wt_freq <- cell_counts %>%
  filter(treatment_simple == "WT") %>%
  dplyr::select(celltype, WT_rel_abundance = rel_abundance)

# Join the WT frequencies with non-WT conditions and compute the relative change
cell_change <- cell_counts %>%
  filter(treatment_simple != "WT") %>%
  left_join(wt_freq, by = "celltype") %>%
  mutate(rel_change = (rel_abundance - WT_rel_abundance) / WT_rel_abundance * 100)

# 5. Plot: Relative Change in Abundance Compared to WT (grouped bar plot)
p_change <- ggplot(cell_change, aes(x = celltype, y = rel_change, fill = treatment_simple)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Cell Type", y = "Relative Change (%)",
       title = "Relative Change in Abundance (Non-WT vs. WT)") +
  theme_minimal()

# Display the plots
print(p_stacked)
print(p_grouped)
print(p_change)


### COMMENTED OUT: SAM signature analysis (run manually if needed) ###
#
# ###########################################
# ###########################################
# ###                                     ###
# ###   SAM SIGNATURE PLOTS               ###
# ###                                     ###
# ###########################################
# ###########################################
#
# ## 0) COMPUTE GLOBAL DENSITY PARAMETERS  --------------------------------
# sig <- c("Fabp5", "Cd63", "Trem2", "Cd9", "Spp1", "Gpnmb")
#
# # Alternative approach: Calculate gene set score first
# immune_named <- AddModuleScore(immune_named, features = list(sig), name = "GeneSet")
#
# # Get UMAP coordinates and gene set scores
# umap_coords <- Embeddings(immune_named, reduction = "umap")
# gene_scores <- immune_named$GeneSet1
#
# # Create a data frame with coordinates and scores
# full_plot_data <- data.frame(
#   UMAP_1 = umap_coords[, 1],
#   UMAP_2 = umap_coords[, 2],
#   gene_score = gene_scores,
#   treatment_simple = immune_named$treatment_simple
# )
#
# # Calculate density using gene scores
# max_score <- max(full_plot_data$gene_score, na.rm = TRUE)
# min_score <- min(full_plot_data$gene_score, na.rm = TRUE)
#
# ## 1) CHOOSE A PALETTE ---------------------------------------------------
# palette_choice <- "viridis"
# scale_color <- switch(
#   palette_choice,
#   viridis = scale_color_viridis_c(option = "viridis",
#                                   limits = c(min_score, max_score),
#                                   oob = squish),
#   magma   = scale_color_viridis_c(option = "magma",
#                                   limits = c(min_score, max_score),
#                                   oob = squish),
#   inferno = scale_color_viridis_c(option = "inferno",
#                                   limits = c(min_score, max_score),
#                                   oob = squish),
#   RdBu    = scale_color_gradientn(
#     colours = rev(brewer.pal(11, "RdBu")),
#     limits  = c(min_score, max_score),
#     oob     = squish),
#   blues   = scale_color_gradient(
#     low = "grey90", high = "navy",
#     limits = c(min_score, max_score),
#     oob = squish)
# )
#
# ## ---------- SHARED SETTINGS --------------------------------------------------
# # One colour scale for every panel -- adjust to taste
# score_range  <- range(full_plot_data$gene_score, na.rm = TRUE)
# scale_colour <- scale_colour_viridis_c(limits = score_range, option = "D")
#
# # Treatments to plot
# treatments <- c("WT", "AIH", "DRC", "SCI", "SNC")
#
# ## ---------- BUILD PLOTS ------------------------------------------------------
# plot_list <- lapply(treatments, function(trt) {
#
#   trt_dat <- full_plot_data[full_plot_data$treatment_simple == trt, ]
#   trt_dat <- trt_dat[order(trt_dat$gene_score), ]        # light-to-dark over-plotting
#
#   ggplot(trt_dat, aes(UMAP_1, UMAP_2, colour = gene_score)) +
#     geom_point(size = 0.7) +
#     scale_colour +
#     labs(colour = "Gene-set\nscore") +
#     ggtitle(trt) +
#     theme_void() +
#     theme(
#       plot.title   = element_text(hjust = 0.5, size = 14),
#       legend.title = element_text(size = 10),
#       legend.text  = element_text(size = 8),
#       legend.position = "right"
#     )
# })
#
# names(plot_list) <- treatments   # convenient for retrieval: plot_list[["WT"]]
#
# ## ---------- OPTIONAL: SHOW ALL TOGETHER -------------------------------------
# # Wrap into a single figure (3x2 grid; empty 6th slot suppressed)
# combined_plt <- wrap_plots(plot_list, nrow = 1, guides = "collect") &
#   theme(legend.position = "right")
# combined_plt         # print or save with ggsave()
#
# # Save the plot
# ggplot2::ggsave("../plots/DRGSNMI_SAM_signature_split.pdf", plot = combined_plt, device = "pdf", width = 15.1, height = 3, dpi = 300)
#
# ##### DRG ONLY
#
# drg_immune_named <- subset(immune_named, subset = tissue == "DRG")
#
# ## 0) COMPUTE GLOBAL DENSITY PARAMETERS  --------------------------------
# sig <- c("Fabp5", "Cd63", "Trem2", "Cd9", "Spp1", "Gpnmb")
#
# # Alternative approach: Calculate gene set score first
# drg_immune_named <- AddModuleScore(drg_immune_named, features = list(sig), name = "GeneSet")
#
# # Get UMAP coordinates and gene set scores
# umap_coords <- Embeddings(drg_immune_named, reduction = "umap")
# gene_scores <- drg_immune_named$GeneSet1
#
# # Create a data frame with coordinates and scores
# full_plot_data <- data.frame(
#   UMAP_1 = umap_coords[, 1],
#   UMAP_2 = umap_coords[, 2],
#   gene_score = gene_scores,
#   treatment_simple = drg_immune_named$treatment_simple
# )
#
# # Calculate density using gene scores
# max_score <- max(full_plot_data$gene_score, na.rm = TRUE)
# min_score <- min(full_plot_data$gene_score, na.rm = TRUE)
#
# ## 1) CHOOSE A PALETTE ---------------------------------------------------
# palette_choice <- "viridis"
# scale_color <- switch(
#   palette_choice,
#   viridis = scale_color_viridis_c(option = "viridis",
#                                   limits = c(min_score, max_score),
#                                   oob = squish),
#   magma   = scale_color_viridis_c(option = "magma",
#                                   limits = c(min_score, max_score),
#                                   oob = squish),
#   inferno = scale_color_viridis_c(option = "inferno",
#                                   limits = c(min_score, max_score),
#                                   oob = squish),
#   RdBu    = scale_color_gradientn(
#     colours = rev(brewer.pal(11, "RdBu")),
#     limits  = c(min_score, max_score),
#     oob     = squish),
#   blues   = scale_color_gradient(
#     low = "grey90", high = "navy",
#     limits = c(min_score, max_score),
#     oob = squish)
# )
#
# ## ---------- SHARED SETTINGS --------------------------------------------------
# # One colour scale for every panel -- adjust to taste
# score_range  <- range(full_plot_data$gene_score, na.rm = TRUE)
# scale_colour <- scale_colour_viridis_c(limits = score_range, option = "D")
#
# # Treatments to plot
# treatments <- c("WT", "AIH", "DRC", "SCI", "SNC")
#
# ## ---------- BUILD PLOTS ------------------------------------------------------
# plot_list <- lapply(treatments, function(trt) {
#
#   trt_dat <- full_plot_data[full_plot_data$treatment_simple == trt, ]
#   trt_dat <- trt_dat[order(trt_dat$gene_score), ]        # light-to-dark over-plotting
#
#   ggplot(trt_dat, aes(UMAP_1, UMAP_2, colour = gene_score)) +
#     geom_point(size = 0.7) +
#     scale_colour +
#     labs(colour = "Gene-set\nscore") +
#     ggtitle(trt) +
#     theme_void() +
#     theme(
#       plot.title   = element_text(hjust = 0.5, size = 14),
#       legend.title = element_text(size = 10),
#       legend.text  = element_text(size = 8),
#       legend.position = "right"
#     )
# })
#
# names(plot_list) <- treatments   # convenient for retrieval: plot_list[["WT"]]
#
# ## ---------- OPTIONAL: SHOW ALL TOGETHER -------------------------------------
# # Wrap into a single figure (3x2 grid; empty 6th slot suppressed)
# combined_plt <- wrap_plots(plot_list, nrow = 1, guides = "collect") &
#   theme(legend.position = "right")
# combined_plt         # print or save with ggsave()
#
# # Save the plot
# ggplot2::ggsave("../plots/DRGSNMI_SAM_signature_split_DRG.pdf", plot = combined_plt, device = "pdf", width = 15.1, height = 3, dpi = 300)
#
# ##### SN ONLY
#
# sn_immune_named <- subset(immune_named, subset = tissue == "SN")
#
# ## 0) COMPUTE GLOBAL DENSITY PARAMETERS  --------------------------------
# sig <- c("Fabp5", "Cd63", "Trem2", "Cd9", "Spp1", "Gpnmb")
#
# # Alternative approach: Calculate gene set score first
# sn_immune_named <- AddModuleScore(sn_immune_named, features = list(sig), name = "GeneSet")
#
# # Get UMAP coordinates and gene set scores
# umap_coords <- Embeddings(sn_immune_named, reduction = "umap")
# gene_scores <- sn_immune_named$GeneSet1
#
# # Create a data frame with coordinates and scores
# full_plot_data <- data.frame(
#   UMAP_1 = umap_coords[, 1],
#   UMAP_2 = umap_coords[, 2],
#   gene_score = gene_scores,
#   treatment_simple = sn_immune_named$treatment_simple
# )
#
# # Calculate density using gene scores
# max_score <- max(full_plot_data$gene_score, na.rm = TRUE)
# min_score <- min(full_plot_data$gene_score, na.rm = TRUE)
#
# ## 1) CHOOSE A PALETTE ---------------------------------------------------
# palette_choice <- "viridis"
# scale_color <- switch(
#   palette_choice,
#   viridis = scale_color_viridis_c(option = "viridis",
#                                   limits = c(min_score, max_score),
#                                   oob = squish),
#   magma   = scale_color_viridis_c(option = "magma",
#                                   limits = c(min_score, max_score),
#                                   oob = squish),
#   inferno = scale_color_viridis_c(option = "inferno",
#                                   limits = c(min_score, max_score),
#                                   oob = squish),
#   RdBu    = scale_color_gradientn(
#     colours = rev(brewer.pal(11, "RdBu")),
#     limits  = c(min_score, max_score),
#     oob     = squish),
#   blues   = scale_color_gradient(
#     low = "grey90", high = "navy",
#     limits = c(min_score, max_score),
#     oob = squish)
# )
#
# ## ---------- SHARED SETTINGS --------------------------------------------------
# # One colour scale for every panel -- adjust to taste
# score_range  <- range(full_plot_data$gene_score, na.rm = TRUE)
# scale_colour <- scale_colour_viridis_c(limits = score_range, option = "D")
#
# # Treatments to plot
# treatments <- c("WT", "AIH", "DRC", "SCI", "SNC")
#
# ## ---------- BUILD PLOTS ------------------------------------------------------
# plot_list <- lapply(treatments, function(trt) {
#
#   trt_dat <- full_plot_data[full_plot_data$treatment_simple == trt, ]
#   trt_dat <- trt_dat[order(trt_dat$gene_score), ]        # light-to-dark over-plotting
#
#   ggplot(trt_dat, aes(UMAP_1, UMAP_2, colour = gene_score)) +
#     geom_point(size = 0.7) +
#     scale_colour +
#     labs(colour = "Gene-set\nscore") +
#     ggtitle(trt) +
#     theme_void() +
#     theme(
#       plot.title   = element_text(hjust = 0.5, size = 14),
#       legend.title = element_text(size = 10),
#       legend.text  = element_text(size = 8),
#       legend.position = "right"
#     )
# })
#
# names(plot_list) <- treatments   # convenient for retrieval: plot_list[["WT"]]
#
# ## ---------- OPTIONAL: SHOW ALL TOGETHER -------------------------------------
# # Wrap into a single figure (3x2 grid; empty 6th slot suppressed)
# combined_plt <- wrap_plots(plot_list, nrow = 1, guides = "collect") &
#   theme(legend.position = "right")
# combined_plt         # print or save with ggsave()
#
# # Save the plot
# ggplot2::ggsave("../plots/DRGSNMI_SAM_signature_split_SN.pdf", plot = combined_plt, device = "pdf", width = 15.1, height = 3, dpi = 300)

message("06_main_immune_clustering complete.")
