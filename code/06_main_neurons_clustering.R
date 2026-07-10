#' ======================================================================
#' 06_main_neurons_clustering.R
#'
#' Subset DRG neurons (nFeature>3000), re-integrate, and cluster at res=0.8.
#'
#' Inputs:  ../r_objects/combined_object_init_named.RDS
#' Outputs: ../r_objects/neurons_object_init.RDS, ../metadata/neurons_*_celltypes*.rds
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

table(Idents(merged_named), merged_named@meta.data$study)

# Exclude SN-only studies (no DRG neurons) and low-quality neuron batches
studies_to_remove <- c("giger_elife_2022", "milbrandt_natneuro_2022", "suter_elife_2021", "horste_pnas_2019", "miller_fd_eneuro_2020")

neurons <- subset(merged_named, idents = c("neuron"))

# Subset by metadata Study slot
neurons <- subset(neurons, subset = study %in% studies_to_remove, invert = T)


# Subset by nFeature_RNA > 3000 — strict QC; neurons require high gene detection depth
neurons <- subset(neurons, subset = nFeature_RNA > 3000)

table(Idents(neurons), neurons@meta.data$study)


options(future.globals.maxSize = 3e+09)
neurons <- SCTransform(neurons)
neurons <- RunPCA(neurons, npcs = 50)
neurons <- IntegrateLayers(
  object = neurons,
  method = RPCAIntegration,
  normalization.method = "SCT",
  k.weight = 40  # reduced from default 100; smallest neuron batch has ~40 cells
)
neurons <- FindNeighbors(neurons, dims = 1:50, reduction = "integrated.dr")
neurons <- FindClusters(neurons, resolution = 0.5)
neurons <- RunUMAP(neurons, reduction = "integrated.dr", dims = 1:50)

# Plot UMAP, group by 'seurat_clusters'
p1 <- DimPlot(neurons, reduction = "umap", label = TRUE, repel = TRUE, raster = F) + NoLegend()
# Display plot
print(p1)
dir.create("../plots/umaps/neurons", recursive = TRUE, showWarnings = FALSE)
ggsave("../plots/umaps/neurons/neurons_initial_umap.png", p1, width = 10, height = 10, units = "in", dpi = 300)

saveRDS(neurons, "../r_objects/neurons_object_init.RDS")


###########################################
###########################################
###                                     ###
###    CLUSTER RESOLUTION + CLUSTREE    ###
###                                     ###
###########################################
###########################################

neurons <- readRDS("../r_objects/neurons_object_init.RDS")

# Do not set RNA assay or normalize - FindClusters requires SCT assay

# Array of resolutions
resolutions <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

dir.create("../plots/clustree/neurons", recursive = TRUE, showWarnings = FALSE)

for (res in resolutions) {
  # Find clusters and create plots
  neurons <- FindClusters(neurons, resolution = res)
  p <- DimPlot(neurons, label = T, raster = F)

  # Save plots with and without legend
  ggplot2::ggsave(paste0("../plots/clustree/neurons/clustree_neurons_", res, "res.png"), plot = p, device = "png", width = 10, height = 10, dpi = 300)
}

# Run clustree to examine clustering layout
p <- clustree(neurons, prefix = "SCT_snn_res.") + theme(text = element_text(size = 20))

# Saving the plot with specified parameters
ggsave(filename = "../plots/clustree/neurons/clustree_neurons_summary_plot.png", plot = p, device = "png", width = 18, height = 9, dpi = 300)


# Resolution 0.8 selected via clustree

neurons <- readRDS("../r_objects/neurons_object_init.RDS")
neurons <- FindClusters(neurons, resolution = 0.8)

DefaultAssay(neurons) <- "RNA"

neurons <- NormalizeData(neurons)
neurons <- JoinLayers(neurons)


# Create a directory for saving output files
output_directory <- "../markers/neurons_init_full"
dir.create(output_directory, showWarnings = FALSE)
find_and_save_all_markers(neurons, output_directory, min.pct = 0.1, logfc.threshold = 0.25, summary = T)


init_abundances <- table(neurons@meta.data$study, Idents(neurons))
write.csv(init_abundances, "../metadata/neurons_initial_cellcounts.csv")


# Extract data
data <- FetchData(neurons, vars = c("nFeature_RNA", "nCount_RNA"))
data$cluster <- Idents(neurons)
average_values <- data %>%
  group_by(cluster) %>%
  summarise(
    Average_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
    Average_nCount_RNA = mean(nCount_RNA, na.rm = TRUE)
  )
# Export the result to a CSV file
write.csv(average_values, "../metadata/neurons_avg_nFeaturenCount.csv", row.names = FALSE)


# Extract the average gene expression for each cluster
avg_exp <- AverageExpression(neurons)

# Restrict analysis to variable genes
neurons <- FindVariableFeatures(object = neurons)
variable_gene_names <- VariableFeatures(neurons, assay = "RNA")

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
ggsave('../plots/correlations/neurons_seuratclusters_correlation_matrix_hclust.pdf', plot=p, device='pdf', width=12, height=12, dpi=300)


### Interactive marker exploration (run manually in RStudio) ###
# FeaturePlot(neurons, features = c("Pvalb", "Trpc1", "Calca", "Ntrk2"), order = T)
# FeaturePlot(neurons, features = c("Calb1", "P2rx3", "Trpm8", "Sst"), order = T)
# FeaturePlot(neurons, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "Syt1"), order = T)
# FeaturePlot(neurons, features = c("Sstr2", "Trpm8", "Moxd1", "Mrgprd"), order = T)
# FeaturePlot(neurons, features = c("P2rx3", "Calca", "Ntrk2", "Th"), order = T)
# FeaturePlot(neurons, features = c("Th", "P2rx3"), order = T)
# FeaturePlot(neurons, features = c("Trpv1", "Trpm8"), order = T)
# FeaturePlot(neurons, features = c("Ntrk2", "Pvalb"), order = T)
# FeaturePlot(neurons, features = c("Calca", "Mrgprd"), order = T)
# FeaturePlot(neurons, features = c("Mrgprb4", "Mrgpra3"), order = T)
# FeaturePlot(neurons, features = c("Bmpr1b", "Sstr2"), order = T)
# FeaturePlot(neurons, features = c("Smr2", "Cysltr2"), order = T)
# FeaturePlot(neurons, features = c("Adra2a", "Ntrk3"), order = T)
# FeaturePlot(neurons, features = c("Ret", "Ntrk2"), order = T)
# FeaturePlot(neurons, features = c("Oprk1", "Tac1"), order = T)
# FeaturePlot(neurons, features = c("Calca", "Agt"), order = T)
# FeaturePlot(neurons, features = c("Insrr", "Adarb2"), order = T)
# FeaturePlot(neurons, features = c("Trpv1", "Pdyn"), order = T)
# FeaturePlot(neurons, features = c("Calb1", "Tnni1"), order = T)
# FeaturePlot(neurons, features = c("Calb1", "Avpr1a"), order = T)
# FeaturePlot(neurons, features = c("Calb1", "Ret"), order = T)
# FeaturePlot(neurons, features = c("Pth2r", "Crhr2"), order = T)
# FeaturePlot(neurons, features = c("Smr2", "Prokr2"), order = T)
# FeaturePlot(neurons, features = c("Crh", "Ecel1"), order = T)
# FeaturePlot(neurons, features = c("Atf3", "Trpm1"), order = T)
# FeaturePlot(neurons, features = c("Ucn", "Trbc2"), order = T)
# FeaturePlot(neurons, features = c("Ptprc", "Cd68"), order = T)
# FeaturePlot(neurons, features = c("Top2a", "nFeature_RNA"), order = T)
# FeaturePlot(neurons, features = c("Cckbr", "Nts"), order = T)
# FeaturePlot(neurons, features = c("Lpar3", "Dcn"), order = T)
# FeaturePlot(neurons, features = c("Avpr1a", "Sstr2"), order = T)
# FeaturePlot(neurons, features = c("Prokr2", "Asic3"), order = T)
# FeaturePlot(neurons, features = c("Gfra3", "S100b"), order = T)
# FeaturePlot(neurons, features = c("Asic3", "Smr2"), order = T)
# FeaturePlot(neurons, features = c("Ikzf1", "Ano1"), order = T)
# FeaturePlot(neurons, features = c("Ret", "Ntrk3"), order = T)
# FeaturePlot(neurons, features = c("Tmc2", "Nxph1"), order = T)
# FeaturePlot(neurons, features = c("Atf3", "Gal"), order = T)

dir.create("../plots/umaps/neurons_highlights", recursive = TRUE, showWarnings = FALSE)
highlightAndSaveAllClusters(neurons, save_directory = "../plots/umaps/neurons_highlights/")

bigdrg_markers <- c(
  "Runx3",
  "Pvalb",
  "Etv1",
  "Pcdh8",
  "Epha3",
  "Aldh1a2",
  "Cacna1h",
  "Skor2",
  "Cckar",
  "Tac1",
  "Plat",
  "Nsg2",
  "Lgi2",
  "Hapln4",
  "S100a16",
  "Kit",
  "Scgn",
  "Adra2c",
  "Slc18a3",
  "Chrna7",
  "Rxfp1",
  "Stum",
  "Foxp2",
  "Trpm8",
  "Sst",
  "Cck",
  "Adora2b",
  "Chrna3",
  "Adora1",
  "Cacng5",
  "Mrgprx1",
  "Mrgprx4",
  "Gfra2",
  "Cdh9",
  "Casq2",
  "Atf3",
  "Ucn",
  "Penk",
  "Trim54"
)

### Interactive DotPlot (run manually in RStudio) ###
# DotPlot(neurons, features = bigdrg_markers, dot.scale = 9, dot.min = 0.05)+
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1), # Rotate the x-axis labels
#     axis.title.x.top = element_blank(),
#     strip.text.x = element_blank(),
#     strip.text.y = element_blank(),
#   ) +
#   labs(x = NULL, y = NULL) +
#   guides(
#     color = guide_colorbar(title = "Avg. Exp.", label = FALSE),
#     size = guide_legend(title = "% Exp.")
#   )

neurons_named <- neurons

# Cluster order based on bigdrg_markers gene order (top to bottom in the marker list)
cluster_order <- c(
  # Pvalb (marker position 2)
  "10",
  # Etv1 (marker position 3)
  "29", "25",
  # Skor2/Lgi2 (marker position 8)
  "20",
  # Cckar (marker position 9)
  "2", "21",
  # S100a16 (marker position 15)
  "1",
  # Kit (marker position 16)
  "3", "18",
  # Rxfp1/Cacng5 (marker position 21)
  "22",
  # Stum (marker position 22)
  "24",
  # Trpm8/Foxp2 (marker positions 23-24)
  "7", "19",
  # Sst (marker position 25)
  "16",
  # Adora2b/Chrna3 (marker position 27)
  "11",
  # Adora1/Gfra2 (marker position 29)
  "4", "8", "12", "14",
  # Mrgprx1/Adora2b/Etv1 (Mrgprx1 at position 31)
  "17",
  # Ucn (marker position 37)
  "28",
  # BATCH clusters at end (including cluster 6 with no annotation)
  "0", "5", "6", "9", "13", "15", "23", "26", "27"
)

# Reorder the cluster factor levels
neurons_named$seurat_clusters <- factor(neurons_named$seurat_clusters, levels = cluster_order)
Idents(neurons_named) <- factor(Idents(neurons_named), levels = cluster_order)

neurons_named <- reverse_idents(neurons_named)

# Generate the DotPlot
### Interactive DotPlot (run manually in RStudio) ###
# DotPlot(neurons_named, features = bigdrg_markers) +
#   RotatedAxis()


neurons_named <- neurons

# Major-level naming
cluster_to_celltype <- c(
  "0" = "BATCH",
  "1" = "Asic3",
  "2" = "Th",
  "3" = "Sstr2",
  "4" = "P2rx3", #
  "5" = "BATCH",
  "6" = "Oprk1",
  "7" = "Trpm8",
  "8" = "P2rx3", #
  "9" = "BATCH",
  "10" = "Pvalb",
  "11" = "Adra2a",
  "12" = "P2rx3", #
  "13" = "BATCH",
  "14" = "P2rx3", #
  "15" = "BATCH",
  "16" = "Cysltr2",
  "17" = "Mrgpra3",
  "18" = "Sstr2",
  "19" = "Trpm8",
  "20" = "Ntrk2",
  "21" = "Th",
  "22" = "Trpv1",
  "23" = "BATCH",
  "24" = "Adra2a",
  "25" = "Mrgprb4",
  "26" = "BATCH",
  "27" = "BATCH",
  "28" = "Ecel1",
  "29" = "BATCH"
)

neurons_named <- RenameIdents(neurons_named, cluster_to_celltype)

major_levels <- c(
    "BATCH",
    "Th",
    "Asic3",
    "P2rx3",
    "Sstr2",
    "Trpv1",
    "Trpm8",
    "Cysltr2",
    "Pvalb",
    "Ntrk2",
    "Ecel1",
    "Adra2a",
    "Oprk1",
    "Mrgpra3",
    "Mrgprb4"
)

Idents(neurons_named) <- factor(Idents(neurons_named), levels = major_levels)


neurons_named <- neurons

# renaming all clusters
cluster_to_celltype <- c(
  "0" = "BATCH",
  "1" = "A",
  "2" = "C",
  "3" = "C",
  "4" = "C", #
  "5" = "BATCH",
  "6" = "C",
  "7" = "C",
  "8" = "C", #
  "9" = "BATCH",
  "10" = "A",
  "11" = "C",
  "12" = "C", #
  "13" = "BATCH",
  "14" = "C", #
  "15" = "BATCH",
  "16" = "C",
  "17" = "C",
  "18" = "C",
  "19" = "C",
  "20" = "A",
  "21" = "C",
  "22" = "C",
  "23" = "BATCH",
  "24" = "A",
  "25" = "C",
  "26" = "BATCH",
  "27" = "BATCH",
  "28" = "Atf3",
  "29" = "BATCH"
)

# NOTES
# Minor was done first
# Major is just A C or U


neurons_named <- RenameIdents(neurons_named, cluster_to_celltype)

major_levels <- c(
    "C",
    "A",
    "Atf3",
    "BATCH"
)

Idents(neurons_named) <- factor(Idents(neurons_named), levels = major_levels)


# Export idents (to be used in main object)
neurons_labels <- Idents(neurons_named)
neurons_cells <- Cells(neurons_named)


# Export major cell type idents and cells
saveRDS(neurons_labels, "../metadata/neurons_major_celltypes.rds")
saveRDS(neurons_cells, "../metadata/neurons_major_cells.rds")

# REMOVE background, technical artifacts
neurons_named <- subset(neurons_named, idents = c("BATCH"), invert = TRUE)

# Quick UMAP rerun after doublet removal
neurons_named <- RunUMAP(neurons_named, reduction = "integrated.dr", dims = 1:50)

# Plot UMAP
p1 <- DimPlot(neurons_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F) + NoLegend()
print(p1)
ggsave(filename = "../plots/umaps/neurons/neurons_umap_major_celltypes.png", plot = p1, device = "png", width = 10, height = 10, dpi = 300)

p3 <- DimPlot(neurons_named, reduction = "umap", label = F, repel = TRUE, raster = F) + NoLegend()
print(p3)
ggsave(filename = "../plots/umaps/neurons/neurons_umap_major_celltypes_unlabeled.png", plot = p3, device = "png", width = 5, height = 5, dpi = 300)

p4 <- DimPlot(neurons_named, reduction = "umap", label = F, repel = TRUE, raster = F, split.by = "tissue") + NoLegend()
print(p4)
ggsave(filename = "../plots/umaps/neurons/neurons_umap_major_celltypes_split_unlabeled.png", plot = p4, device = "png", width = 10, height = 5, dpi = 300)


#### FIGURE S1 EXPORTS

p_supp_major <- DimPlot(neurons_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F, split.by = "tissue") +
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
ggsave(filename = "../plots/s1_umaps/neurons_umap_major_celltypes_split.png", plot = p_supp_major, device = "png", width = 3.5, height = 4, dpi = 300)

ggsave(filename = "../plots/s1_umaps/neurons_umap_major_celltypes_split.pdf", plot = p_supp_major, device = "pdf", width = 3.5, height = 4, dpi = 300)


####### MINOR SUBDIVISIONS

neurons_named <- neurons

# Minor-level naming
cluster_to_celltype <- c(
  "0" = "BATCH",
  "1" = "A-PEP",
  "2" = "C-LTMR",
  "3" = "C-PEP.Sstr2",
  "4" = "C-NP.Mrgprd", #
  "5" = "BATCH",
  "6" = "C-PEP.Oprk1",
  "7" = "C-Thermo.Trpm8",
  "8" = "C-NP.Mrgprd", #
  "9" = "BATCH",
  "10" = "A-Propr.Pvalb",
  "11" = "C-PEP.Dcn",
  "12" = "C-NP.Mrgprd", #
  "13" = "BATCH",
  "14" = "C-NP.Mrgprd", #
  "15" = "BATCH",
  "16" = "C-NP.Sst",
  "17" = "C-NP.Mrgpra3",
  "18" = "C-PEP.Sstr2",
  "19" = "C-Thermo.Trpm8",
  "20" = "A-LTMR.Ntrk2",
  "21" = "C-LTMR",
  "22" = "C-Thermo.Rxfp1",
  "23" = "BATCH",
  "24" = "A-PEP",
  "25" = "C-NP.Mrgpra3",
  "26" = "BATCH",
  "27" = "BATCH",
  "28" = "U.Atf3",
  "29" = "BATCH"
)

# NOTES
# What are we missing from the BIGDRG Conversion Table (Fig. 6)?

# - Distinction between Trpm8/Rxfp1 Thermo [minor issue]

# - A-PEP.SCGN/ADRA2C aka Calca+Bmpr1b
# - A-PEP.KIT aka Calca+Bmpr1b
# - Above issues are/were solved by removing S100a16 distinction from A.PEP
# - Cluster 1 seems to have Kit/S100a16/Bmpr1b expression and likely represents all A.PEP

# - Distinctions between A/AB/AB/AB-LTMR subsets are generally limited
# - For now all just Ntrk2 A-LTMR (they are Ntrk2high Ntrk3low it seems)

# - No prop distinctions [typical in mouse]

# - seems like Cluster 1 (+Cluster 24) are A-PEP neurons
# - these are all collapsed right now, but based on plotting markers
# - there is definitely substructure to these clusters
# - for now, will leave collapsed as just A-PEP though
# - but check with Valeria before final version to see what she prefers


neurons_named <- RenameIdents(neurons_named, cluster_to_celltype)

major_levels <- c(
  "C-Thermo.Trpm8",
  "C-Thermo.Rxfp1",
  "C-LTMR",
  "C-NP.Mrgpra3",
  "C-NP.Mrgprd",
  "C-NP.Sst",
  "C-PEP.Sstr2",
  "C-PEP.Oprk1",
  "C-PEP.Dcn",
  "A-PEP",
  "A-Propr.Pvalb",
  "A-LTMR.Ntrk2",
  "U.Atf3",
  "BATCH"
)

Idents(neurons_named) <- factor(Idents(neurons_named), levels = major_levels)


# Export idents (to be used in main object)
neurons_labels <- Idents(neurons_named)
neurons_cells <- Cells(neurons_named)


# Export minor cell type idents and cells
saveRDS(neurons_labels, "../metadata/neurons_minor_celltypes.rds")
saveRDS(neurons_cells, "../metadata/neurons_minor_cells.rds")


# REMOVE background, technical artifacts
neurons_named <- subset(neurons_named, idents = c("BATCH"), invert = TRUE)

# Quick UMAP rerun after doublet removal
neurons_named <- RunUMAP(neurons_named, reduction = "integrated.dr", dims = 1:50)

# Plot UMAP
p1 <- DimPlot(neurons_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F) + NoLegend()
print(p1)
ggsave(filename = "../plots/umaps/neurons/neurons_umap_minor_celltypes.png", plot = p1, device = "png", width = 10, height = 10, dpi = 300)

p3 <- DimPlot(neurons_named, reduction = "umap", label = F, repel = TRUE, raster = F) + NoLegend()
print(p3)
ggsave(filename = "../plots/umaps/neurons/neurons_umap_minor_celltypes_unlabeled.png", plot = p3, device = "png", width = 5, height = 5, dpi = 300)

p4 <- DimPlot(neurons_named, reduction = "umap", label = F, repel = TRUE, raster = F, split.by = "tissue") + NoLegend()
print(p4)
ggsave(filename = "../plots/umaps/neurons/neurons_umap_minor_celltypes_split_unlabeled.png", plot = p4, device = "png", width = 10, height = 5, dpi = 300)


neuron_dotplot_genes <- c(
  "Trpm8",
  "Rxfp1",
  "Th",
  "Mrgpra3",
  #"Mrgprb4",
  "Mrgprd",
  "Sst",
  "Kit",
  "Sstr2",
  "Oprk1",
  "Dcn",
  "S100a16",
  #"Scgn",
  #"Smr2",
  #"Ntrk3",
  "Prokr2",
  "Bmpr1b",
  "Pvalb",
  "Ntrk2",
  "Atf3"
)

neurons_named <- reverse_idents(neurons_named)

pdot <- DotPlot(neurons_named, features = neuron_dotplot_genes, dot.scale = 9, dot.min = 0.05)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate the x-axis labels
    axis.title.x.top = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
  ) +
  labs(x = NULL, y = NULL) +
  guides(
    color = guide_colorbar(title = "Avg. Exp.", label = FALSE),
    size = guide_legend(title = "% Exp.")
  )

print(pdot)
ggsave(filename = "../plots/s1_umaps/neurons_markers_dotplot.pdf", plot = pdot, device = "pdf", width = 7, height = 4.5, dpi = 300)


#### FIGURE S1 EXPORTS

p_supp_minor <- DimPlot(neurons_named, reduction = "umap", label = TRUE, repel = TRUE, raster = F, split.by = "tissue") +
  labs(x = "UMAP-1", y = "UMAP-2") +
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  NoLegend()


print(p_supp_minor)
ggsave(filename = "../plots/s1_umaps/neurons_umap_minor_celltypes_split.png", plot = p_supp_minor, device = "png", width = 3.5, height = 4, dpi = 300)

ggsave(filename = "../plots/s1_umaps/neurons_umap_minor_celltypes_split.pdf", plot = p_supp_minor, device = "pdf", width = 3.5, height = 4, dpi = 300)


### ARCHIVE: Previous neuron subsetting approaches (commented out for review) ###
#
# # Assuming your Seurat object is named 'seurat_obj'
# study_counts <- table(neurons$Study)
# small_batches <- names(study_counts[study_counts < 10])
#
# print(study_counts)
# print(small_batches)
#
# # Define the studies you want to remove
# studies_to_remove <- c("giger_elife_2022", "milbrandt_natneuro_2022", "suter_elife_2021")
#
# # Remove these studies
# keep_cells <- !(neurons$Study %in% studies_to_remove)
# remove_cells <- (neurons$Study %in% studies_to_remove)
#
# # Subset by metadata Study slot
# neurons <- subset(neurons, subset = Study %in% small_batches, invert = T)
#
# #neurons <- neurons[, keep_cells]
#
# new_study_counts <- table(neurons$Study)
# print(new_study_counts)
#
#
# # Load necessary libraries
# library(Seurat)
# library(SeuratObject)
#
# # Load the merged Seurat object
# merged_named <- readRDS("../r_objects/combined_object_init_named.RDS")
#
# # Subset to include only "Neurons"
# neurons <- subset(merged_named, idents = "Neurons")
#
# # Verify the subset
# print("Number of cells after subsetting to Neurons:")
# print(dim(neurons))
#
# # Identify study counts
# study_counts <- table(neurons$Study)
# small_batches <- names(study_counts[study_counts < 10])
#
# print("Study Counts Before Removal:")
# print(study_counts)
# print("Small Batches to Remove:")
# print(small_batches)
#
# # Define the studies you want to remove
# studies_to_remove <- c("giger_elife_2022", "milbrandt_natneuro_2022", "suter_elife_2021")
#
# # Remove these studies using subset with 'invert = TRUE'
# neurons <- subset(neurons, subset = Study %in% small_batches, invert = T)
#
# # Verify the removal
# new_study_counts <- table(neurons$Study)
# print("Study Counts After Removal:")
# print(new_study_counts)
#
#
# # Load DietSeurat function if not already loaded
# # DietSeurat is part of Seurat package starting from version 3.2
# # Ensure your Seurat package is up to date within v5 constraints
#
# # Create a DietSeurat object retaining only the "RNA" assay and essential metadata
# # Specify assays to keep; here, only "RNA"
# # Remove other components like dimensional reductions, graphs, and misc data
#
# neurons_raw <- DietSeurat(
#   object = neurons,
#   assays = "RNA",
#   dimreducs = NULL,   # Remove dimensional reductions
#   graphs = NULL,      # Remove graph information
#   misc = NULL,        # Remove miscellaneous data
#   features = NULL,    # Keep all features in "RNA"
# )
#
# neurons_raw@meta.data$SCT_snn_res.0.5 <- NULL
# neurons_raw@meta.data$SCT_snn_res.1.1 <- NULL
# neurons_raw@meta.data$nCount_SCT <- NULL
# neurons_raw@meta.data$nFeature_SCT <- NULL
# neurons_raw@meta.data$seurat_clusters <- NULL
#
# # Verify the DietSeurat object
# print("Assays in DietSeurat object:")
# print(Assays(neurons_raw))
#
# print("Metadata columns:")
# print(colnames(neurons_raw@meta.data))
#
#
# # Set future globals max size
# options(future.globals.maxSize = 3e+09)
# neurons <- SCTransform(neurons)
# neurons <- RunPCA(neurons, npcs = 50)
# neurons <- IntegrateLayers(
#   object = neurons,
#   method = RPCAIntegration,
#   normalization.method = "SCT",
#   k.weight = 45
# )
# neurons <- FindNeighbors(neurons, dims = 1:50, reduction = "integrated.dr")
# neurons <- FindClusters(neurons, resolution = 0.5)
# neurons <- RunUMAP(neurons, reduction = "integrated.dr", dims = 1:50)
#
# # Plot UMAP, group by 'seurat_clusters'
# p1 <- DimPlot(neurons, reduction = "umap", label = TRUE, repel = TRUE, raster = F) + NoLegend()
# # Display plot
# p1
#
# saveRDS(neurons, "../r_objects/neurons_object_init.RDS")

message("06_main_neurons_clustering complete.")
