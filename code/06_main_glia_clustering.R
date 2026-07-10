#' ======================================================================
#' 06_main_glia_clustering.R
#'
#' Subset, re-integrate, and cluster glial populations (115K cells) at res=0.9.
#'
#' Inputs:  ../r_objects/combined_object_init_named.RDS
#' Outputs: ../r_objects/glia_object_init.RDS, ../metadata/glia_*_celltypes.rds
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
###    GLIA SUBSET AND INTEGRATION      ###
###                                     ###
###########################################
###########################################

merged_named <- readRDS("../r_objects/combined_object_init_named.RDS")

glia <- subset(merged_named, idents = c("glia"))
# 115518 cells

options(future.globals.maxSize = 3e+09)
glia <- SCTransform(glia)
glia <- RunPCA(glia, npcs = 50)
glia <- IntegrateLayers(
  object = glia,
  method = RPCAIntegration,
  normalization.method = "SCT",
  k.weight = 80  # reduced from default 100; smallest glial batch has ~80 cells
)
glia <- FindNeighbors(glia, dims = 1:50, reduction = "integrated.dr")
glia <- FindClusters(glia, resolution = 0.5)
glia <- RunUMAP(glia, reduction = "integrated.dr", dims = 1:50)

# Plot UMAP, group by 'seurat_clusters'
p1 <- DimPlot(glia, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, raster = F) + NoLegend()

# Display plot
print(p1)

dir.create("../plots/umaps/glia", recursive = TRUE, showWarnings = FALSE)
ggsave("../plots/umaps/glia/glia_initial_umap.png", p1, width = 10, height = 10, units = "in", dpi = 300)

saveRDS(glia, "../r_objects/glia_object_init.RDS")

###########################################
###########################################
###                                     ###
###    CLUSTER RESOLUTION + CLUSTREE    ###
###                                     ###
###########################################
###########################################

glia <- readRDS("../r_objects/glia_object_init.RDS")

# Do not set RNA assay or normalize - FindClusters requires SCT assay

# Array of resolutions
resolutions <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

dir.create("../plots/clustree/glia", recursive = TRUE, showWarnings = FALSE)

for (res in resolutions) {
  # Find clusters and create plots
  glia <- FindClusters(glia, resolution = res)
  p <- DimPlot(glia, label = T, raster = F)

  # Save plot
  ggplot2::ggsave(paste0("../plots/clustree/glia/clustree_glia_", res, "res.png"), plot = p, device = "png", width = 10, height = 10, dpi = 300)
}

# Run clustree to examine clustering layout
p <- clustree(glia, prefix = "SCT_snn_res.") + theme(text = element_text(size = 20))

# Saving the plot with specified parameters
ggsave(filename = "../plots/clustree/glia/clustree_glia_summary_plot.png", plot = p, device = "png", width = 18, height = 9, dpi = 300)

# Resolution 0.9 selected via clustree

###########################################
###########################################
###                                     ###
###      GLIA MARKER GENES & SCORES     ###
###        0.9 RESOLUTION CHOSEN        ###
###           & INITIAL UMAPS           ###
###                                     ###
###########################################
###########################################

glia <- readRDS("../r_objects/glia_object_init.RDS")
glia <- FindClusters(glia, resolution = 0.9)
DefaultAssay(glia) <- "RNA"
glia <- NormalizeData(glia)
glia <- JoinLayers(glia)

# Create a directory for saving output files
# 20240626 GLIA FULL run with min.pct 0.05 and diff.pct -Inf
output_directory <- "../markers/glia_findmarkers_full"
dir.create(output_directory, showWarnings = FALSE)

# Identify and save marker genes for each cluster
# Metadata will print to console and be saved to a CSV file
find_and_save_all_markers(glia, output_directory, min.pct = 0.1, logfc.threshold = 0.25, summary = T)

# Plot initial umap plots
p1 <- DimPlot(glia, reduction = "umap", label = TRUE, repel = TRUE, raster = F) + NoLegend()
p2 <- DimPlot(glia, reduction = "umap", split.by = "study", raster = F)
p3 <- DimPlot(glia, reduction = "umap", split.by = "tissue", raster = F)

dir.create("../plots/umaps/glia", recursive = TRUE, showWarnings = FALSE)
ggsave(filename = "../plots/umaps/glia/glia_initial_umap.png", plot = p1, device = "png", width = 5, height = 5, dpi = 300)
ggsave(filename = "../plots/umaps/glia/glia_initial_umap_condition.png", plot = p2, device = "png", width = 40, height = 5, dpi = 300)
ggsave(filename = "../plots/umaps/glia/glia_initial_umap_tissue.png", plot = p3, device = "png", width = 10, height = 5, dpi = 300)

init_abundances <- table(glia@meta.data$study, Idents(glia))
write.csv(init_abundances, "../metadata/glia_initial_cellcounts.csv")

# Extract data
data <- FetchData(glia, vars = c("nFeature_RNA", "nCount_RNA"))
data$cluster <- Idents(glia)
average_values <- data %>%
  group_by(cluster) %>%
  summarise(
    Average_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
    Average_nCount_RNA = mean(nCount_RNA, na.rm = TRUE)
  )
# Export the result to a CSV file
write.csv(average_values, "../metadata/glia_avg_nFeaturenCount.csv", row.names = FALSE)

save_directory = "../plots/umaps/glia_highlights/"
dir.create(save_directory, recursive = TRUE, showWarnings = FALSE)
highlightAndSaveAllClusters(glia, save_directory = save_directory)

### Interactive marker exploration (run manually in RStudio) ###
# FeaturePlot(glia, features = c("Ptprz1", "Ctnnd2", "Slc4a4", "Entpd1"), order = T, raster = F)
# FeaturePlot(glia, features = c("Fabp7", "Ednrb", "Bgn", "Ogn"), order = T, raster = F)
# FeaturePlot(glia, features = c("Pdpn", "Bcan", "Cdh15", "Agt"), order = T, raster = F)
# FeaturePlot(glia, features = c("Gadd45b", "Nr4a3", "Nr4a2", "Ccn3"), order = T, raster = F)
# FeaturePlot(glia, features = c("Mrgprf", "Cldn11", "Entpd2", "Ngfr"), order = T, raster = F)
# FeaturePlot(glia, features = c("Ngfr", "Gfra3", "Gfap", "L1cam"), order = T, raster = F)
# FeaturePlot(glia, features = c("Fabp7", "Cdh15", "Morrbid", "Crym"), order = T, raster = F)
# FeaturePlot(glia, features = c("Heyl", "Hes5", "Aqp4", "Mlc1"), order = T, raster = F)
# FeaturePlot(glia, features = c("Ncmap", "Pllp", "Prx", "Cldn19"), order = T, raster = F)
# FeaturePlot(glia, features = c("Fxyd6", "Ogn", "Pmp2", "Ecscr"), order = T, raster = F)
# FeaturePlot(glia, features = c("Ngf", "Hcn1", "Tmcc3", "Gldn"), order = T, raster = F)
# FeaturePlot(glia, features = c("Sept3", "Ddn", "Epha6", "Emid1"), order = T, raster = F)
# FeaturePlot(glia, features = c("Scn7a", "Cadm2", "Ank3", "Ntrk3"), order = T, raster = F)
# FeaturePlot(glia, features = c("Adgb", "Glp2r", "Kcnq5", "Xkr4"), order = T, raster = F)
# FeaturePlot(glia, features = c("Aldh1a1", "Opcml", "Ephb2", "Etl4"), order = T, raster = F)
# FeaturePlot(glia, features = c("Scn7a", "Atf3", "Socs3", "Egr1"), order = T, raster = F)
# FeaturePlot(glia, features = c("Junb", "Fos", "Ier2", "Cxcl1"), order = T, raster = F)
# FeaturePlot(glia, features = c("Scn7a", "G0s2", "Cyr61", "Casp4"), order = T, raster = F)
# FeaturePlot(glia, features = c("Tll1", "Nhlh2", "Gucy1b2", "Hdac9"), order = T, raster = F)
# FeaturePlot(glia, features = c("Afm", "Trim30a", "Epha6", "Plcxd3"), order = T, raster = F)
# FeaturePlot(glia, features = c("Cdkn1c", "Csrp2", "Pou3f1", "Nkd1"), order = T, raster = F)
# FeaturePlot(glia, features = c("Dusp4", "Sh3gl3", "Nipal1", "Gnb2l1"), order = T, raster = F)
# FeaturePlot(glia, features = c("Aldh1a1", "Gfra3", "Fcgr2b", "C3"), order = T, raster = F)
# FeaturePlot(glia, features = c("Igfbp7", "Hrct1", "Rassf4", "C1ra"), order = T, raster = F)
# FeaturePlot(glia, features = c("Scn7a", "Klk8", "Mrgprf", "Ephx2"), order = T, raster = F)
# FeaturePlot(glia, features = c("Pde10a", "Pdgfc", "Sema3d", "Hspa1b"), order = T, raster = F)
# FeaturePlot(glia, features = c("Ephb1", "Kif26b", "Egr3", "Fap"), order = T, raster = F)
# FeaturePlot(glia, features = c("Shc3", "Atg4a", "Fbn2", "Mrtfa"), order = T, raster = F)
# FeaturePlot(glia, features = c("Kcnq5", "Adamts12", "Adam23", "Xylt1"), order = T, raster = F)
# FeaturePlot(glia, features = c("Ednrb", "Fos", "Jun", "Ier2"), order = T, raster = F)
# FeaturePlot(glia, features = c("Hspa1a", "Tnfaip6", "Atf3", "Socs3"), order = T, raster = F)
# FeaturePlot(glia, features = c("Lipg", "Cst3", "Ptgds", "Eva1a"), order = T, raster = F)
# FeaturePlot(glia, features = c("Gdpd2", "Fzd2", "Ntm", "Rgma"), order = T, raster = F)
# FeaturePlot(glia, features = c("Cxcl10", "Vcam1", "Isg15", "Irf1"), order = T, raster = F)
# FeaturePlot(glia, features = c("Relb", "Nfkbia", "Apoe", "Birc3"), order = T, raster = F)
# FeaturePlot(glia, features = c("Ccl2", "Gem", "Tnfaip3", "Nfkbie"), order = T, raster = F)
# FeaturePlot(glia, features = c("Scn7a", "Cxcl14", "Ucn2", "Btc"), order = T, raster = F)
# FeaturePlot(glia, features = c("Clcf1", "Tnfaip2", "Lgals3", "Hmga1"), order = T, raster = F)
# FeaturePlot(glia, features = c("Hmga2", "Timp1", "Gm32255", "Sdc1"), order = T, raster = F)
# FeaturePlot(glia, features = c("Tgfbi", "Tnc", "Ngfr", "Anxa1"), order = T, raster = F)

###########################################
###########################################
###                                     ###
###          GLIA CORRELATIONS          ###
###                                     ###
###########################################
###########################################

glia <- readRDS("../r_objects/glia_object_init.RDS")
glia <- FindClusters(glia, resolution = 0.9)
DefaultAssay(glia) <- "RNA"
glia <- NormalizeData(glia)
glia <- JoinLayers(glia)

# Extract the average gene expression for each cluster
avg_exp <- AverageExpression(glia)

# Restrict analysis to variable genes
glia <- FindVariableFeatures(object = glia)
variable_gene_names <- VariableFeatures(glia, assay = "RNA")

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
ggsave('../plots/correlations/glia_seuratclusters_correlation_matrix_hclust.pdf', plot=p, device='pdf', width=12, height=12, dpi=300)

###########################################
###########################################
###                                     ###
###     GLIA CLUSTER NAMING & UMAPS     ###
###                                     ###
###########################################
###########################################

glia <- readRDS("../r_objects/glia_object_init.RDS")
glia <- FindClusters(glia, resolution = 0.9)

DefaultAssay(glia) <- "RNA"
glia <- JoinLayers(glia)
glia <- NormalizeData(glia)

glia_named <- glia

# Minor-level naming
cluster_to_celltype_minor <- c(
  "0" = "nmSC",
  "1" = "SGC",
  "2" = "SGC",
  "3" = "mSC",
  "4" = "SGC",
  "5" = "SGC IEG",
  "6" = "SGC",
  "7" = "BATCH",
  "8" = "mSC IEG",
  "9" = "SGC",
  "10" = "nmSC IEG",
  "11" = "BATCH",
  "12" = "nmSC",
  "13" = "nmSC",
  "14" = "G_PROG",
  "15" = "SGC IFN",
  "16" = "mSC",
  "17" = "nmSC Ngfr+",
  "18" = "BATCH",
  "19" = "SGC",
  "20" = "G_CYC",
  "21" = "mSC",
  "22" = "mSC REP",
  "23" = "nmSC IEG",
  "24" = "BATCH",
  "25" = "BATCH"
)

glia_named <- RenameIdents(glia_named, cluster_to_celltype_minor)

minor_levels <- c(
  "SGC",
  "SGC IEG",
  "SGC IFN",
  "nmSC",
  "nmSC IEG",
  "nmSC Ngfr+",
  "mSC",
  "mSC IEG",
  "mSC REP",
  "G_CYC",
  "G_PROG",
  "BATCH"
)

Idents(glia_named) <- factor(Idents(glia_named), levels = minor_levels)

# Export minor cell type idents (to be used in main object)
glia_labels <- Idents(glia_named)
glia_cells <- Cells(glia_named)

# Export minor cell type idents and cells
saveRDS(glia_labels, "../metadata/glia_minor_celltypes.rds")
saveRDS(glia_cells, "../metadata/glia_minor_cells.rds")

print(DimPlot(glia_named, label=T))

#### ASIDE FOR MAJOR CELL TYPE EXPORTS ####

glia_named <- glia

# Major-level naming
cluster_to_celltype_major <- c(
  "0" = "nmSC",
  "1" = "SGC",
  "2" = "SGC",
  "3" = "mSC",
  "4" = "SGC",
  "5" = "SGC",
  "6" = "SGC",
  "7" = "BATCH",
  "8" = "mSC",
  "9" = "SGC",
  "10" = "nmSC",
  "11" = "BATCH",
  "12" = "nmSC",
  "13" = "nmSC",
  "14" = "G_CYC",
  "15" = "SGC",
  "16" = "mSC",
  "17" = "nmSC",
  "18" = "BATCH",
  "19" = "SGC",
  "20" = "G_CYC",
  "21" = "mSC",
  "22" = "mSC",
  "23" = "nmSC",
  "24" = "BATCH",
  "25" = "BATCH"
)

glia_named <- RenameIdents(glia_named, cluster_to_celltype_major)

major_levels <- c(
  "SGC",
  "nmSC",
  "mSC",
  "G_CYC",
  "BATCH"
)

Idents(glia_named) <- factor(Idents(glia_named), levels = major_levels)

# Export major cell type idents (to be used in main object)
glia_labels <- Idents(glia_named)
glia_cells <- Cells(glia_named)

# Export major cell type idents and cells
saveRDS(glia_labels, "../metadata/glia_major_celltypes.rds")
saveRDS(glia_cells, "../metadata/glia_major_cells.rds")

print(DimPlot(glia_named, label=T))

#### END ASIDE ####

# Re-apply minor-level idents for BATCH removal, plotting, and downstream export
glia_named <- glia
glia_named <- RenameIdents(glia_named, cluster_to_celltype_minor)
Idents(glia_named) <- factor(Idents(glia_named), levels = minor_levels)

# REMOVE background, technical artifacts
glia_named <- subset(glia_named, idents = c("BATCH"), invert = TRUE)

# Quick UMAP rerun after doublet removal
glia_named <- RunUMAP(glia_named, reduction = "integrated.dr", dims = 1:50)

# Minor cell type colors — matches the canonical glia_cell_type_colors in
# utils/variables.R and the current Figure 3 (SGC = purple, nmSC = yellow,
# mSC = blue, G_CYC / G_PROG = green). Retired clusters nmSC II / mSC III removed.
minor_cell_type_colors <- c(
  "SGC"        = "#b74ad3",
  "SGC IEG"    = "#7d4ad3",
  "SGC IFN"    = "#9901a3",
  "nmSC"       = "#ffff70",
  "nmSC IEG"   = "#f7f197",
  "nmSC Ngfr+" = "#ffc20a",
  "mSC"        = "#1a85ff",
  "mSC IEG"    = "#005FDB",
  "mSC REP"    = "#00C1FD",
  "G_CYC"      = "#00F7C5",
  "G_PROG"     = "#00D5AF"
)

glia_named <- reverse_idents(glia_named)

# Export cell type abundances
print(table(Idents(glia_named), glia_named@meta.data$tissue))
# (copied to Excel spreadsheet for sharing)

# Plot UMAP and split UMAP for quick visualization
# Figure panel exports found in figure_panel_export.R
# Plot using the minor cell type colors
p <- DimPlot(
  object = glia_named,
  reduction = "umap",
  cols = minor_cell_type_colors,
  raster = FALSE
) +
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  NoLegend()

psplit <- DimPlot(
  object = glia_named,
  reduction = "umap",
  cols = minor_cell_type_colors,
  raster = FALSE,
  split.by = "tissue"
) +
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  NoLegend()

# Set the plot label
p2 <- LabelClusters(
  plot = p,
  id = "ident",
  repel = TRUE,
  size = 2.7,
  fontface = "bold",
  colour = "black",
  box.padding = unit(0.35, "lines"),
  segment.color = "black",
  max.overlaps = Inf
)

p2split <- LabelClusters(
  plot = psplit,
  id = "ident",
  repel = TRUE,
  size = 2.7,
  fontface = "bold",
  colour = "black",
  box.padding = unit(0.35, "lines"),
  segment.color = "black",
  max.overlaps = Inf
)

print(p2)
print(p2split)

# Save glia object with minor-level idents for downstream analysis
saveRDS(glia_named, "../r_objects/DRGSNMI_glia.RDS")

### COMMENTED OUT: Markers, FeaturePlots, and GO enrichment ###
### These analyses are handled by dedicated downstream scripts ###
# merged_named <- readRDS("../r_objects/DRGSNMI_glia.RDS")
#
# # Create a directory for saving output files
# output_directory <- "../markers/glia_named_findmarkers_full"
# dir.create(output_directory, showWarnings = FALSE)
#
# # Identify and save marker genes for each cluster
# # Metadata will print to console and be saved to a CSV file
# find_and_save_all_markers(merged_named,
#                           output_directory,
#                           min.pct = -Inf,
#                           summary = T)
#
# markers <- c(
#   "Scn7a",
#   "Kcna1",
#   "Bcan",
#   "Pou3f1",
#   "Cdkn1c",
#   "Nkd1",
#   "Cxcl14",
#   "Ngfr",
#   "Shh",
#   "Ucn2",
#   "Mlc1",
#   "Cdh15",
#   "Fos",
#   "Egr1",
#   "Egr2",
#   "Hspa1a",
#   "Atf3",
#   "Gfap",
#   "Cst3",
#   "Crym",
#   "Fabp7",
#   "Ptprz1",
#   "Spp1",
#   "Mfap5",
#   "Ncmap",
#   "Gjb1",
#   "Gldn",
#   "Emid1",
#   "Gadd45b",
#   "Kcnk1",
#   "Tll1"
# )
#
# dir.create("../plots/umaps/glia", recursive = TRUE, showWarnings = FALSE)
#
# PlotGene <- function(obj, gene, output_directory = "../plots/umaps/glia/") {
#   p <- FeaturePlot(obj, features = gene, order = T, raster = F)
#   p <- p + theme(
#     plot.title = element_text(size = 20, face = "italic", hjust = 0.5),
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title.x = element_text(size = 10, hjust = 0),
#     axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
#   )
#   ggplot2::ggsave(paste0(output_directory, "glia_", gene, ".pdf"), plot = p, device = "pdf", width = 5, height = 5, dpi = 300)
#   ggplot2::ggsave(paste0(output_directory, "glia_", gene, ".png"), plot = p, device = "png", width = 5, height = 5, dpi = 300)
#   ggplot2::ggsave(paste0(output_directory, "glia_", gene, ".tiff"), plot = p, device = "tiff", width = 5, height = 5, dpi = 300)
# }
#
# for (marker in markers) {
#   PlotGene(merged_named, marker)
# }
#
# markers <- c(
#   "Mki67",
#   "Top2a",
#   "Plk1",
#   "Pou3f1",
#   "Nkd1",
#   "Sox4",
#   "Gap43",
#   "Ngfr",
#   "Nes"
# )
#
# for (marker in markers) {
#   PlotGene(glia, marker)
# }
#
# ###########################################
# ###########################################
# ###                                     ###
# ###     GLIA NAMED MARKER GO ENR.       ###
# ###                                     ###
# ###########################################
# ###########################################
#
# # Load glia object
# merged_named <- readRDS("../r_objects/DRGSNMI_glia.RDS")
#
# dir.create("../markers/glia_named_findmarkers_full/GO/", showWarnings = FALSE)
#
# # Process each cluster's DE genes for GO analysis
# unique_idents <- unique(Idents(merged_named))
# for (ident in unique_idents) {
#     # Load DE genes from CSV files
#     marker_file <- paste0('../markers/glia_named_findmarkers_full/FindMarkers_', ident, '.csv')
#     marker_data <- read.csv(marker_file, stringsAsFactors = FALSE,  row.names = 1)
#
#     # Extract up or down genes (you can set the thresholds)
#     upregulated_genes <- row.names(marker_data[marker_data$p_val_adj < 0.05 & marker_data$avg_log2FC > 1, ])
#     #downregulated_genes <- row.names(marker_data[marker_data$p_val_adj < 0.05 & marker_data$avg_log2FC < 1, ])
#     # Use all DGE genes
#     #de_genes <- row.names(marker_data[marker_data$p_val_adj < 0.05 & abs(marker_data$avg_log2FC) > 1, ])
#
#     # Up genes only
#     go_enrich_up <- enrichGO(
#       gene         = upregulated_genes,
#       OrgDb        = org.Mm.eg.db,
#       keyType      = 'SYMBOL',
#       ont          = 'BP' # Biological Process
#     )
#     # Save results
#     write.csv(go_enrich_up, file = paste0('../markers/glia_named_findmarkers_full/GO/GOBP_', ident, '_up.csv'))
#
# }
#
# # Initialize a named list to store depths for memoization
# depth_cache <- list()
#
# # Function to compute enrichment
# compute_enrichment <- function(gene_ratio, bg_ratio) {
#   # Convert "x/y" to numeric ratio
#   parse_ratio <- function(ratio_str) {
#     parts <- strsplit(ratio_str, "/")[[1]]
#     return(as.numeric(parts[1]) / as.numeric(parts[2]))
#   }
#
#   gene_ratio_num <- sapply(gene_ratio, parse_ratio)
#   bg_ratio_num <- sapply(bg_ratio, parse_ratio)
#
#   return(gene_ratio_num / bg_ratio_num)
# }
#
# # Function to compute depth of a GO term with memoization
# compute_go_depth <- function(go_id, ontology = "BP") {
#   # Check if depth is already computed
#   if (go_id %in% names(depth_cache)) {
#     return(depth_cache[[go_id]])
#   }
#
#   # Define root GO terms based on ontology
#   if (ontology == "BP") {
#     root_term <- "GO:0008150"
#     ancestors <- GOBPANCESTOR
#   } else if (ontology == "MF") {
#     root_term <- "GO:0003674"
#     ancestors <- GOMFANCESTOR
#   } else if (ontology == "CC") {
#     root_term <- "GO:0005575"
#     ancestors <- GOCCANCESTOR
#   } else {
#     stop("Invalid ontology. Choose from 'BP', 'MF', or 'CC'.")
#   }
#
#   # If the GO term is the root, depth is 1
#   if (go_id == root_term) {
#     depth_cache[[go_id]] <<- 1
#     return(1)
#   }
#
#   # Get parents
#   parents <- as.list(ancestors[[go_id]])
#
#   if (is.null(parents[[1]])) {
#     # If no parents found, assign NA
#     depth_cache[[go_id]] <<- NA
#     return(NA)
#   }
#
#   # Compute depths of all parents recursively
#   parent_depths <- sapply(parents[[1]], compute_go_depth, ontology = ontology)
#
#   if (all(is.na(parent_depths))) {
#     depth_cache[[go_id]] <<- NA
#     return(NA)
#   }
#
#   # Depth is 1 + minimum parent depth
#   depth <- 1 + min(parent_depths, na.rm = TRUE)
#
#   # Store in cache
#   depth_cache[[go_id]] <<- depth
#
#   return(depth)
# }
#
# # Process each cluster's DE genes for GO analysis
# unique_idents <- unique(Idents(merged_named))
#
# for (ident in unique_idents) {
#   # Load DE genes from CSV files
#   marker_file <- paste0('../markers/glia_named_findmarkers_full/FindMarkers_', ident, '.csv')
#
#   if (!file.exists(marker_file)) {
#     warning(paste("Marker file does not exist for cluster:", ident))
#     next
#   }
#
#   marker_data <- read.csv(marker_file, stringsAsFactors = FALSE, row.names = 1)
#
#   # Extract up or down genes (you can set the thresholds)
#   upregulated_genes <- row.names(marker_data[marker_data$p_val_adj < 0.05 & marker_data$avg_log2FC > 0 & marker_data$pct.1 > 0.100, ])
#   #downregulated_genes <- row.names(marker_data[marker_data$p_val_adj < 0.05 & marker_data$avg_log2FC < 1, ])
#   # Use all DGE genes
#   #de_genes <- row.names(marker_data[marker_data$p_val_adj < 0.05, ])
#
#   # Function to perform GO enrichment and enhance results
#   perform_go_enrichment <- function(gene_list, label, ontology = "BP") {
#     go_enrich <- enrichGO(
#       gene         = gene_list,
#       OrgDb        = org.Mm.eg.db,
#       keyType      = 'SYMBOL',
#       ont          = ontology, # Biological Process
#       pAdjustMethod = "BH",
#       qvalueCutoff  = 0.05,
#       readable      = TRUE
#     )
#
#     if (is.null(go_enrich) || nrow(go_enrich@result) == 0) {
#       warning(paste("No GO enrichment results for cluster:", ident, "(", label, ")"))
#       return(NULL)
#     }
#
#     # Convert enrichment result to data frame
#     go_enrich_df <- as.data.frame(go_enrich)
#
#     # Calculate enrichment
#     go_enrich_df$enrichment <- compute_enrichment(go_enrich_df$GeneRatio, go_enrich_df$BgRatio)
#
#     # Calculate depth
#     go_enrich_df$depth <- sapply(go_enrich_df$ID, compute_go_depth, ontology = ontology)
#
#     return(go_enrich_df)
#   }
#
#   # GO enrichment for all DE genes
#   #go_enrich_all <- perform_go_enrichment(de_genes, label = "all", ontology = "BP")
#
#   #if (!is.null(go_enrich_all)) {
#   #  # Save enhanced results
#   #  output_file_all <- paste0('../markers/glia_named_findmarkers_full/GO/GOBP_', ident, '_all.csv')
#   #  write.csv(go_enrich_all, file = output_file_all, row.names = FALSE)
#   #  message(paste("Saved GO enrichment results for cluster:", ident, " (all DE genes)"))
#   #}
#
#   # GO enrichment for upregulated genes only
#   go_enrich_up <- perform_go_enrichment(upregulated_genes, label = "upregulated", ontology = "BP")
#
#   if (!is.null(go_enrich_up)) {
#     # Save enhanced results
#     output_file_up <- paste0('../markers/glia_named_findmarkers_full/GO/GOBP_', ident, '_up.csv')
#     write.csv(go_enrich_up, file = output_file_up, row.names = FALSE)
#     message(paste("Saved GO enrichment results for cluster:", ident, " (upregulated genes)"))
#   }
# }
#
# ### Duplicate GO summary loop ###
# # unique_idents <- unique(Idents(merged_named))
# #
# # for (ident in unique_idents) {
# #   marker_file <- paste0('../markers/glia_named_findmarkers_full/FindMarkers_', ident, '.csv')
# #
# #   if (!file.exists(marker_file)) {
# #     warning(paste("Marker file does not exist for cluster:", ident))
# #     next
# #   }
# #
# #   marker_data <- read.csv(marker_file, stringsAsFactors = FALSE, row.names = 1)
# #
# #   upregulated_genes <- row.names(marker_data[marker_data$p_val_adj < 0.05 & marker_data$avg_log2FC > 0 & marker_data$pct.1 > 0.100, ])
# #   de_genes <- row.names(marker_data[marker_data$p_val_adj < 0.05, ])
# #
# #   print(paste0("Cluster: ", ident))
# #   print(paste0("Upregulated genes: ", length(upregulated_genes)))
# #   print(paste0("DE genes: ", length(de_genes)))
# #   print("")
# # }

message("06_main_glia_clustering complete.")
