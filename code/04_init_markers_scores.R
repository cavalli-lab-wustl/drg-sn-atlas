#' ======================================================================
#' 04_init_markers_scores.R
#'
#' FindAllMarkers on initial clusters; QC violin plots and cell count exports.
#'
#' Inputs:  ../r_objects/combined_init_object.RDS
#' Outputs: ../markers/combined_init/ (CSVs), ../plots/qc_*.pdf
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
###    COMBINED MARKER GENES & SCORES   ###
###        0.9 RESOLUTION CHOSEN        ###
###           & INITIAL UMAPS           ###
###                                     ###
###########################################
###########################################

combined_object <- readRDS("../r_objects/combined_init_object.RDS")
# Resolution 0.9 selected in 02_init_clustree.R
combined_object <- FindClusters(combined_object, resolution = 0.9)
DefaultAssay(combined_object) <- "RNA"
combined_object <- NormalizeData(combined_object)
combined_object <- JoinLayers(combined_object)

# Create a directory for saving output files
output_directory <- "../markers/combined_init_findmarkers_full"
dir.create(output_directory)

# Identify and save marker genes for each cluster
# min.pct = 0.10: gene must be detected in >= 10% of cells in a cluster
# logfc.threshold = 0.5: require >= 0.5 log2FC (stricter than default 0.25)
find_and_save_all_markers(combined_object, output_directory, min.pct = 0.10, logfc.threshold = 0.5, summary = TRUE)

# Plot initial umap plots
dir.create("../plots/umaps/combined_init", recursive = TRUE, showWarnings = FALSE)
p1 <- DimPlot(combined_object, reduction = "umap", label = TRUE, raster = F) + NoLegend()
p2 <- DimPlot(combined_object, reduction = "umap", split.by = "study", raster = F)
p3 <- DimPlot(combined_object, reduction = "umap", split.by = "tissue", raster = F)

ggsave(filename = "../plots/umaps/combined_init/combined_initial_umap.png", plot = p1, device = "png", width = 5, height = 5, dpi = 300)
ggsave(filename = "../plots/umaps/combined_init/combined_initial_umap_condition.png", plot = p2, device = "png", width = 40, height = 5, dpi = 300)
ggsave(filename = "../plots/umaps/combined_init/combined_initial_umap_tissue.png", plot = p3, device = "png", width = 10, height = 5, dpi = 300)

# Create a directory for saving output files
output_directory <- "../metadata"
dir.create(output_directory)

init_abundances <- table(combined_object@meta.data$study, Idents(combined_object))
write.csv(init_abundances, "../metadata/combined_initial_cellcounts_study.csv")

init_abundances <- table(combined_object@meta.data$orig.ident, Idents(combined_object))
write.csv(init_abundances, "../metadata/combined_initial_cellcounts_batch.csv")

# Extract data
data <- FetchData(combined_object, vars = c("nFeature_RNA", "nCount_RNA"))

# Add cluster information
data$cluster <- Idents(combined_object)

# Calculate the average 'nFeature_RNA' and 'nCount_RNA' per cluster
average_values <- data %>%
  group_by(cluster) %>%
  summarise(
    Average_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
    Average_nCount_RNA = mean(nCount_RNA, na.rm = TRUE)
  )

# Export the result to a CSV file
write.csv(average_values, "../metadata/combined_avg_nFeaturenCount.csv", row.names = FALSE)

message("04_init_markers_scores complete.")
