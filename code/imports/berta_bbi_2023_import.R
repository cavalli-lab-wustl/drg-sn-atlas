#' ======================================================================
#' berta_bbi_2023_import.R
#'
#' Import Berta BBI 2023 DRG dataset (WT). GEO: GSE236914.
#'
#' Inputs:  GEO count matrices
#' Outputs: ../../r_objects/batch_objects/berta_*.RDS
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
source('../utils/libraries.R')
source('../utils/helper_functions.R')


# Directory Paths

# Input
input_dir <- "../../raw_data/berta_bbi_2023/"

# RDS Outputs
batch_obj_dir <- "../../r_objects/batch_objects/"




# Read each file and save Seurat object
# Add full project/sample name to each object
# Add percent.pt metadata to each object


file <- "GSE236914_RAW"

# Generate full sample name
tissue <- "DRG"
condition <- "WT"
condition2 <- "WT_CD1"
method <- "10X"
cellnuc <- "CELL"
species <- "MM"
investigator <- "Berta"
id <- "GSE236914"

sample_name <- paste(tissue, condition, condition2, method, cellnuc, species, investigator, id, sep = "_")


seurat_obj <- Read10X(paste0(input_dir, file))



# No min.features filter applied — published dataset assumed pre-filtered by authors
seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = sample_name)

# Expected: 27998 features across 13281 samples

# Inspect available metadata slots
head(seurat_obj@meta.data)

# Inspect orig.ident layer information
unique(seurat_obj@meta.data$orig.ident)

# This study has only one batch/library/layer, all WT

# Add metadata slots "age", "study_sex"
seurat_obj@meta.data$age <- "Adult"
seurat_obj@meta.data$study_sex <- "M"
seurat_obj@meta.data$study <- "berta_bbi_2023"
seurat_obj@meta.data$study_simple <- "Tonello, 2023"
seurat_obj@meta.data$tissue <- "DRG"
seurat_obj@meta.data$treatment <- "WT_CD1"
seurat_obj@meta.data$treatment_simple <- "WT"
seurat_obj@meta.data$tech <- "10X_RNA_V3"
seurat_obj@meta.data$tech_simple <- "10X"
seurat_obj@meta.data$cellnuc <- "CELL"


seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^mt-")

# Inspect updated metadata slots
head(seurat_obj@meta.data)


# Save Seurat Object
saveRDS(seurat_obj, file = paste0(batch_obj_dir, "berta_bbi_2023.RDS"))

message("berta_bbi_2023 import complete.")

###############################################################################
### SOLO VERIFICATION (not part of pipeline — for independent QC exploration)
### Uncomment and run interactively to examine this dataset in isolation.
###############################################################################

# # Load the batch object
# combined <- readRDS(paste0(batch_obj_dir, "berta_bbi_2023.RDS"))
#
# # SCTransform → PCA → UMAP → clustering
# combined <- SCTransform(combined)
# combined <- RunPCA(combined, npcs = 50)
# combined <- RunUMAP(combined, reduction = "pca", dims = 1:50)
# combined <- FindNeighbors(combined, reduction = "pca", dims = 1:50)
# combined <- FindClusters(combined, resolution = 0.9)
#
# # Plot UMAP
# p1 <- DimPlot(combined, reduction = "umap", label = TRUE, raster = FALSE) + NoLegend()
# print(p1)
#
# # Switch to RNA assay for markers
# DefaultAssay(combined) <- "RNA"
# combined <- NormalizeData(combined)
# combined <- JoinLayers(combined)
#
# # Find markers and QC metrics
# find_and_save_all_markers(combined, "solo_markers/", min.pct = 0.1, min.diff.pct = 0.1, summary = TRUE)
# write.csv(table(Idents(combined)), "solo_cellcounts.csv")
#
# # Cluster correlation heatmap
# avg_exp <- AverageExpression(combined)
# combined <- FindVariableFeatures(combined)
# avg_exp_var <- as.matrix(avg_exp$RNA[VariableFeatures(combined, assay = "RNA"), ])
# cor_matrix <- cor(avg_exp_var)
# hc <- hclust(as.dist(1 - cor_matrix), method = "average")
# cor_ordered <- cor_matrix[order.dendrogram(as.dendrogram(hc)),
#                           order.dendrogram(as.dendrogram(hc))]
# melted <- reshape2::melt(cor_ordered)
# ggplot(melted, aes(Var1, Var2, fill = value)) +
#   geom_tile() + scale_fill_viridis_c(limits = c(-1, 1)) +
#   theme_minimal() + coord_fixed()
#
# # Cluster-to-celltype mapping (res 0.9, 29 clusters)
# cluster_to_celltype <- c(
#   "0" = "neuron_lq_1", "1" = "sgc_1", "2" = "sgc_2", "3" = "msc",
#   "4" = "nmsc_1", "5" = "neuron_lq_2", "6" = "mural", "7" = "doublets",
#   "8" = "vec", "9" = "nmsc_2", "10" = "neuron_lq_3", "11" = "neuron_1",
#   "12" = "neuron_2", "13" = "neuron_3", "14" = "neuron_4", "15" = "neuron_lq_4",
#   "16" = "fibroblast_1", "17" = "sgc_3", "18" = "neuron_5", "19" = "immune_1",
#   "20" = "erythrocytes", "21" = "neuron_6", "22" = "neuron_7", "23" = "neuron_8",
#   "24" = "mitotic_1", "25" = "immune_2", "26" = "neuron_9",
#   "27" = "fibroblast_2", "28" = "neuron_10"
# )
#
# # Major cell type mapping
# major_mapping <- c(
#   "0" = "Neuron", "1" = "SGC", "2" = "SGC", "3" = "mSC",
#   "4" = "nmSC", "5" = "Neuron", "6" = "Mural", "7" = "Doublet",
#   "8" = "VEC", "9" = "nmSC", "10" = "Neuron", "11" = "Neuron",
#   "12" = "Neuron", "13" = "Neuron", "14" = "Neuron", "15" = "Neuron",
#   "16" = "Fibroblast", "17" = "SGC", "18" = "Neuron", "19" = "Immune",
#   "20" = "Erythrocyte", "21" = "Neuron", "22" = "Neuron", "23" = "Neuron",
#   "24" = "Mitotic", "25" = "Immune", "26" = "Neuron",
#   "27" = "Fibroblast", "28" = "Neuron"
# )
# combined <- RenameIdents(combined, major_mapping)
# combined <- subset(combined, idents = "Doublet", invert = TRUE)
# DimPlot(combined, label = TRUE, repel = TRUE, raster = FALSE) + NoLegend()
