#' ======================================================================
#' miller_re_natcomm_2023_import.R
#'
#' Import Miller Nat Comms 2023 DRG dataset (1 sample). GEO: GSE158892.
#'
#' Inputs:  GEO count matrices
#' Outputs: ../../r_objects/batch_objects/miller_re_*.RDS
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
input_dir <- "../../raw_data/miller_re_natcomm_2023/GSE198485_RAW/"

# RDS Outputs
batch_obj_dir <- "../../r_objects/batch_objects/"


# File Names
file <- "GSM5949409_filtered_feature_bc_matrix.h5"


# Generate full sample name
tissue <- "DRG"
condition <- "10X_DRG_18W"
method <- "10X"
cellnuc <- "CELL"
species <- "MM"
investigator <- "Miller_RE"
id <- str_extract(file, "GSM\\d+")

sample_name <- paste(tissue, condition, method, cellnuc, species, investigator, id, sep = "_")

message("Processing ", file)
message("Sample Name: ", sample_name)

# seurat_obj <- Read10X(paste0(input_dir, file))
# This dataset is published in .h5 format
seurat_obj <- Read10X_h5(paste0(input_dir, file), use.names = TRUE, unique.features = TRUE)

seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = sample_name, min.features = 500)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^mt-")


seurat_obj@meta.data$age <- "Adult"
seurat_obj@meta.data$study_sex <- "M"
seurat_obj@meta.data$study <- "miller_re_natcomm_2023"
seurat_obj@meta.data$study_simple <- "Obeidat, 2023"
seurat_obj@meta.data$tissue <- "DRG"
seurat_obj@meta.data$treatment <- "WT_C57"
seurat_obj@meta.data$treatment_simple <- "WT"
seurat_obj@meta.data$tech <- "10X_RNA_V3_1"
seurat_obj@meta.data$tech_simple <- "10X"
seurat_obj@meta.data$cellnuc <- "CELL"


seurat_obj
head(seurat_obj@meta.data)
unique(seurat_obj@meta.data$orig.ident)


# Save merged  Seurat objects
saveRDS(seurat_obj, file = paste0(batch_obj_dir, "miller_re_natcomm_2023_seurat.RDS"))

message("miller_re_natcomm_2023 import complete.")

###############################################################################
### SOLO VERIFICATION (not part of pipeline — for independent QC exploration)
### Uncomment and run interactively to examine this dataset in isolation.
###############################################################################

# combined <- seurat_obj
# options(future.globals.maxSize = 3e+09)
# combined <- SCTransform(combined)
# combined <- RunPCA(combined, npcs = 50)
# combined <- RunUMAP(combined, reduction = "pca", dims = 1:50)
# combined <- FindNeighbors(combined, reduction = "pca", dims = 1:50)
# combined <- FindClusters(combined, resolution = 0.5)
# DimPlot(combined, reduction = "umap", label = TRUE, raster = FALSE) + NoLegend()
# DefaultAssay(combined) <- "RNA"
# combined <- NormalizeData(combined)
