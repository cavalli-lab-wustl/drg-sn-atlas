#' ======================================================================
#' suter_elife_2021_import.R
#'
#' Import Suter eLife 2021 sciatic nerve dataset (P1+P60). GEO: GSE142541.
#'
#' Inputs:  GEO count matrices
#' Outputs: ../../r_objects/batch_objects/suter_*.RDS
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
input_dir <- "../../raw_data/suter_elife_2021/GSE138577_RAW/"

# RDS Outputs
batch_obj_dir <- "../../r_objects/batch_objects/"

# File Names
files <- c(
    "GSM4113877",
    "GSM4113878",
    "GSM4113879",
    "GSM4113880",
    "GSM4113881",
    "GSM4113882"
)

conditions <- c(
    "10X_P1_1",
    "10X_P1_2",
    "10X_P1_3",
    "10X_P60_1",
    "10X_P60_2",
    "10X_P60_3"
)

simple_conditions <- c(
    "P1",
    "P1",
    "P1",
    "P60",
    "P60",
    "P60"
)

file_conditions <- setNames(conditions, files)
file_simple_conditions <- setNames(simple_conditions, files)

# Read each file and save Seurat object
# Add full project/sample name to each object
# Add percent.pt metadata to each object

# Initialize list to store Seurat objects for merging into study-level objects
P1_objects <- list()
P60_objects <- list()

# Loop through each file, read it, and save it as a Seurat object
for (file in files) {

  # Generate full sample name
  tissue <- "SN"
  condition <- file_conditions[file]
  method <- "10X"
  cellnuc <- "CELL"
  species <- "MM"
  investigator <- "Suter"
  id <- str_extract(file, "GSM\\d+")

  sample_name <- paste(tissue, condition, method, cellnuc, species, investigator, id, sep = "_")

  message("Processing ", file)
  message("Sample Name: ", sample_name)

  seurat_obj <- Read10X(paste0(input_dir, file))

  seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = sample_name, min.features = 500)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^mt-")

  simple_condition <- file_simple_conditions[file]

  # Add file to appropriate list
    if (simple_condition == "P1") {
        P1_objects <- c(P1_objects, seurat_obj)
    } else if (simple_condition == "P60") {
        P60_objects <- c(P60_objects, seurat_obj)
    } else {
        message("ERROR: Unknown condition: ", simple_condition)
    }
}

# Combine Seurat objects
P1_seurat <- merge(P1_objects[[1]], y = P1_objects[-1])
P60_seurat <- merge(P60_objects[[1]], y = P60_objects[-1])

P1_seurat
P60_seurat

P1_seurat@meta.data$age <- "P1"
P1_seurat@meta.data$study_sex <- "Mixed"
P1_seurat@meta.data$study <- "suter_elife_2021"
P1_seurat@meta.data$study_simple <- "Gerber, 2021"
P1_seurat@meta.data$tissue <- "SN"
P1_seurat@meta.data$treatment <- "WT_0Cre_ReYFPTg_bActDsRedTg"
P1_seurat@meta.data$treatment_simple <- "WT"

P1_seurat@meta.data$tech <- "10X_RNA_V2"
P1_seurat@meta.data$tech_simple <- "10X"
P1_seurat@meta.data$cellnuc <- "CELL"

P1_seurat
head(P1_seurat@meta.data)
unique(P1_seurat@meta.data$orig.ident)

P60_seurat@meta.data$age <- "P60"
P60_seurat@meta.data$study_sex <- "Mixed"
P60_seurat@meta.data$study <- "suter_elife_2021"
P60_seurat@meta.data$study_simple <- "Gerber, 2021"
P60_seurat@meta.data$tissue <- "SN"
P60_seurat@meta.data$treatment <- "WT_0Cre_ReYFPTg_bActDsRedTg"
P60_seurat@meta.data$treatment_simple <- "WT"

P60_seurat@meta.data$tech <- "10X_RNA_V2"
P60_seurat@meta.data$tech_simple <- "10X"
P60_seurat@meta.data$cellnuc <- "CELL"

P60_seurat
head(P60_seurat@meta.data)
unique(P60_seurat@meta.data$orig.ident)

# Save merged  Seurat objects
saveRDS(P1_seurat, file = paste0(batch_obj_dir, "suter_elife_2021_P1_seurat.RDS"))
saveRDS(P60_seurat, file = paste0(batch_obj_dir, "suter_elife_2021_P60_seurat.RDS"))

message("suter_elife_2021 import complete.")

###############################################################################
### SOLO VERIFICATION (not part of pipeline — for independent QC exploration)
### Uncomment and run interactively to examine this dataset in isolation.
###############################################################################

# combined <- P60_seurat
# options(future.globals.maxSize = 3e+09)
# combined <- SCTransform(combined)
# combined <- RunPCA(combined, npcs = 50)
# combined <- RunUMAP(combined, reduction = "pca", dims = 1:50)
# combined <- FindNeighbors(combined, reduction = "pca", dims = 1:50)
# combined <- FindClusters(combined, resolution = 0.5)
# DimPlot(combined, reduction = "umap", label = TRUE, raster = FALSE) + NoLegend()
# DefaultAssay(combined) <- "RNA"
# combined <- NormalizeData(combined)
