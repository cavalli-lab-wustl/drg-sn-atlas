#' ======================================================================
#' kaminker_natcomm_2023_import.R
#'
#' Import Kaminker Nat Comms 2023 DRG nuclei dataset (10 samples). GEO: GSE201654.
#'
#' Inputs:  GEO count matrices
#' Outputs: ../../r_objects/batch_objects/kaminker_*.RDS
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
input_dir <- "../../raw_data/kaminker_natcomm_2023/GSE201654_RAW/"

# RDS Outputs
single_obj_dir <- "../../r_objects/single_objects/"
batch_obj_dir <- "../../r_objects/batch_objects/"


# File Names
files <- c(
  "GSM6069069_SAM24385747_raw_feature_bc_matrix.h5",
  "GSM6069070_SAM24385748_raw_feature_bc_matrix.h5",
  "GSM6069071_SAM24395829_raw_feature_bc_matrix.h5",
  "GSM6069072_SAM24395830_raw_feature_bc_matrix.h5",
  "GSM6069073_SAM24395831_raw_feature_bc_matrix.h5",
  "GSM6069074_SAM24395832_raw_feature_bc_matrix.h5",
  "GSM6069075_SAM24395833_raw_feature_bc_matrix.h5",
  "GSM6069076_SAM24395834_raw_feature_bc_matrix.h5",
  "GSM6069077_SAM24395835_raw_feature_bc_matrix.h5",
  "GSM6069078_SAM24395836_raw_feature_bc_matrix.h5"
)

conditions <- c(
  "mouseDRG_sample1",
  "mouseDRG_sample2",
  "mouseDRG_sample3",
  "mouseDRG_sample4",
  "mouseDRG_sample5",
  "mouseDRG_sample6",
  "mouseDRG_sample7",
  "mouseDRG_sample8",
  "mouseDRG_sample9",
  "mouseDRG_sample10"
)

simple_conditions <- c(
  "WT",
  "WT",
  "WT",
  "WT",
  "WT",
  "WT",
  "WT",
  "WT",
  "WT",
  "WT"
)

file_conditions <- setNames(conditions, files)
file_simple_conditions <- setNames(simple_conditions, files)

# Read each file and save Seurat object
# Add full project/sample name to each object
# Add percent.pt metadata to each object

# Initialize list to store Seurat objects for merging into study-level objects
kaminker_natcomm_2023_objects <- list()


# Loop through each file, read it, and save it as a Seurat object
for (file in files) {

  # Generate full sample name
  tissue <- "DRG"
  condition <- file_simple_conditions[file]
  condition2 <- file_conditions[file]
  method <- "10X"
  cellnuc <- "NUC"
  species <- "MM"
  investigator <- "Kaminker"
  id <- str_extract(file, "GSM\\d+")

  sample_name <- paste(tissue, condition, condition2, method, cellnuc, species, investigator, id, sep = "_")

  message("Processing ", file)
  message("Sample Name: ", sample_name)

  seurat_obj <- Read10X_h5(paste0(input_dir, file))

  # 500 nFeature_RNA - needs basic trim
  seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = sample_name, min.features = 500)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^mt-")

  # Save Seurat Object
  saveRDS(seurat_obj, file = paste0(single_obj_dir, sample_name, ".RDS"))

  # Add file to list
  kaminker_natcomm_2023_objects <- c(kaminker_natcomm_2023_objects, seurat_obj)
}

# Combine Seurat objects
kaminker_natcomm_2023_seurat <- merge(kaminker_natcomm_2023_objects[[1]], y = kaminker_natcomm_2023_objects[-1])


kaminker_natcomm_2023_seurat
head(kaminker_natcomm_2023_seurat@meta.data)
unique(kaminker_natcomm_2023_seurat@meta.data$orig.ident)



kaminker_natcomm_2023_seurat$age <- "Adult"

kaminker_natcomm_2023_seurat$study <- "kaminker_natcomm_2023"

kaminker_natcomm_2023_seurat$study_sex <- "Mixed"

kaminker_natcomm_2023_seurat$study_simple <- "Jung, 2023"

kaminker_natcomm_2023_seurat$tissue <- "DRG"

kaminker_natcomm_2023_seurat$treatment <- "WT"

kaminker_natcomm_2023_seurat$treatment_simple <- "WT"

kaminker_natcomm_2023_seurat@meta.data$tech <- "10X_RNA_V3_1"

kaminker_natcomm_2023_seurat@meta.data$tech_simple <- "10X"

kaminker_natcomm_2023_seurat@meta.data$cellnuc <- "NUC"





# Save merged Seurat objects
saveRDS(kaminker_natcomm_2023_seurat, file = paste0(batch_obj_dir, "kaminker_natcomm_2023_WT_seurat.RDS"))

message("kaminker_natcomm_2023 import complete.")

###############################################################################
### SOLO VERIFICATION (not part of pipeline — for independent QC exploration)
### Uncomment and run interactively to examine this dataset in isolation.
###############################################################################

# combined <- kaminker_natcomm_2023_seurat
# combined <- SCTransform(combined)
# combined <- RunPCA(combined, npcs = 50)
# combined <- RunUMAP(combined, reduction = "pca", dims = 1:50)
# combined <- FindNeighbors(combined, reduction = "pca", dims = 1:50)
# combined <- FindClusters(combined, resolution = 0.5)
# DimPlot(combined, reduction = "umap", label = TRUE, raster = FALSE) + NoLegend()
# DefaultAssay(combined) <- "RNA"
# combined <- NormalizeData(combined)
