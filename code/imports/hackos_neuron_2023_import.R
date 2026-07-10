#' ======================================================================
#' hackos_neuron_2023_import.R
#'
#' Import Hackos Neuron 2023 DRG dataset (NaV1.7 variants). GEO: GSE213825.
#'
#' Inputs:  GEO count matrices
#' Outputs: ../../r_objects/batch_objects/hackos_*.RDS
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
input_dir <- "../../raw_data/hackos_neuron_2023/"


# RDS Outputs
single_obj_dir <- "../../r_objects/single_objects/"
batch_obj_dir <- "../../r_objects/batch_objects/"

# File Names
files <- c(
  "GSM6594872_SAM24349928_raw_feature_bc_matrix.h5",
  "GSM6594873_SAM24349929_raw_feature_bc_matrix.h5",
  "GSM6594874_SAM24352442_raw_feature_bc_matrix.h5",
  "GSM6594875_SAM24352443_raw_feature_bc_matrix.h5",
  "GSM6594876_SAM24352444_raw_feature_bc_matrix.h5",
  "GSM6594877_SAM24352445_raw_feature_bc_matrix.h5",
  "GSM6594878_SAM24353297_raw_feature_bc_matrix.h5",
  "GSM6594879_SAM24353298_raw_feature_bc_matrix.h5",
  "GSM6594880_SAM24353299_raw_feature_bc_matrix.h5",
  "GSM6594881_SAM24353300_raw_feature_bc_matrix.h5",
  "GSM6594882_SAM24374036_raw_feature_bc_matrix.h5",
  "GSM6594883_SAM24374037_raw_feature_bc_matrix.h5",
  "GSM6594884_SAM24374038_raw_feature_bc_matrix.h5",
  "GSM6594885_SAM24374039_raw_feature_bc_matrix.h5",
  "GSM6594886_SAM24374040_raw_feature_bc_matrix.h5",
  "GSM6594887_SAM24374041_raw_feature_bc_matrix.h5"
)

conditions <- c(
  "Scn9a_loxP_loxP_pCAAG_Cre_ER",
  "Scn9a_loxP_loxP_pCAAG_WT",
  "Scn9a_loxP_loxP_pCAAG_Cre_ER",
  "Scn9a_loxP_loxP_pCAAG_WT",
  "Scn9a_loxP_loxP_pCAAG_Cre_ER",
  "Scn9a_loxP_loxP_pCAAG_WT",
  "Scn9a_loxP_loxP_pCAAG_Cre_ER",
  "Scn9a_loxP_loxP_pCAAG_WT",
  "Scn9a_loxP_loxP_pCAAG_WT",
  "Scn9a_loxP_loxP_pCAAG_Cre_ER",
  "Scn9a_loxP_loxP_pCAAG_Cre_ER",
  "Scn9a_loxP_wt_Slc17a8_WT",
  "Scn9a_loxP_loxP_Slc17a8_Cre",
  "Scn9a_loxP_loxP_pCAAG_Cre_ER",
  "Scn9a_loxP_loxP_Slc17a8_WT",
  "Scn9a_loxP_loxP_Slc17a8_Cre"
)

simple_conditions <- c(
  "NaV1_7_KO",
  "NaV1_7_WT",
  "NaV1_7_KO",
  "NaV1_7_WT",
  "NaV1_7_KO",
  "NaV1_7_WT",
  "NaV1_7_KO",
  "NaV1_7_WT",
  "NaV1_7_WT",
  "NaV1_7_KO",
  "NaV1_7_KO",
  "NaV1_7_WT",
  "NaV1_7_cKO",
  "NaV1_7_KO",
  "NaV1_7_WT",
  "NaV1_7_cKO"
)

file_conditions <- setNames(conditions, files)
file_simple_conditions <- setNames(simple_conditions, files)

# Read each file and save Seurat object
# Add full project/sample name to each object
# Add percent.pt metadata to each object

# Initialize list to store Seurat objects for merging into study-level objects
NaV1_7_WT_objects <- list()
NaV1_7_KO_objects <- list()
NaV1_7_cKO_objects <- list()

# Loop through each file, read it, and save it as a Seurat object
for (file in files) {

  # Generate full sample name
  tissue <- "DRG"
  condition <- file_simple_conditions[file]
  condition2 <- file_conditions[file]
  method <- "10X"
  cellnuc <- "CELL"
  species <- "MM"
  investigator <- "Hackos"
  id <- str_extract(file, "GSM\\d+")

  sample_name <- paste(tissue, condition, condition2, method, cellnuc, species, investigator, id, sep = "_")

  message("Processing ", file)
  message("Sample Name: ", sample_name)

  seurat_obj <- Read10X_h5(paste0(input_dir, file))

  seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = sample_name, min.features = 500)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^mt-")

  # Save Seurat Object
  saveRDS(seurat_obj, file = paste0(single_obj_dir, sample_name, ".RDS"))

  # Add file to appropriate list
  if (condition == "NaV1_7_WT") {
    NaV1_7_WT_objects <- c(NaV1_7_WT_objects, seurat_obj)
  } else if (condition == "NaV1_7_KO") {
    NaV1_7_KO_objects <- c(NaV1_7_KO_objects, seurat_obj)
  } else if (condition == "NaV1_7_cKO") {
    NaV1_7_cKO_objects <- c(NaV1_7_cKO_objects, seurat_obj)
  } else {
    message("ERROR: Unknown condition: ", condition)
  }
}

# Combine Seurat objects
NaV1_7_WT_seurat <- merge(NaV1_7_WT_objects[[1]], y = NaV1_7_WT_objects[-1])
NaV1_7_KO_seurat <- merge(NaV1_7_KO_objects[[1]], y = NaV1_7_KO_objects[-1])
NaV1_7_cKO_seurat <- merge(NaV1_7_cKO_objects[[1]], y = NaV1_7_cKO_objects[-1])



NaV1_7_WT_seurat
NaV1_7_KO_seurat
NaV1_7_cKO_seurat

head(NaV1_7_WT_seurat@meta.data)
head(NaV1_7_KO_seurat@meta.data)
head(NaV1_7_cKO_seurat@meta.data)

# add age, study_sex, study, study_simple, tissue, treatment, treatment_simple to each object

NaV1_7_WT_seurat$age <- "Adult"
NaV1_7_KO_seurat$age <- "Adult"
NaV1_7_cKO_seurat$age <- "Adult"

NaV1_7_WT_seurat$study <- "hackos_neuron_2023"
NaV1_7_KO_seurat$study <- "hackos_neuron_2023"
NaV1_7_cKO_seurat$study <- "hackos_neuron_2023"

NaV1_7_WT_seurat$study_sex <- "Mixed"
NaV1_7_KO_seurat$study_sex <- "Mixed"
NaV1_7_cKO_seurat$study_sex <- "Mixed"

NaV1_7_WT_seurat$study_simple <- "Deng, 2023"
NaV1_7_KO_seurat$study_simple <- "Deng, 2023"
NaV1_7_cKO_seurat$study_simple <- "Deng, 2023"

NaV1_7_WT_seurat$tissue <- "DRG"
NaV1_7_KO_seurat$tissue <- "DRG"
NaV1_7_cKO_seurat$tissue <- "DRG"

NaV1_7_WT_seurat$treatment <- "WT_NaV1_7"
NaV1_7_KO_seurat$treatment <- "NaV1_7_KO"
NaV1_7_cKO_seurat$treatment <- "NaV1_7_cKO"

NaV1_7_WT_seurat$treatment_simple <- "WT"
NaV1_7_KO_seurat$treatment_simple <- "NaV1_KO"
NaV1_7_cKO_seurat$treatment_simple <- "NaV1_cKO"

NaV1_7_WT_seurat@meta.data$tech <- "10X_RNA_V3"
NaV1_7_KO_seurat@meta.data$tech <- "10X_RNA_V3"
NaV1_7_cKO_seurat@meta.data$tech <- "10X_RNA_V3"

NaV1_7_WT_seurat@meta.data$tech_simple <- "10X"
NaV1_7_KO_seurat@meta.data$tech_simple <- "10X"
NaV1_7_cKO_seurat@meta.data$tech_simple <- "10X"

NaV1_7_WT_seurat@meta.data$cellnuc <- "CELL"
NaV1_7_KO_seurat@meta.data$cellnuc <- "CELL"
NaV1_7_cKO_seurat@meta.data$cellnuc <- "CELL"


# Save merged  Seurat objects
saveRDS(NaV1_7_WT_seurat, file = paste0(batch_obj_dir, "hackos_neuron_2023_NaV1_7_WT_seurat.RDS"))
saveRDS(NaV1_7_KO_seurat, file = paste0(batch_obj_dir, "hackos_neuron_2023_NaV1_7_KO_seurat.RDS"))
saveRDS(NaV1_7_cKO_seurat, file = paste0(batch_obj_dir, "hackos_neuron_2023_NaV1_7_cKO_seurat.RDS"))

message("hackos_neuron_2023 import complete.")

###############################################################################
### SOLO VERIFICATION (not part of pipeline — for independent QC exploration)
### Uncomment and run interactively to examine this dataset in isolation.
###############################################################################

# combined <- NaV1_7_WT_seurat
# options(future.globals.maxSize = 3e+09)
# combined <- SCTransform(combined)
# combined <- RunPCA(combined, npcs = 50)
# combined <- RunUMAP(combined, reduction = "pca", dims = 1:50)
# combined <- FindNeighbors(combined, reduction = "pca", dims = 1:50)
# combined <- FindClusters(combined, resolution = 0.5)
# DimPlot(combined, reduction = "umap", label = TRUE, raster = FALSE) + NoLegend()
# DefaultAssay(combined) <- "RNA"
# combined <- NormalizeData(combined)
