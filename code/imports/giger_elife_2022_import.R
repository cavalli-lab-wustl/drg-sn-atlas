#' ======================================================================
#' giger_elife_2022_import.R
#'
#' Import Giger eLife 2022 sciatic nerve dataset (4 timepoints). GEO: GSE198582.
#'
#' Inputs:  GEO count matrices
#' Outputs: ../../r_objects/batch_objects/giger_*.RDS
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
input_dir <- "../../raw_data/giger_elife_2022/"

# RDS Outputs
batch_obj_dir <- "../../r_objects/batch_objects/"

# Read each file and save Seurat object
# Add full project/sample name to each object
# Add percent.pt metadata to each object

# Giger et al. 2022 data is supplied in processed Seurat .RDS objects with comprehensive metadata
# All that is required for atlas prep is addition of study metadata and saving of the named Seurat objects

SN_0DPI <- readRDS(paste0(input_dir, "GSE198582_post_injury_day_0.RDS"))
SN_1DPI <- readRDS(paste0(input_dir, "GSE198582_post_injury_day_1.RDS"))
SN_3DPI <- readRDS(paste0(input_dir, "GSE198582_post_injury_day_3.RDS"))
SN_7DPI <- readRDS(paste0(input_dir, "GSE198582_post_injury_day_7.RDS"))

head(SN_0DPI@meta.data)

# note - these contain rich metadata, including cell type annotations
# can be used for comparison later (by barcode) but for now just add orig.ident and DRGSNMI-format metadata

# Generate full sample name
# This template is used to generate new orig.idents, but here we will modify the existing IDs
# commented out but kept as reference for filename structure
#tissue <- "SN"
#condition <- orig.ident
#method <- "10X"
#cellnuc <- "CELL"
#species <- "MM"
#investigator <- "Giger"
#sample_name <- paste(tissue, condition, method, cellnuc, species, investigator, sep = "_")

SN_0DPI$orig.ident <- paste("SN", SN_0DPI$orig.ident, "10X", "CELL", "MM", "Giger", sep = "_")
SN_1DPI$orig.ident <- paste("SN", SN_1DPI$orig.ident, "10X", "CELL", "MM", "Giger", sep = "_")
SN_3DPI$orig.ident <- paste("SN", SN_3DPI$orig.ident, "10X", "CELL", "MM", "Giger", sep = "_")
SN_7DPI$orig.ident <- paste("SN", SN_7DPI$orig.ident, "10X", "CELL", "MM", "Giger", sep = "_")

# add age, study_sex, study, study_simple, tissue, treatment, treatment_simple to each object

SN_0DPI$age <- "Adult"
SN_1DPI$age <- "Adult"
SN_3DPI$age <- "Adult"
SN_7DPI$age <- "Adult"

SN_0DPI$study <- "giger_elife_2022"
SN_1DPI$study <- "giger_elife_2022"
SN_3DPI$study <- "giger_elife_2022"
SN_7DPI$study <- "giger_elife_2022"

SN_0DPI$study_sex <- "Mixed"
SN_1DPI$study_sex <- "Mixed"
SN_3DPI$study_sex <- "Mixed"
SN_7DPI$study_sex <- "Mixed"

SN_0DPI$study_simple <- "Zhao, 2022"
SN_1DPI$study_simple <- "Zhao, 2022"
SN_3DPI$study_simple <- "Zhao, 2022"
SN_7DPI$study_simple <- "Zhao, 2022"

SN_0DPI$tissue <- "SN"
SN_1DPI$tissue <- "SN"
SN_3DPI$tissue <- "SN"
SN_7DPI$tissue <- "SN"

SN_0DPI$treatment <- "WT_SNC_0DPI"
SN_1DPI$treatment <- "SNC_1DPI"
SN_3DPI$treatment <- "SNC_3DPI"
SN_7DPI$treatment <- "SNC_7DPI"

SN_0DPI$treatment_simple <- "WT"
SN_1DPI$treatment_simple <- "SNC"
SN_3DPI$treatment_simple <- "SNC"
SN_7DPI$treatment_simple <- "SNC"

SN_0DPI@meta.data$tech <- "10X_RNA_V3_1"
SN_1DPI@meta.data$tech <- "10X_RNA_V3_1"
SN_3DPI@meta.data$tech <- "10X_RNA_V3_1"
SN_7DPI@meta.data$tech <- "10X_RNA_V3_1"

SN_0DPI@meta.data$tech_simple <- "10X"
SN_1DPI@meta.data$tech_simple <- "10X"
SN_3DPI@meta.data$tech_simple <- "10X"
SN_7DPI@meta.data$tech_simple <- "10X"

SN_0DPI@meta.data$cellnuc <- "CELL"
SN_1DPI@meta.data$cellnuc <- "CELL"
SN_3DPI@meta.data$cellnuc <- "CELL"
SN_7DPI@meta.data$cellnuc <- "CELL"


head(SN_0DPI@meta.data)
head(SN_1DPI@meta.data)
head(SN_3DPI@meta.data)
head(SN_7DPI@meta.data)

SN_0DPI
SN_1DPI
SN_3DPI
SN_7DPI

# Save objects with study metadata and naming convention

saveRDS(SN_0DPI, file = paste0(batch_obj_dir, "giger_elife_2022_SN_0DPI_seurat.RDS"))
saveRDS(SN_1DPI, file = paste0(batch_obj_dir, "giger_elife_2022_SN_1DPI_seurat.RDS"))
saveRDS(SN_3DPI, file = paste0(batch_obj_dir, "giger_elife_2022_SN_3DPI_seurat.RDS"))
saveRDS(SN_7DPI, file = paste0(batch_obj_dir, "giger_elife_2022_SN_7DPI_seurat.RDS"))

message("giger_elife_2022 import complete.")
