#' ======================================================================
#' kuruvilla_cellrep_2022_import.R
#'
#' Import Kuruvilla Cell Reports 2022 DRG dataset (Drop-seq). GEO: GSE182099.
#'
#' Inputs:  GEO count matrices
#' Outputs: ../../r_objects/batch_objects/kuruvilla_*.RDS
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
input_dir <- "../../raw_data/kuruvilla_cellrep_2022/"

# RDS Outputs
batch_obj_dir <- "../../r_objects/batch_objects/"




# Kuruvilla et al. 2022 data is supplied in processed Seurat .RDS objects with comprehensive metadata
# All that is required for atlas prep is extraction of DRG data, addition of study metadata, and saving of the named Seurat object


combined_obj <- readRDS(paste0(input_dir, "complete-data-object.RDS"))
combined_obj

# Extract DRG data by "tissue" metadata slot
combined_obj <- subset(combined_obj, subset = tissue == "DRG")
combined_obj


combined_obj$orig.ident <- paste("DRG", combined_obj$orig.ident, "DS", "CELL", "MM", "Kuruvilla", sep = "_")

combined_obj@meta.data$age <- "Adult"
combined_obj@meta.data$study_sex <- "Mixed"
combined_obj@meta.data$study <- "kuruvilla_cellrep_2022"
combined_obj@meta.data$study_simple <- "Mapps, 2022"
#combined_obj@meta.data$tissue <- "DRG" # already present in source object, no need to overwrite
combined_obj@meta.data$treatment <- "WT"
combined_obj@meta.data$treatment_simple <- "WT"
combined_obj@meta.data$tech <- "DS_2_0"
combined_obj@meta.data$tech_simple <- "DS"
combined_obj@meta.data$cellnuc <- "CELL"


combined_obj
head(combined_obj@meta.data)
unique(combined_obj@meta.data$orig.ident)



# Save objects with study metadata and naming convention

saveRDS(combined_obj, file = paste0(batch_obj_dir, "kuruvilla_cellrep_2022_DRG_seurat.RDS"))


message("kuruvilla_cellrep_2022 import complete.")

###############################################################################
### SOLO VERIFICATION (not part of pipeline — for independent QC exploration)
### Uncomment and run interactively to examine this dataset in isolation.
### Pre-processed data — only plotting required.
###############################################################################

# DimPlot(combined_obj, reduction = "umap", label = TRUE, raster = FALSE) + NoLegend()
