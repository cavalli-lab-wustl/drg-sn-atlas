#' ======================================================================
#' milbrandt_natneuro_2022_import.R
#'
#' Import Milbrandt Nat Neurosci 2022 sciatic nerve nuclei dataset. GEO: GSE175421.
#'
#' Inputs:  GEO count matrices
#' Outputs: ../../r_objects/batch_objects/milbrandt_*.RDS
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
input_dir <- "../../raw_data/milbrandt_natneuro_2022/"

# RDS Outputs
batch_obj_dir <- "../../r_objects/batch_objects/"

# Milbrandt et al. 2022 data is supplied in processed Seurat .RDS objects with comprehensive metadata
# All that is required for atlas prep is addition of study metadata and saving of the named Seurat objects
#
# GEO provides 5 objects:
#   sn  = sciatic nerve (SN1/SN2/SN3 + SNC1/SNC2 — all WT, "SNC" = batch label not crush)
#   sc  = Schwann cell subset of sn
#   pn  = peripheral nerves (multiple nerve types)
#   pn_sc = Schwann cell subset of pn
#   ic  = immune cells (sorted, no UMAP embeddings)
#
# After examination (barcode overlap, cell counts, metadata):
#   - sc is a strict subset of sn (redundant)
#   - pn_sc is a strict subset of pn (redundant)
#   - ic has no UMAP embeddings, outside scope of DRGSNMI
#   - pn contains non-sciatic nerves, fewer SN cells than sn object
#   --> Only sn is used for the DRGSNMI atlas
#
# SN vs SNC in orig.ident: NOT crush vs control — represents batch (meta.data$tech column).
# All cells are WT. Cell counts match numbers reported in paper.

###############################################################################
### DATASET EXPLORATION (not part of pipeline — used to evaluate objects)
### Uncomment to load all 5 objects and compare barcodes, stats, etc.
###############################################################################

# sn <- readRDS(paste0(input_dir, "GSE182098_sciatic-nerve_single-nuclei-atlas-16Dec2020.RDS"))
# sc <- readRDS(paste0(input_dir, "GSE182098_schwann_cell_single-nuclei-atlas-12Jan2021.RDS"))
# pn <- readRDS(paste0(input_dir, "GSE182098_peripheral-nerves_single-nuclei-atlas-22Jan2021.RDS"))
# pn_sc <- readRDS(paste0(input_dir, "GSE182098_peripheral-nerves-Schwann-cell-specific_single-nuclei-atlas-8Feb2021.RDS"))
# ic <- readRDS(paste0(input_dir, "GSE182098_immune-cell_single-cell-atlas-15Jan2021.RDS"))
#
# # Utility: extract basic stats from a Seurat object
# get_seurat_stats <- function(obj, name) {
#   list(name = name, n_cells = ncol(obj), n_features = nrow(obj),
#        assays = names(obj@assays), cell_identities = levels(Idents(obj)),
#        metadata_cols = colnames(obj@meta.data))
# }
#
# # Check barcode overlaps between all 5 objects
# seurat_objects <- list(sn = sn, sc = sc, pn = pn, pn_sc = pn_sc, ic = ic)
# barcodes_list <- lapply(seurat_objects, colnames)
# overlap_matrix <- matrix(0, 5, 5, dimnames = list(names(seurat_objects), names(seurat_objects)))
# for (i in 1:5) for (j in 1:5) if (i != j)
#   overlap_matrix[i, j] <- length(intersect(barcodes_list[[i]], barcodes_list[[j]]))
# print(overlap_matrix)
#
# # Print stats for each object
# for (name in names(seurat_objects)) {
#   s <- get_seurat_stats(seurat_objects[[name]], name)
#   cat("\n=====", s$name, "=====\n")
#   cat("Cells:", s$n_cells, " Features:", s$n_features, "\n")
#   cat("Assays:", paste(s$assays, collapse = ", "), "\n")
# }

# Load only the sciatic nerve object for atlas import
sn <- readRDS(paste0(input_dir, "GSE182098_sciatic-nerve_single-nuclei-atlas-16Dec2020.RDS"))

sn
head(sn@meta.data)
unique(sn@meta.data$orig.ident)

sn$orig.ident <- paste("SN", sn$orig.ident, "10X", "NUC", "MM", "Milbrandt", sep = "_")

sn@meta.data$age <- "Adult"
sn@meta.data$study_sex <- "Mixed"
sn@meta.data$study <- "milbrandt_natneuro_2022"
sn@meta.data$study_simple <- "Yim, 2022"
sn@meta.data$tissue <- "SN"
sn@meta.data$treatment <- "WT"
sn@meta.data$treatment_simple <- "WT"

sn@meta.data$tech <- "10X_RNA_V2"
sn@meta.data$tech_simple <- "10X"
sn@meta.data$cellnuc <- "NUC"


sn
head(sn@meta.data)
unique(sn@meta.data$orig.ident)

saveRDS(sn, file = paste0(batch_obj_dir, "milbrandt_natneuro_2022_sciatic_nerve_seurat.RDS"))


message("milbrandt_natneuro_2022 import complete.")

###############################################################################
### SOLO VERIFICATION (not part of pipeline — for independent QC exploration)
### Uncomment and run interactively. Pre-processed data — only plotting required.
### Requires loading all 5 objects (see DATASET EXPLORATION block above).
###############################################################################

# DimPlot(sn, reduction = "umap", group.by = "res.0.6", label = TRUE, raster = FALSE) + NoLegend()
# DimPlot(sc, reduction = "umap", group.by = "res.0.6", label = TRUE, raster = FALSE) + NoLegend()
# DimPlot(pn, reduction = "umap", group.by = "integrated_snn_res.0.5", label = TRUE, raster = FALSE) + NoLegend()
# DimPlot(pn_sc, reduction = "umap", group.by = "integrated_snn_res.0.5", label = TRUE, raster = FALSE) + NoLegend()
# DimPlot(ic, reduction = "umap", group.by = "res.0.8", label = TRUE, raster = FALSE) + NoLegend()
