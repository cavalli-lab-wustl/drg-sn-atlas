#' ======================================================================
#' 20_main_qc_export.R
#'
#' Final QC assessment; export mixed identity annotations and cell metadata.
#'
#' Inputs:  ../r_objects/DRGSNMI_main_multilevel.RDS
#' Outputs: ../metadata/*.csv
#'
#' Author: Michael B. Thomsen, Ph.D.
#' Created: 2024-07-01
#' ======================================================================


# UPDATED 20251110 with new object (IMG and MLC annotations)
seurat_obj <- readRDS("../r_objects/DRGSNMI_main_multilevel.RDS")

# Create a mixed identity column: cluster level for neurons, minor level for others
seurat_obj$mixed_ident <- ifelse(
  seurat_obj$Major %in% c("NEU_A", "NEU_C", "NEU_ATF3"),
  as.character(seurat_obj$Minor),
  as.character(seurat_obj$Major)
)

Idents(seurat_obj) <- "mixed_ident"

DimPlot(seurat_obj)
unique(Idents(seurat_obj))

# Define logical order: Neurons -> Glia -> Fibro -> Vasc -> Other
ordered_levels <- c(
  # Neurons (organized by functional groups)
  "C-Thermo.Trpm8",
  "C-Thermo.Rxfp1",
  "C-LTMR",
  "C-NP.Mrgpra3",
  "C-NP.Mrgprd",
  "C-NP.Sst",
  "C-PEP.Sstr2",
  "C-PEP.Oprk1",
  "C-PEP.Dcn",
  "A-PEP",
  "A-Propr.Pvalb",
  "A-LTMR.Ntrk2",
  "U.Atf3",

  # Glia
  "SGC",                                               # Satellite glia
  "NMSC",                                              # Non-myelinating Schwann
  "MYSC",                                              # Myelinating Schwann
  "GSC",                                               # Glial stem

  # Fibroblasts/Connective
  "EF",                                                #
  "EM",                                                #
  "PM",                                                #
  "FSC",                                               #

  # Vascular
  "VEC_A",                                              # Arterial endothelial
  "VEC_V",                                              # Venous endothelial
  "LEC",                                               # Lymphatic endothelial
  "PERI",                                              # Pericytes
  "VSMC",                                              # Vascular smooth muscle
  "MSC",
  "VSC",

  # Other/Immune
  "MAC", "IMG", "MONO",                         # Macrophage lineage
  "DC",                                                # Dendritic cells
  "LEUK",                                              # Leukocytes (general)
  "MAST",                                              # Mast cells
  "EOS",                                               # Eosinophils
  "ERY",                                               # Erythrocytes
  "NT", "HSC"                           # Stem/other
)

Idents(seurat_obj) <- factor(Idents(seurat_obj), levels = ordered_levels)

p1 <- VlnPlot(seurat_obj, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"), pt.size = 0)
p3 <- VlnPlot(seurat_obj, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"), group.by = "study_simple", pt.size = 0)

ggsave("../plots/qc/qc_vlnplots_celltype.pdf", plot = p1, width = 12, height = 3.5, dpi = 300)
ggsave("../plots/qc/qc_vlnplots_study_simple.pdf", plot = p3, width = 5, height = 3.5, dpi = 300)

# PERCENT MT EXPORT

p1 <- VlnPlot(seurat_obj, features = c("percent.mt"), pt.size = 0) + NoLegend()
p2 <- VlnPlot(seurat_obj, features = c("percent.mt"), group.by = "study_simple", pt.size = 0) + NoLegend()

ggsave("../plots/qc/qc_vlnplot_percentmt_celltype.pdf", plot = p1, width = 12, height = 3.5, dpi = 300)
ggsave("../plots/qc/qc_vlnplot_percentmt_study_simple.pdf", plot = p2, width = 5, height = 5, dpi = 300)


# NFEATURE EXPORT

p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA"), pt.size = 0) + NoLegend()
p2 <- VlnPlot(seurat_obj, features = c("nFeature_RNA"), group.by = "study_simple", pt.size = 0) + NoLegend()

ggsave("../plots/qc/qc_vlnplot_nFeature_RNA_celltype.pdf", plot = p1, width = 12, height = 3.5, dpi = 300)
ggsave("../plots/qc/qc_vlnplot_nFeature_RNA_study_simple.pdf", plot = p2, width = 5, height = 5, dpi = 300)


# NCOUNT EXPORT

p1 <- VlnPlot(seurat_obj, features = c("nCount_RNA"), pt.size = 0) + NoLegend()
p2 <- VlnPlot(seurat_obj, features = c("nCount_RNA"), group.by = "study_simple", pt.size = 0) + NoLegend()

ggsave("../plots/qc/qc_vlnplot_nCount_RNA_celltype.pdf", plot = p1, width = 12, height = 3.5, dpi = 300)
ggsave("../plots/qc/qc_vlnplot_nCount_RNA_study_simple.pdf", plot = p2, width = 5, height = 5, dpi = 300)

### MAJOR-LEVEL QC EXPORT ###

# UPDATED 20260115 with new object, rerun to export data with correct Major/Minor divisions
seurat_obj <- readRDS("../r_objects/DRGSNMI_main_multilevel.RDS")

# Define logical order: Neurons -> Glia -> Fibro -> Vasc -> Other
ordered_levels <- c(
  # Neurons (organized by functional groups)

  "NEU_A",
  "NEU_C",
  "NEU_ATF3",

  # Glia
  "SGC",                                               # Satellite glia
  "NMSC",                                              # Non-myelinating Schwann
  "MYSC",                                              # Myelinating Schwann
  "GSC",                                               # Glial stem

  # Fibroblasts/Connective
  "EF",                                                #
  "EM",                                                #
  "PM",                                                #
  "FSC",                                               #

  # Vascular
  "VEC_A",                                              # Arterial endothelial
  "VEC_V",                                              # Venous endothelial
  "LEC",                                               # Lymphatic endothelial
  "PERI",                                              # Pericytes
  "VSMC",                                              # Vascular smooth muscle
  "MSC",
  "VSC",

  # Other/Immune
  "MAC", "IMG", "MONO",                         # Macrophage lineage
  "DC",                                                # Dendritic cells
  "LEUK",                                              # Leukocytes (general)
  "MAST",                                              # Mast cells
  "EOS",                                               # Eosinophils
  "ERY",                                               # Erythrocytes
  "NT", "HSC"                           # Stem/other
)


Idents(seurat_obj) <- factor(Idents(seurat_obj), levels = ordered_levels)



DimPlot(seurat_obj)

# Load required libraries
library(Seurat)
library(openxlsx)
library(dplyr)
library(tidyr)

# ------------------------------
# Functions for QC summary statistics
# ------------------------------

# Function to calculate cell counts by cluster and a given grouping variable (e.g., condition or orig.ident)
calculate_cell_counts <- function(seurat_obj, group_var = "study_simple") {
  data.frame(
    Cluster = Idents(seurat_obj),
    Group = seurat_obj[[group_var]][, 1]
  ) %>%
    count(Cluster, Group) %>%
    pivot_wider(names_from = Group, values_from = n, values_fill = 0)
}

# Function to calculate detailed summary statistics by a given grouping variable (e.g., condition or orig.ident)
calculate_metrics_summary <- function(seurat_obj, group_var = "study_simple") {
  seurat_obj@meta.data %>%
    group_by_at(group_var) %>%
    summarize(
      count = n(),
      ## nFeature_RNA stats
      mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
      median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
      sd_nFeature_RNA = sd(nFeature_RNA, na.rm = TRUE),
      min_nFeature_RNA = min(nFeature_RNA, na.rm = TRUE),
      q1_nFeature_RNA = quantile(nFeature_RNA, 0.25, na.rm = TRUE),
      q3_nFeature_RNA = quantile(nFeature_RNA, 0.75, na.rm = TRUE),
      max_nFeature_RNA = max(nFeature_RNA, na.rm = TRUE),
      ## nCount_RNA stats
      mean_nCount_RNA = mean(nCount_RNA, na.rm = TRUE),
      median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
      sd_nCount_RNA = sd(nCount_RNA, na.rm = TRUE),
      min_nCount_RNA = min(nCount_RNA, na.rm = TRUE),
      q1_nCount_RNA = quantile(nCount_RNA, 0.25, na.rm = TRUE),
      q3_nCount_RNA = quantile(nCount_RNA, 0.75, na.rm = TRUE),
      max_nCount_RNA = max(nCount_RNA, na.rm = TRUE),
      ## percent.mt stats
      mean_percent_mt = mean(percent.mt, na.rm = TRUE),
      median_percent_mt = median(percent.mt, na.rm = TRUE),
      sd_percent_mt = sd(percent.mt, na.rm = TRUE),
      min_percent_mt = min(percent.mt, na.rm = TRUE),
      q1_percent_mt = quantile(percent.mt, 0.25, na.rm = TRUE),
      q3_percent_mt = quantile(percent.mt, 0.75, na.rm = TRUE),
      max_percent_mt = max(percent.mt, na.rm = TRUE)
    ) %>%
    ungroup()
}

# Function to calculate detailed summary statistics by cluster (for QC metrics)
calculate_metrics_by_cluster <- function(seurat_obj) {
  data.frame(
    Cluster = Idents(seurat_obj),
    nFeature_RNA = seurat_obj$nFeature_RNA,
    nCount_RNA = seurat_obj$nCount_RNA,
    percent.mt = seurat_obj$percent.mt
  ) %>%
    group_by(Cluster) %>%
    summarize(
      count = n(),
      ## nFeature_RNA stats
      mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
      median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
      sd_nFeature_RNA = sd(nFeature_RNA, na.rm = TRUE),
      min_nFeature_RNA = min(nFeature_RNA, na.rm = TRUE),
      q1_nFeature_RNA = quantile(nFeature_RNA, 0.25, na.rm = TRUE),
      q3_nFeature_RNA = quantile(nFeature_RNA, 0.75, na.rm = TRUE),
      max_nFeature_RNA = max(nFeature_RNA, na.rm = TRUE),
      ## nCount_RNA stats
      mean_nCount_RNA = mean(nCount_RNA, na.rm = TRUE),
      median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
      sd_nCount_RNA = sd(nCount_RNA, na.rm = TRUE),
      min_nCount_RNA = min(nCount_RNA, na.rm = TRUE),
      q1_nCount_RNA = quantile(nCount_RNA, 0.25, na.rm = TRUE),
      q3_nCount_RNA = quantile(nCount_RNA, 0.75, na.rm = TRUE),
      max_nCount_RNA = max(nCount_RNA, na.rm = TRUE),
      ## percent.mt stats
      mean_percent_mt = mean(percent.mt, na.rm = TRUE),
      median_percent_mt = median(percent.mt, na.rm = TRUE),
      sd_percent_mt = sd(percent.mt, na.rm = TRUE),
      min_percent_mt = min(percent.mt, na.rm = TRUE),
      q1_percent_mt = quantile(percent.mt, 0.25, na.rm = TRUE),
      q3_percent_mt = quantile(percent.mt, 0.75, na.rm = TRUE),
      max_percent_mt = max(percent.mt, na.rm = TRUE)
    ) %>%
    ungroup()
}

# ------------------------------
# Calculate summary statistics for different groupings
# ------------------------------

# Example using "study_simple" as the grouping variable
cell_counts_condition <- calculate_cell_counts(seurat_obj, group_var = "study_simple")
metrics_by_condition <- calculate_metrics_summary(seurat_obj, group_var = "study_simple")

# Example using "orig.ident" as the grouping variable
cell_counts_origident <- calculate_cell_counts(seurat_obj, group_var = "orig.ident")
metrics_by_origident <- calculate_metrics_summary(seurat_obj, group_var = "orig.ident")

# Metrics by cluster are common regardless of the grouping variable
metrics_by_cluster <- calculate_metrics_by_cluster(seurat_obj)

# ------------------------------
# Write outputs into a single Excel workbook with multiple sheets
# ------------------------------

# Create a new workbook
wb <- createWorkbook()

# Function to write data and auto-adjust column widths in each worksheet
write_data_and_set_widths <- function(wb, sheet_name, data) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data)
  setColWidths(wb, sheet_name, cols = 1:ncol(data), widths = "auto")
}

# Write each data frame to its own worksheet
write_data_and_set_widths(wb, "Major Cell Counts (Study)", cell_counts_condition)
write_data_and_set_widths(wb, "Major Metrics by Study", metrics_by_condition)
write_data_and_set_widths(wb, "Major Cell Counts (Batch)", cell_counts_origident)
write_data_and_set_widths(wb, "Major Metrics by Batch", metrics_by_origident)
write_data_and_set_widths(wb, "Major Metrics by Cluster", metrics_by_cluster)

# Define the output file path for the Excel workbook
output_file <- "../metadata/seurat_metadata_export_detailed.xlsx"

# Save the workbook (overwrite if the file already exists)
saveWorkbook(wb, output_file, overwrite = TRUE)

print(paste("Detailed metadata exported to", output_file))

### MINOR-LEVEL QC EXPORT ###



Idents(seurat_obj) <- seurat_obj$Minor
DimPlot(seurat_obj)

ordered_levels <- c(
  # Sensory Neurons - A-fiber
  "A-PEP",
  "A-LTMR.Ntrk2",
  "A-Propr.Pvalb",
  # Sensory Neurons - C-fiber
  "C-LTMR",
  "C-NP.Mrgpra3",
  "C-NP.Mrgprd",
  "C-NP.Sst",
  "C-Thermo.Trpm8",
  "C-Thermo.Rxfp1",
  "C-PEP.Oprk1",
  "C-PEP.Sstr2",
  "C-PEP.Dcn",
  "U.Atf3",
  # Satellite Glial Cells
  "SGC",
  "SGC IEG",
  "SGC IFN",
  # Schwann Cells - Non-myelinating
  "nmSC",
  "nmSC IEG",
  "nmSC Ngfr+",
  # Schwann Cells - Myelinating
  "mSC",
  "mSC IEG",
  "mSC REP",
  # Boundary Cap
  "BC",
  # Endoneurial Fibroblasts
  "EF I",
  "EF II",
  "EF III",
  "EF IV",
  "EF V",
  "EF VI",
  "EF VII",
  "EF VIII",
  "EF IMM",
  # Epineurial/Meningeal
  "EM I",
  "EM II",
  "EM III",
  "EM IV",
  "EM IMM",
  # Other Mesenchymal
  "PM I",
  "PM II",
  "PM III",
  "FSC",
  "MSC",
  "MLC",
  # Vascular Endothelial - Arterial
  "VEC A i",
  "VEC A ii",
  "VEC A iii",
  "VEC A iv",
  "VEC A v",
  # Vascular Endothelial - Venous
  "VEC V i",
  "VEC V ii",
  "VEC V iii",
  "VEC V iv",
  # Vascular Endothelial - Lymphatic
  "VEC LA",
  "VEC LV",
  "LEC",
  # Vascular Endothelial - Other
  "VEC IMM",
  # Pericyte
  "PERI I",
  "PERI II",
  "PERI IMM",
  # Vascular Smooth Muscle / Mural
  "VSMC I",
  "VSMC II",
  "VSMC III",
  "VSC",
  # Macrophages
  "MAC A",
  "MAC B",
  "MAC M1",
  "MAC M2",
  "MAC FIB",
  "epiMAC",
  # Monocytes & Dendritic Cells
  "MONO",
  "MoDC",
  "cDC",
  "pDC",
  "DCx",
  # Other Immune
  "T NK",
  "MAST",
  "EOS",
  # Blood / Progenitors
  "ERY",
  "HSC",
  "PROG I",
  "PROG II",
  # Other / Unassigned
  "NT A",
  "NT B",
  "IMG"
)

Idents(seurat_obj) <- factor(Idents(seurat_obj), levels = ordered_levels)

# ------------------------------
# Calculate summary statistics for different groupings
# ------------------------------

# Example using "study_simple" as the grouping variable
cell_counts_condition <- calculate_cell_counts(seurat_obj, group_var = "study_simple")
metrics_by_condition <- calculate_metrics_summary(seurat_obj, group_var = "study_simple")

# Example using "orig.ident" as the grouping variable
cell_counts_origident <- calculate_cell_counts(seurat_obj, group_var = "orig.ident")
metrics_by_origident <- calculate_metrics_summary(seurat_obj, group_var = "orig.ident")

# Metrics by cluster are common regardless of the grouping variable
metrics_by_cluster <- calculate_metrics_by_cluster(seurat_obj)

# ------------------------------
# Write outputs into a single Excel workbook with multiple sheets
# ------------------------------

# Create a new workbook
wb <- createWorkbook()

# Function to write data and auto-adjust column widths in each worksheet
write_data_and_set_widths <- function(wb, sheet_name, data) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data)
  setColWidths(wb, sheet_name, cols = 1:ncol(data), widths = "auto")
}

# Write each data frame to its own worksheet
write_data_and_set_widths(wb, "Minor Cell Counts (Study)", cell_counts_condition)
write_data_and_set_widths(wb, "Minor Metrics by Study", metrics_by_condition)
write_data_and_set_widths(wb, "Minor Cell Counts (Batch)", cell_counts_origident)
write_data_and_set_widths(wb, "Minor Metrics by Batch", metrics_by_origident)
write_data_and_set_widths(wb, "Minor Metrics by Cluster", metrics_by_cluster)

# Define the output file path for the Excel workbook
output_file <- "../metadata/seurat_metadata_export_detailed_MINOR.xlsx"

# Save the workbook (overwrite if the file already exists)
saveWorkbook(wb, output_file, overwrite = TRUE)

print(paste("Detailed metadata exported to", output_file))

### FINAL MARKERS EXPORT ###

# UPDATED 20260115 with new object, rerun to export data with correct Major/Minor divisions
seurat_obj <- readRDS("../r_objects/DRGSNMI_main_multilevel.RDS")

# Define logical order: Neurons -> Glia -> Fibro -> Vasc -> Other
ordered_levels <- c(
  # Neurons (organized by functional groups)

  "NEU_A",
  "NEU_C",
  "NEU_ATF3",

  # Glia
  "SGC",                                               # Satellite glia
  "NMSC",                                              # Non-myelinating Schwann
  "MYSC",                                              # Myelinating Schwann
  "GSC",                                               # Glial stem

  # Fibroblasts/Connective
  "EF",                                                #
  "EM",                                                #
  "PM",                                                #
  "FSC",                                               #

  # Vascular
  "VEC_A",                                              # Arterial endothelial
  "VEC_V",                                              # Venous endothelial
  "LEC",                                               # Lymphatic endothelial
  "PERI",                                              # Pericytes
  "VSMC",                                              # Vascular smooth muscle
  "MSC",
  "VSC",

  # Other/Immune
  "MAC", "IMG", "MONO",                         # Macrophage lineage
  "DC",                                                # Dendritic cells
  "LEUK",                                              # Leukocytes (general)
  "MAST",                                              # Mast cells
  "EOS",                                               # Eosinophils
  "ERY",                                               # Erythrocytes
  "NT", "HSC"                           # Stem/other
)


Idents(seurat_obj) <- factor(Idents(seurat_obj), levels = ordered_levels)

output_directory <- "../markers/DRGSNMI_main_named_final_markers"
dir.create(output_directory)
# min.diff.pct = 0.25: require 25% difference in detection rate between groups
find_and_save_all_markers(seurat_obj, output_directory, min.pct = 0.1, min.diff.pct = 0.25, summary = T)

# Load required libraries
library(readr)
library(dplyr)
library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Function to add a sheet to the workbook
add_sheet <- function(cell_type, wb) {
  # Convert cell_type to character in case it's not already
  cell_type <- as.character(cell_type)

  # Construct the filename based on the cluster name
  file <- paste0("../markers/DRGSNMI_main_named_final_markers/FindMarkers_", cell_type, ".csv")

  # Check if file exists
  if (!file.exists(file)) {
    warning(paste("File not found:", file))
    return(NULL)
  }

  # Use the cluster name directly as the sheet name
  sheet_name <- cell_type

  # Read the CSV file
  data <- read_csv(file)

  # Rename the first column to "Gene"
  colnames(data)[1] <- "Gene"

  # Add a new worksheet with the sheet name
  addWorksheet(wb, sheet_name)

  # Write the data to the worksheet
  writeData(wb, sheet_name, data)

  # Apply a text format style to the first column to prevent Excel from converting gene names to dates
  cellStyle <- createStyle(numFmt = "TEXT")
  addStyle(wb, sheet_name, style = cellStyle, rows = 1:(nrow(data)+1), cols = 1, gridExpand = TRUE)
}

# Apply the function to each cluster using the unique identities from the Seurat object
lapply(unique(Idents(seurat_obj)), add_sheet, wb = wb)

# Save the workbook
saveWorkbook(wb, "../markers/DRGSNMI_main_named_final_markers/supplemental-file-2.xlsx", overwrite = TRUE)
