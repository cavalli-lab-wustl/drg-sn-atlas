#' ======================================================================
#' 07_main_namexfer.R
#'
#' Transfer major/minor cell type labels from subsets back to the main object.
#'
#' Inputs:  ../r_objects/combined_object_init_named.RDS, ../metadata/*_celltypes.rds
#' Outputs: ../r_objects/DRGSNMI_main_multilevel.RDS
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

# ===== LOAD BASE OBJECT =====
merged_named <- readRDS("../r_objects/combined_object_init_named.RDS")
cat("Total cells in object:", ncol(merged_named), "\n\n")

# ===== LOAD ALL CELL TYPES AND CELL IDs =====

cat("Loading MAJOR level files...\n")
major_fibroblasts_labels <- readRDS("../metadata/fibroblasts_major_celltypes.rds")
major_fibroblasts_cells <- readRDS("../metadata/fibroblasts_major_cells.rds")
major_glia_labels <- readRDS("../metadata/glia_major_celltypes.rds")
major_glia_cells <- readRDS("../metadata/glia_major_cells.rds")
major_immune_labels <- readRDS("../metadata/immune_major_celltypes.rds")
major_immune_cells <- readRDS("../metadata/immune_major_cells.rds")
major_neurons_labels <- readRDS("../metadata/neurons_major_celltypes.rds")
major_neurons_cells <- readRDS("../metadata/neurons_major_cells.rds")
major_mural_labels <- readRDS("../metadata/mural_major_celltypes.rds")
major_mural_cells <- readRDS("../metadata/mural_major_cells.rds")
major_vecs_labels <- readRDS("../metadata/vecs_major_celltypes.rds")
major_vecs_cells <- readRDS("../metadata/vecs_major_cells.rds")
major_imoon_labels <- readRDS("../metadata/imoonglia_major_celltypes.rds")
major_imoon_cells <- readRDS("../metadata/imoonglia_major_cells.rds")

cat("Loading MINOR level files...\n")
minor_fibroblasts_labels <- readRDS("../metadata/fibroblasts_minor_celltypes.rds")
minor_fibroblasts_cells <- readRDS("../metadata/fibroblasts_minor_cells.rds")
minor_glia_labels <- readRDS("../metadata/glia_minor_celltypes.rds")
minor_glia_cells <- readRDS("../metadata/glia_minor_cells.rds")
minor_immune_labels <- readRDS("../metadata/immune_minor_celltypes.rds")
minor_immune_cells <- readRDS("../metadata/immune_minor_cells.rds")
minor_neurons_labels <- readRDS("../metadata/neurons_minor_celltypes.rds")
minor_neurons_cells <- readRDS("../metadata/neurons_minor_cells.rds")
minor_mural_labels <- readRDS("../metadata/mural_minor_celltypes.rds")
minor_mural_cells <- readRDS("../metadata/mural_minor_cells.rds")
minor_vecs_labels <- readRDS("../metadata/vecs_minor_celltypes.rds")
minor_vecs_cells <- readRDS("../metadata/vecs_minor_cells.rds")
minor_imoon_labels <- readRDS("../metadata/imoonglia_minor_celltypes.rds")
minor_imoon_cells <- readRDS("../metadata/imoonglia_minor_cells.rds")

# ===== IDENTIFY CELLS NOT IN SUBSETS =====
cat("\nIdentifying cells not in subsets...\n")
all_subset_cells <- c(major_fibroblasts_cells, major_glia_cells, major_immune_cells,
                      major_neurons_cells, major_mural_cells, major_vecs_cells,
                      major_imoon_cells)
cat("Cells in subsets:", length(all_subset_cells), "\n")

remaining_cells <- setdiff(colnames(merged_named), all_subset_cells)
cat("Cells NOT in subsets:", length(remaining_cells), "\n")

# Get original idents for remaining cells
if (length(remaining_cells) > 0) {
  remaining_labels <- as.character(Idents(merged_named)[remaining_cells])
  cat("Unique labels in remaining cells:", unique(remaining_labels), "\n")
}

# ===== CREATE METADATA COLUMN FOR MINOR LEVEL =====
cat("\nCreating Minor annotations...\n")
merged_named$Minor <- as.character(NA)

merged_named$Minor[minor_immune_cells] <- as.character(minor_immune_labels)
merged_named$Minor[minor_vecs_cells] <- as.character(minor_vecs_labels)
merged_named$Minor[minor_fibroblasts_cells] <- as.character(minor_fibroblasts_labels)
merged_named$Minor[minor_glia_cells] <- as.character(minor_glia_labels)
merged_named$Minor[minor_neurons_cells] <- as.character(minor_neurons_labels)
merged_named$Minor[minor_mural_cells] <- as.character(minor_mural_labels)
merged_named$Minor[minor_imoon_cells] <- as.character(minor_imoon_labels)

# Add remaining cells
if (length(remaining_cells) > 0) {
  merged_named$Minor[remaining_cells] <- remaining_labels
}

print(unique(merged_named$Minor))
print(unique(Idents(merged_named)))

# Apply renaming
merged_named$Minor[merged_named$Minor == "neuron"] <- "BATCH"
merged_named$Minor[merged_named$Minor == "ery"] <- "ERY"

cat("Minor level - Unique labels:", length(unique(merged_named$Minor[!is.na(merged_named$Minor)])), "\n")

# ===== CREATE METADATA COLUMN FOR MAJOR LEVEL =====
cat("\nCreating Major annotations...\n")
merged_named$Major <- as.character(NA)

merged_named$Major[major_immune_cells] <- as.character(major_immune_labels)
merged_named$Major[major_vecs_cells] <- as.character(major_vecs_labels)
merged_named$Major[major_fibroblasts_cells] <- as.character(major_fibroblasts_labels)
merged_named$Major[major_glia_cells] <- as.character(major_glia_labels)
merged_named$Major[major_neurons_cells] <- as.character(major_neurons_labels)
merged_named$Major[major_mural_cells] <- as.character(major_mural_labels)
merged_named$Major[major_imoon_cells] <- as.character(major_imoon_labels)

# Add remaining cells
if (length(remaining_cells) > 0) {
  merged_named$Major[remaining_cells] <- remaining_labels
}

# Apply first round of renaming
merged_named$Major[merged_named$Major == "ery"] <- "ERY"
merged_named$Major[merged_named$Major == "nmSC"] <- "NMSC"
merged_named$Major[merged_named$Major == "mSC"] <- "MYSC"
merged_named$Major[merged_named$Major == "VEC A"] <- "VEC_A"
merged_named$Major[merged_named$Major == "VEC V"] <- "VEC_V"
merged_named$Major[merged_named$Major == "C"] <- "NEU_C"
merged_named$Major[merged_named$Major == "A"] <- "NEU_A"
merged_named$Major[merged_named$Major == "Atf3"] <- "NEU_ATF3"

cat("Major level before neuron collapse - Unique labels:", length(unique(merged_named$Major[!is.na(merged_named$Major)])), "\n")

# ===== CREATE METADATA COLUMN FOR CLASS LEVEL =====
cat("\nCreating Class annotations...\n")
merged_named$Class <- merged_named$Major

# Apply MAJOR renaming
major_to_class_remap <- c(
  "PERI" = "MURAL", "VSMC" = "MURAL", "M_CYC" = "MURAL",
  "NEU_A" = "NEU", "NEU_C" = "NEU", "NEU_ATF3" = "NEU",
  "SGC" = "GLIA", "MYSC" = "GLIA", "NMSC" = "GLIA", "G_CYC" = "GLIA",
  "EM" = "FIBRO", "EF" = "FIBRO", "PM" = "FIBRO", "F_CYC" = "FIBRO",
  "VEC_A" = "EC", "VEC_V" = "EC", "LEC" = "EC", "E_CYC" = "EC",
  "IMG" = "IMM", "MAC" = "IMM", "MONO" = "IMM", "DC" = "IMM",
  "NT" = "IMM", "EOS" = "IMM", "MAST" = "IMM", "LEUK" = "IMM", "H_CYC" = "IMM",
  "neuron" = "BATCH", "ERY" = "ERY"
)

for (old_name in names(major_to_class_remap)) {
  merged_named$Class[merged_named$Class == old_name] <- major_to_class_remap[old_name]
}

cat("Class level - Unique labels:", length(unique(merged_named$Class[!is.na(merged_named$Class)])), "\n")

# ===== SET MAJOR AS DEFAULT IDENTS AND REMOVE BATCH/NEURON =====
cat("\nSetting Major as idents and removing BATCH/neuron...\n")
Idents(merged_named) <- merged_named$Major

cells_before <- ncol(merged_named)
cells_to_remove <- c("BATCH", "neuron")
if (any(cells_to_remove %in% Idents(merged_named))) {
  merged_named <- subset(merged_named, idents = cells_to_remove, invert = TRUE)
  cat("BATCH/neuron cells removed:", cells_before - ncol(merged_named), "\n")
}

cat("Major level after neuron collapse - Unique labels:", length(unique(merged_named$Major)), "\n")

# ===== RERUN UMAP =====
cat("\nRerunning UMAP...\n")
merged_named <- RunUMAP(merged_named, reduction = "integrated.dr", dims = 1:50)

# ===== FINAL QC =====
cat("\n===== FINAL QC =====\n")
cat("Final cells:", ncol(merged_named), "\n")
cat("Final Class labels:", length(unique(merged_named$Class)), "\n")
cat("Final Major labels:", length(unique(merged_named$Major)), "\n")
cat("Final Minor labels:", length(unique(merged_named$Minor)), "\n")

cat("Class labels:", sort(unique(merged_named$Class)), "\n")
cat("Major labels:", sort(unique(merged_named$Major)), "\n")
cat("Minor labels:", sort(unique(merged_named$Minor)), "\n")

cat("\nMajor distribution:\n")
print(table(merged_named$Major))

# ===== VISUALIZE =====
DimPlot(merged_named, raster = FALSE, label = TRUE)

# ===== NORMALIZE =====
DefaultAssay(merged_named) <- "RNA"
merged_named <- NormalizeData(merged_named)
merged_named <- JoinLayers(merged_named)

# ===== SAVE OBJECT =====
saveRDS(merged_named, "../r_objects/DRGSNMI_main_multilevel.RDS")
cat("\nSaved to: ../r_objects/DRGSNMI_main_multilevel.RDS\n")

message("07_main_namexfer complete.")
