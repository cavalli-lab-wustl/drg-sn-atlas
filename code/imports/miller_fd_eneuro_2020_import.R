#' ======================================================================
#' miller_fd_eneuro_2020_import.R
#'
#' Import Miller eNeuro 2020 sciatic nerve dataset (adult+P1). GEO: GSE139103.
#'
#' Inputs:  GEO count matrices
#' Outputs: ../../r_objects/batch_objects/miller_fd_*.RDS
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
input_dir <- "../../raw_data/miller_fd_eneuro_2020/GSE147285_RAW/"

# RDS Outputs
batch_obj_dir <- "../../r_objects/batch_objects/"


# File Names
files <- c(
    "GSM4423506_Inj_Sciatic_3d.txt.gz",
    "GSM4423507_Neo_Sciatic_FACS.txt.gz",
    "GSM4423508_Neo_Sciatic_Beads.txt.gz",
    "GSM4423509_Uninj_Sciatic.txt.gz"
)

simple_conditions <- c(
    "INJ_SN_3d",
    "JUV_SN_FACS",
    "JUV_SN_BEAD",
    "CTRL_SN"
)

file_simple_conditions <- setNames(simple_conditions, files)

# Read each file and save Seurat object
# Add full project/sample name to each object
# Add percent.pt metadata to each object

file <- "GSM4423506_Inj_Sciatic_3d.txt.gz"
# Generate full sample name
tissue <- "SN"
condition <- file_simple_conditions[file]
method <- "DS"
cellnuc <- "CELL"
species <- "MM"
investigator <- "Miller_FD"
id <- str_extract(file, "GSM\\d+")
sample_name <- paste(tissue, condition, method, cellnuc, species, investigator, id, sep = "_")

message("Processing ", file)
message("Sample Name: ", sample_name)

count_data_1 <- read.delim(gzfile(paste0(input_dir, file)), header = TRUE, row.names = 1, check.names = FALSE)

seurat_obj_1 <- CreateSeuratObject(counts = count_data_1, project = sample_name)
seurat_obj_1[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj_1, pattern = "^mt-")


seurat_obj_1@meta.data$age <- "Adult"
seurat_obj_1@meta.data$study_sex <- "Mixed"
seurat_obj_1@meta.data$study <- "miller_fd_eneuro_2020"
seurat_obj_1@meta.data$study_simple <- "Toma, 2020"
seurat_obj_1@meta.data$tissue <- "SN"
seurat_obj_1@meta.data$treatment <- "SNC_3DPI"
seurat_obj_1@meta.data$treatment_simple <- "SNC"

seurat_obj_1@meta.data$tech <- "DS_1_2"
seurat_obj_1@meta.data$tech_simple <- "DS"
seurat_obj_1@meta.data$cellnuc <- "CELL"


seurat_obj_1
head(seurat_obj_1@meta.data)
unique(seurat_obj_1@meta.data$orig.ident)


# Save Seurat Object
saveRDS(seurat_obj_1, file = paste0(batch_obj_dir, sample_name, ".RDS"))


file <- "GSM4423507_Neo_Sciatic_FACS.txt.gz"
# Generate full sample name
tissue <- "SN"
condition <- file_simple_conditions[file]
method <- "DS"
cellnuc <- "CELL"
species <- "MM"
investigator <- "Miller_FD"
id <- str_extract(file, "GSM\\d+")
sample_name <- paste(tissue, condition, method, cellnuc, species, investigator, id, sep = "_")

message("Processing ", file)
message("Sample Name: ", sample_name)

count_data_2 <- read.delim(gzfile(paste0(input_dir, file)), header = TRUE, row.names = 1, check.names = FALSE)

seurat_obj_2 <- CreateSeuratObject(counts = count_data_2, project = sample_name)
seurat_obj_2[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj_2, pattern = "^mt-")


seurat_obj_2@meta.data$age <- "P1"
seurat_obj_2@meta.data$study_sex <- "Mixed"
seurat_obj_2@meta.data$study <- "miller_fd_eneuro_2020"
seurat_obj_2@meta.data$study_simple <- "Toma, 2020"
seurat_obj_2@meta.data$tissue <- "SN"
seurat_obj_2@meta.data$treatment <- "WT_CD1"
seurat_obj_2@meta.data$treatment_simple <- "WT"

seurat_obj_2@meta.data$tech <- "DS_1_2"
seurat_obj_2@meta.data$tech_simple <- "DS"
seurat_obj_2@meta.data$cellnuc <- "CELL"


seurat_obj_2
head(seurat_obj_2@meta.data)
unique(seurat_obj_2@meta.data$orig.ident)


# Save Seurat Object
saveRDS(seurat_obj_2, file = paste0(batch_obj_dir, sample_name, ".RDS"))


file <- "GSM4423508_Neo_Sciatic_Beads.txt.gz"
# Generate full sample name
tissue <- "SN"
condition <- file_simple_conditions[file]
method <- "DS"
cellnuc <- "CELL"
species <- "MM"
investigator <- "Miller_FD"
id <- str_extract(file, "GSM\\d+")
sample_name <- paste(tissue, condition, method, cellnuc, species, investigator, id, sep = "_")

message("Processing ", file)
message("Sample Name: ", sample_name)

count_data_3 <- read.delim(gzfile(paste0(input_dir, file)), header = TRUE, row.names = 1, check.names = FALSE)

seurat_obj_3 <- CreateSeuratObject(counts = count_data_3, project = sample_name)
seurat_obj_3[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj_3, pattern = "^mt-")


seurat_obj_3@meta.data$age <- "P1"
seurat_obj_3@meta.data$study_sex <- "Mixed"
seurat_obj_3@meta.data$study <- "miller_fd_eneuro_2020"
seurat_obj_3@meta.data$study_simple <- "Toma, 2020"
seurat_obj_3@meta.data$tissue <- "SN"
seurat_obj_3@meta.data$treatment <- "WT_CD1"
seurat_obj_3@meta.data$treatment_simple <- "WT"

seurat_obj_3@meta.data$tech <- "DS_1_2"
seurat_obj_3@meta.data$tech_simple <- "DS"
seurat_obj_3@meta.data$cellnuc <- "CELL"


seurat_obj_3
head(seurat_obj_3@meta.data)
unique(seurat_obj_3@meta.data$orig.ident)


# Save Seurat Object
saveRDS(seurat_obj_3, file = paste0(batch_obj_dir, sample_name, ".RDS"))


file <- "GSM4423509_Uninj_Sciatic.txt.gz"
# Generate full sample name
tissue <- "SN"
condition <- file_simple_conditions[file]
method <- "DS"
cellnuc <- "CELL"
species <- "MM"
investigator <- "Miller_FD"
id <- str_extract(file, "GSM\\d+")
sample_name <- paste(tissue, condition, method, cellnuc, species, investigator, id, sep = "_")

message("Processing ", file)
message("Sample Name: ", sample_name)

count_data_4 <- read.delim(gzfile(paste0(input_dir, file)), header = TRUE, row.names = 1, check.names = FALSE)

seurat_obj_4 <- CreateSeuratObject(counts = count_data_4, project = sample_name)
seurat_obj_4[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj_4, pattern = "^mt-")


seurat_obj_4@meta.data$age <- "Adult"
seurat_obj_4@meta.data$study_sex <- "Mixed"
seurat_obj_4@meta.data$study <- "miller_fd_eneuro_2020"
seurat_obj_4@meta.data$study_simple <- "Toma, 2020"
seurat_obj_4@meta.data$tissue <- "SN"
seurat_obj_4@meta.data$treatment <- "WT_CD1"
seurat_obj_4@meta.data$treatment_simple <- "WT"

seurat_obj_4@meta.data$tech <- "DS_1_2"
seurat_obj_4@meta.data$tech_simple <- "DS"
seurat_obj_4@meta.data$cellnuc <- "CELL"


seurat_obj_4
head(seurat_obj_4@meta.data)
unique(seurat_obj_4@meta.data$orig.ident)


# Save Seurat Object
saveRDS(seurat_obj_4, file = paste0(batch_obj_dir, sample_name, ".RDS"))


message("miller_fd_eneuro_2020 import complete.")

###############################################################################
### SOLO VERIFICATION (not part of pipeline — for independent QC exploration)
### Uncomment and run interactively to examine this dataset in isolation.
### Note: uses RPCA integration across the 4 samples (different ages/treatments).
###############################################################################

# combined <- merge(seurat_obj_1, y = list(seurat_obj_2, seurat_obj_3, seurat_obj_4))
# options(future.globals.maxSize = 3e+09)
# combined <- SCTransform(combined)
# combined <- RunPCA(combined, npcs = 50)
# combined <- IntegrateLayers(combined, method = RPCAIntegration, normalization.method = "SCT")
# combined <- FindNeighbors(combined, dims = 1:50, reduction = "integrated.dr")
# combined <- FindClusters(combined, resolution = 0.5)
# combined <- RunUMAP(combined, reduction = "integrated.dr", dims = 1:50)
# DimPlot(combined, reduction = "umap", label = TRUE, raster = FALSE) + NoLegend()
# DefaultAssay(combined) <- "RNA"
# combined <- NormalizeData(combined)
