#' ======================================================================
#' horste_pnas_2019_import.R
#'
#' Import Horste PNAS 2019 sciatic nerve dataset (10X V2). GEO: GSE137870.
#'
#' Inputs:  GEO count matrices
#' Outputs: ../../r_objects/batch_objects/horste_*.RDS
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
input_dir <- "../../raw_data/horste_pnas_2019/"

# RDS Outputs
batch_obj_dir <- "../../r_objects/batch_objects/"

# GSE142541_mouse_type.csv contains 2 columns with sample metadata
# col 1 = barcode
# col 2 = tissue type
# NOTE lower input threshold (nGene 200 instead of 500) for this study, from V2 but high quality



# Read each file and save Seurat object
# Add full project/sample name to each object
# Add percent.pt metadata to each object

file <- c("GSE142541")

# Generate full sample name
tissue <- "SN"
condition <- "10X_V2_SN"
method <- "10X"
cellnuc <- "CELL"
species <- "MM"
investigator <- "Horste"
id <- str_extract(file, "GSE\\d+")

sample_name <- paste(tissue, condition, method, cellnuc, species, investigator, id, sep = "_")

message("Processing ", file)
message("Sample Name: ", sample_name)

seurat_obj <- Read10X(paste0(input_dir, file))

seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = sample_name, min.features = 200)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^mt-")

# Read metadata and filter for subset "SC"
metadata_path <- paste0(input_dir, "GSE142541_mouse_type.csv")
metadata <- read.csv(metadata_path, header = FALSE, col.names = c("barcode", "subset"))
sc_subset_barcodes <- metadata$barcode[metadata$subset == "SC"]
nodctrl_subset_barcodes <- metadata$barcode[metadata$subset == "NODCTRL"]
nodicam_subset_barcodes <- metadata$barcode[metadata$subset == "NODICAM"]

# Subset Seurat object to contain only "SC" subset cells
sc_subset <- subset(seurat_obj, cells = sc_subset_barcodes)
nodctrl_subset <- subset(seurat_obj, cells = nodctrl_subset_barcodes)
nodicam_subset <- subset(seurat_obj, cells = nodicam_subset_barcodes)

sc_subset
head(sc_subset@meta.data)
nodctrl_subset
head(nodctrl_subset@meta.data)
nodicam_subset
head(nodicam_subset@meta.data)

# add age, study_sex, study, study_simple, tissue, treatment, treatment_simple to each object


sc_subset$age <- "Adult"
nodctrl_subset$age <- "Adult"
nodicam_subset$age <- "Adult"

sc_subset$study <- "horste_pnas_2019"
nodctrl_subset$study <- "horste_pnas_2019"
nodicam_subset$study <- "horste_pnas_2019"

sc_subset$study_sex <- "F"
nodctrl_subset$study_sex <- "F"
nodicam_subset$study_sex <- "F"

sc_subset$study_simple <- "Wolbert, 2019"
nodctrl_subset$study_simple <- "Wolbert, 2019"
nodicam_subset$study_simple <- "Wolbert, 2019"

sc_subset$tissue <- "SN"
nodctrl_subset$tissue <- "SN"
nodicam_subset$tissue <- "SN"

sc_subset$treatment <- "WT"
nodctrl_subset$treatment <- "WT_NODCTRL"
nodicam_subset$treatment <- "NODICAM"

sc_subset$treatment_simple <- "WT"
nodctrl_subset$treatment_simple <- "WT"
nodicam_subset$treatment_simple <- "NODICAM"

sc_subset@meta.data$tech <- "10X_RNA_V2"
nodctrl_subset@meta.data$tech <- "10X_RNA_V2"
nodicam_subset@meta.data$tech <- "10X_RNA_V2"

sc_subset@meta.data$tech_simple <- "10X"
nodctrl_subset@meta.data$tech_simple <- "10X"
nodicam_subset@meta.data$tech_simple <- "10X"

sc_subset@meta.data$cellnuc <- "CELL"
nodctrl_subset@meta.data$cellnuc <- "CELL"
nodicam_subset@meta.data$cellnuc <- "CELL"





# Save Seurat Object
saveRDS(sc_subset, file = paste0(batch_obj_dir,"horste_pnas_2019_WT_seurat.RDS"))
saveRDS(nodctrl_subset, file = paste0(batch_obj_dir,"horste_pnas_2019_WT_NODCTRL_seurat.RDS"))
saveRDS(nodicam_subset, file = paste0(batch_obj_dir,"horste_pnas_2019_NODICAM_seurat.RDS"))

message("horste_pnas_2019 import complete.")

###############################################################################
### SOLO VERIFICATION (not part of pipeline — for independent QC exploration)
### Uncomment and run interactively to examine this dataset in isolation.
###############################################################################

# combined <- sc_subset
# options(future.globals.maxSize = 3e+09)
# combined <- SCTransform(combined)
# combined <- RunPCA(combined, npcs = 50)
# combined <- RunUMAP(combined, reduction = "pca", dims = 1:50)
# combined <- FindNeighbors(combined, reduction = "pca", dims = 1:50)
# combined <- FindClusters(combined, resolution = 0.5)
# DimPlot(combined, reduction = "umap", label = TRUE, raster = FALSE) + NoLegend()
# DefaultAssay(combined) <- "RNA"
# combined <- NormalizeData(combined)
