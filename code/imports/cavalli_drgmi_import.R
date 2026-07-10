#' ======================================================================
#' cavalli_drgmi_import.R
#'
#' Import Cavalli lab DRG/MI dataset (17 samples, 6 conditions). GEO: GSE317728.
#'
#' Inputs:  Raw 10X count matrices
#' Outputs: ../../r_objects/batch_objects/cavalli_*.RDS
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

# Input
input_dir <- "../../raw_data/"

# RDS Outputs
single_obj_dir <- "../../r_objects/single_objects/"
batch_obj_dir <- "../../r_objects/batch_objects/"

# File Names
files <- c(
    "contra_n1_nc20",
    "contra_n2_nc20",
    "snc_n1_nc20",
    "snc_n2_nc20",
    "aih_n1_el21",
    "aih_n2_el21",
    "contra_n3_el21",
    "contra_n4_el21",
    "drc_n1_el21",
    "drc_n2_el21",
    "sci_n1_el21",
    "sci_n2_el21",
    "snc_n3_el21",
    "contra_n5_unpb",
    "naive_n1_unpb",
    "naive_n2_unpb",
    "snc_n4_unpb"
)

simple_conditions <- c(
    "contra",
    "contra",
    "snc",
    "snc",
    "aih",
    "aih",
    "contra",
    "contra",
    "drc",
    "drc",
    "sci",
    "sci",
    "snc",
    "contra",
    "naive",
    "naive",
    "snc"
)

study <- c(
    "cavalli_natcomm_2020",
    "cavalli_natcomm_2020",
    "cavalli_natcomm_2020",
    "cavalli_natcomm_2020",
    "cavalli_elife_2021",
    "cavalli_elife_2021",
    "cavalli_elife_2021",
    "cavalli_elife_2021",
    "cavalli_elife_2021",
    "cavalli_elife_2021",
    "cavalli_elife_2021",
    "cavalli_elife_2021",
    "cavalli_elife_2021",
    "cavalli_unpub",
    "cavalli_unpub",
    "cavalli_unpub",
    "cavalli_unpub"
)


# used to load data in from correct study directory
file_study <- setNames(study, files)

# used to map object to condition for batch combining
file_simple_conditions <- setNames(simple_conditions, files)

# Read each file and save Seurat object
# Add full project/sample name to each object
# Add percent.pt metadata to each object

# Initialize list to store Seurat objects for merging into study-level objects
contra_objects <- list()
snc_objects <- list()
aih_objects <- list()
drc_objects <- list()
sci_objects <- list()
naive_objects <- list()

# Loop through each file, read it, and save it as a Seurat object
for (file in files) {

  # Generate full sample name
  tissue <- "DRG"
  condition <- file
  method <- "10X"
  cellnuc <- "CELL"
  species <- "MM"
  investigator <- "Cavalli"
  sample_name <- paste(tissue, condition, method, cellnuc, species, investigator, sep = "_")

  message("Processing ", file)
  message("Sample Name: ", sample_name)

  seurat_obj <- Read10X(paste0(input_dir, file_study[file], "/", file, "/raw_feature_bc_matrix"))

  seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = sample_name, min.features = 500)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^mt-")

  # Save Seurat Object
  saveRDS(seurat_obj, file = paste0(single_obj_dir, sample_name, ".RDS"))

  file_condition <- file_simple_conditions[file]

  # Add file to appropriate list
    if (file_condition == "contra") {
            contra_objects <- c(contra_objects, seurat_obj)
        } else if (file_condition == "snc") {
            snc_objects <- c(snc_objects, seurat_obj)
        } else if (file_condition == "aih") {
            aih_objects <- c(aih_objects, seurat_obj)
        } else if (file_condition == "drc") {
            drc_objects <- c(drc_objects, seurat_obj)
        } else if (file_condition == "sci") {
            sci_objects <- c(sci_objects, seurat_obj)
        } else if (file_condition == "naive") {
            naive_objects <- c(naive_objects, seurat_obj)
        } else {
            message("ERROR: Unknown file_condition: ", file_condition)
    }

}

# Combine Seurat objects
contra_seurat <- merge(contra_objects[[1]], y = contra_objects[-1])
snc_seurat <- merge(snc_objects[[1]], y = snc_objects[-1])
aih_seurat <- merge(aih_objects[[1]], y = aih_objects[-1])
drc_seurat <- merge(drc_objects[[1]], y = drc_objects[-1])
sci_seurat <- merge(sci_objects[[1]], y = sci_objects[-1])
naive_seurat <- merge(naive_objects[[1]], y = naive_objects[-1])

# Remove individual objects from memory (optional)


# add age, study_sex, study, study_simple, tissue, treatment, treatment_simple to each object

contra_seurat@meta.data$age <- "Adult"
snc_seurat@meta.data$age <- "Adult"
aih_seurat@meta.data$age <- "Adult"
drc_seurat@meta.data$age <- "Adult"
sci_seurat@meta.data$age <- "Adult"
naive_seurat@meta.data$age <- "Adult"


contra_seurat@meta.data$study_sex <- "F"
snc_seurat@meta.data$study_sex <- "F"
aih_seurat@meta.data$study_sex <- "F"
drc_seurat@meta.data$study_sex <- "F"
sci_seurat@meta.data$study_sex <- "F"
naive_seurat@meta.data$study_sex <- "F"


contra_seurat@meta.data$study <- "cavalli_contra"
snc_seurat@meta.data$study <- "cavalli_snc"
aih_seurat@meta.data$study <- "cavalli_aih"
drc_seurat@meta.data$study <- "cavalli_drc"
sci_seurat@meta.data$study <- "cavalli_sci"
naive_seurat@meta.data$study <- "cavalli_naive"

contra_seurat@meta.data$study_simple <- "Avraham, 2020, 2021, & Unpublished"
snc_seurat@meta.data$study_simple <- "Avraham, 2020, 2021, & Unpublished"
aih_seurat@meta.data$study_simple <- "Avraham, 2020, 2021, & Unpublished"
drc_seurat@meta.data$study_simple <- "Avraham, 2020, 2021, & Unpublished"
sci_seurat@meta.data$study_simple <- "Avraham, 2020, 2021, & Unpublished"
naive_seurat@meta.data$study_simple <- "Avraham, 2020, 2021, & Unpublished"

contra_seurat@meta.data$tissue <- "DRG"
snc_seurat@meta.data$tissue <- "DRG"
aih_seurat@meta.data$tissue <- "DRG"
drc_seurat@meta.data$tissue <- "DRG"
sci_seurat@meta.data$tissue <- "DRG"
naive_seurat@meta.data$tissue <- "DRG"

contra_seurat@meta.data$treatment <- "WT_CONTRA"
snc_seurat@meta.data$treatment <- "SNC_3DPI"
aih_seurat@meta.data$treatment <- "AIH_3DPI"
drc_seurat@meta.data$treatment <- "DRC_3DPI"
sci_seurat@meta.data$treatment <- "SCI_3DPI"
naive_seurat@meta.data$treatment <- "WT_NAIVE"

contra_seurat@meta.data$treatment_simple <- "WT"
snc_seurat@meta.data$treatment_simple <- "SNC"
aih_seurat@meta.data$treatment_simple <- "AIH"
drc_seurat@meta.data$treatment_simple <- "DRC"
sci_seurat@meta.data$treatment_simple <- "SCI"
naive_seurat@meta.data$treatment_simple <- "WT"

contra_seurat@meta.data$tech <- "10X_RNA"
snc_seurat@meta.data$tech <- "10X_RNA"
aih_seurat@meta.data$tech <- "10X_RNA"
drc_seurat@meta.data$tech <- "10X_RNA"
sci_seurat@meta.data$tech <- "10X_RNA"
naive_seurat@meta.data$tech <- "10X_RNA"

contra_seurat@meta.data$tech_simple <- "10X"
snc_seurat@meta.data$tech_simple <- "10X"
aih_seurat@meta.data$tech_simple <- "10X"
drc_seurat@meta.data$tech_simple <- "10X"
sci_seurat@meta.data$tech_simple <- "10X"
naive_seurat@meta.data$tech_simple <- "10X"

contra_seurat@meta.data$cellnuc <- "CELL"
snc_seurat@meta.data$cellnuc <- "CELL"
aih_seurat@meta.data$cellnuc <- "CELL"
drc_seurat@meta.data$cellnuc <- "CELL"
sci_seurat@meta.data$cellnuc <- "CELL"
naive_seurat@meta.data$cellnuc <- "CELL"


# Save merged  Seurat objects
saveRDS(contra_seurat, file = paste0(batch_obj_dir, "cavalli_contra_seurat.RDS"))
saveRDS(snc_seurat, file = paste0(batch_obj_dir, "cavalli_snc_seurat.RDS"))
saveRDS(aih_seurat, file = paste0(batch_obj_dir, "cavalli_aih_seurat.RDS"))
saveRDS(drc_seurat, file = paste0(batch_obj_dir, "cavalli_drc_seurat.RDS"))
saveRDS(sci_seurat, file = paste0(batch_obj_dir, "cavalli_sci_seurat.RDS"))
saveRDS(naive_seurat, file = paste0(batch_obj_dir, "cavalli_naive_seurat.RDS"))

# Notify Completion
message("All files processed and combined Seurat object saved.")
